#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Method Manager for Magnetics."""
import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
import pygimli.meshtools as mt
from pygimli.viewer import pv
# from pygimli.frameworks import MeshMethodManager
from .MagneticsModelling import MagneticsModelling
from .tools import depthWeighting


class MagManager():  # MeshMethodManager):
    """Magnetics Manager."""

    def __init__(self, data=None, **kwargs):
        """Create Magnetics Manager instance."""
        self.DATA = kwargs.pop("DATA", None)
        self.x = kwargs.pop("x", None)
        self.y = kwargs.pop("y", None)
        self.z = kwargs.pop("z", None)
        self.igrf = kwargs.pop("igrf", None)
        self.mesh = kwargs.pop("mesh", None)
        self.cmp = kwargs.pop("cmp", ["TFA"])
        # super().__init__()
        self.inv = None #
        self.fwd = None
        if isinstance(data, str):
            self.DATA = np.genfromtxt("fault-synth.mag", names=True)
            self.x = self.DATA["x"]
            self.y = self.DATA["y"]
            self.z = np.abs(self.DATA["z"])
            self.cmp = [t for t in self.DATA.dtype.names
                        if t.startswith("B") or t.startswith("T")]

    def showData(self, cmp=None):
        """Show data."""
        cmp = cmp or self.cmp
        nc = 2 if len(cmp) > 1 else 1
        nr = (len(cmp)+1) // 2
        fig, ax = plt.subplots(nr, nc, sharex=True, sharey=True, squeeze=False,
                               figsize=(7, len(self.cmp)*1.2))
        axs = np.atleast_1d(ax.flat)
        for i, c in enumerate(cmp):
            fld = self.DATA[c]
            vv = max(-np.min(fld)*1., np.max(fld)*1.)
            sc = axs[i].scatter(self.x, self.y, c=fld, cmap="bwr", vmin=-vv, vmax=vv)
            axs[i].set_title(c)
            axs[i].set_aspect(1.0)
            fig.colorbar(sc, ax=ax.flat[i])

    def createGrid(self, dx=50, depth=800, bnd=0):
        """Create a grid."""
        x = np.arange(min(self.x)-bnd, max(self.x)+bnd+.1, dx)
        y = np.arange(min(self.y)-bnd, max(self.y)+bnd+.1, dx)
        z = np.arange(-depth, .1, dx)
        self.mesh = mt.createGrid(x=x, y=y, z=z)

    def createMesh(self, area=1e5, quality=1.3, addPLC=None):
        """Create an unstructured mesh."""
        geo = mt.createCube(start=[-600, -600, -800], end=[600, 600, 0])
        for xi, yi in zip(self.x, self.y):
            geo.createNode([xi, yi, 0])
        if addPLC:
            geo += addPLC
        self.mesh = mt.createMesh(geo, quality=quality, area=area)

    def createFOP(self):
        """Create forward operator (computationally extensive!)."""
        points = np.column_stack([self.x, self.y, -np.abs(self.z)])
        self.fwd = MagneticsModelling(self.mesh.NED(), points,
                                      cmp=self.cmp, igrf=self.igrf)

    def inversion(self, noise_level=2, noisify=False, **kwargs):
        """Inversion."""
        datavec = np.concatenate([self.DATA[c] for c in self.cmp])
        if noisify:
            datavec += np.random.randn(len(datavec)) * noise_level

        self.inv = pg.Inversion(fop=self.fwd, verbose=True)
        # self.inv.setForwardOperator(self.fwd)
        kwargs.setdefault("startModel", 0.001)
        kwargs.setdefault("relativeError", 0.001)
        kwargs.setdefault("lam", 10)
        kwargs.setdefault("verbose", True)
        thrs = kwargs.pop("symlogThreshold", 0)
        if thrs > 0:
            self.inv.dataTrans = pg.trans.TransSymLog(thrs)

        limits = kwargs.pop("limits", [0, 0.1])
        self.inv.setRegularization(limits=limits)
        C = kwargs.pop("C", 1)
        cType = kwargs.pop("cType", C)
        if hasattr(C, "__iter__"):
            self.inv.setRegularization(correlationLengths=C)
            cType = -1
        elif isinstance(C, pg.core.MatrixBase):
            self.inv.setRegularization(C=C)
        else:
            self.inv.setRegularization(cType=C)

        z0 = kwargs.pop("z0", 25)  # Oldenburg&Li(1996)
        if kwargs.pop("depthWeighting", False):
            cw = self.fwd.regionManager().constraintWeights()
            dw = depthWeighting(self.mesh, cell=not(cType==1), z0=z0)
            if len(dw) == len(cw):
                dw *= cw

            dw *= kwargs.pop("mul", 1)
            self.inv.setConstraintWeights(dw)

        model = self.inv.run(datavec, absoluteError=noise_level, **kwargs)
        return model

    def showDataFit(self):
        """Show data, model response and misfit."""
        nc = len(self.cmp)
        _, ax = pg.plt.subplots(ncols=3, nrows=nc, figsize=(12, 3*nc), sharex=True, sharey=True, squeeze=False)
        vals = np.reshape(self.inv.dataVals, [nc, -1])
        mm = np.max(np.abs(vals))
        resp = np.reshape(self.inv.response, [nc, -1])
        errs = np.reshape(self.inv.errorVals, [nc, -1])  # relative!
        misf = (vals - resp) / np.abs(errs *  vals)
        fkw = dict(cmap="bwr", vmin=-mm, vmax=mm)
        mkw = dict(cmap="bwr", vmin=-3, vmax=3)
        for i in range(nc):
            ax[i, 0].scatter(self.x, self.y, c=vals[i], **fkw)
            ax[i, 1].scatter(self.x, self.y, c=resp[i], **fkw)
            ax[i, 2].scatter(self.x, self.y, c=misf[i], **mkw)

    def show3DModel(self, label=None, trsh=0.025, synth=None, cMin=0, cMax=0.03, cMap="Spectral_r",
                position="yz", elevation=10, azimuth=25, zoom=1.2, invert=False):
        """Standard 3D view."""
        if label is None:
            label = self.inv.model
        if not isinstance(label, str):
            self.mesh["bla"] = np.array(label)
            label = "bla"

        flt = None
        pl, _ = pg.show(self.mesh, style="wireframe", hold=True,
                        alpha=0.1, backend="trame")
        # mm = [min(self.mesh[label]), min(self.mesh[label])]
        if trsh > 0:
            flt = {"threshold": dict(value=trsh, scalars=label, invert=invert)}
            pv.drawModel(pl, self.mesh, label=label, style="surface",
                        cMin=cMin, cMax=cMax, cMap=cMap, filter=flt)

        pv.drawMesh(pl, self.mesh, label=label, style="surface",
                    cMap=cMap, cMin=cMin, cMax=cMax,
                    filter={"slice": dict(normal=[-1, 0, 0], origin=[0, 0, 0])})

        if synth:
            pv.drawModel(pl, synth, style="wireframe")

        pl.camera_position = position
        pl.camera.azimuth = azimuth
        pl.camera.elevation = elevation
        pl.camera.zoom(zoom)
        pl.show()
        return pl


if __name__ == "__main__":
    pass
