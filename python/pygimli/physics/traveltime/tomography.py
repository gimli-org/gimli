#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tomography Method Manager"""

import sys
import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.mplviewer.dataview import plotVecMatrix
from . refraction import Refraction
from . importData import readTOMfile


class Tomography(Refraction):

    """Traveltime tomography for tomographic (e.g. crosshole) measurements"""

    def __init__(self, data=None, tcorr=0, name='new', **kwargs):
        """Init function with optional data load

        Parameters
        ----------
        data : pg.DataContainer or string

        tcorr : float [0]
            correct travel times by common shift

        name : str [data if being string, otherwise 'new']
            basename for saving Figures, results etc.

        ndig : int [2]
            number of digits to round positions (e.g. 2=cm), alternatively:

        roundto : float [0]
            unit spacing to round positions on
        """
        if isinstance(data, str):
            name = data[:data.rfind('.')]
            if data.lower()[-4:] == '.tom':
                data = readTOMfile(data, **kwargs)
            else:
                data = pg.DataContainer(data, 's g')
        if tcorr != 0:
            data.set('t', data('t') + tcorr)
        super(Tomography, self).__init__(data, name=name, **kwargs)

    def createMesh(self, quality=34.6, maxarea=0.1, addpoints=None):
        """Create (inversion) mesh by circumventing PLC"""
        data = self.dataContainer
        sx = list(pg.x(data.sensorPositions()))
        sz = list(pg.y(data.sensorPositions()))

        if addpoints is not None:
            for po in addpoints:
                sx.append(po[0])
                sz.append(po[1])

        iS = np.argsort(np.arctan2(sx-np.mean(sx), sz-np.mean(sz)))
        plc = pg.Mesh(2)
        nodes = [plc.createNode(sx[i], sz[i], 0) for i in iS]
        for i in range(len(nodes)-1):
            plc.createEdge(nodes[i], nodes[i+1])

        plc.createEdge(nodes[-1], nodes[0])
        tri = pg.TriangleWrapper(plc)
        tri.setSwitches("-pzFq"+str(quality)+"a"+str(maxarea))
        self.setMesh(tri.generate())

    def offset(self):
        """return shot-geophone distance"""
        data = self.dataContainer
        return np.array([data.sensorPosition(int(data('g')[i])).distance(
          data.sensorPosition(int(data('s')[i]))) for i in range(data.size())])

    def getVA(self, t=None):
        """return apparent velocity"""
        if t is None:
            t = self.dataContainer('t')
        return self.offset() / t

    def createStartModel(self):
        """create (gradient) starting model with vtop/vbottom bounds"""
        va = self.getVA()
        nModel = self.fop.regionManager().parameterCount()
        return pg.RVector(nModel, 1./np.mean(va))

    def showVA(self, t=None, ax=None, usepos=True, name='va', squeeze=True):
        """show apparent velocity as image plot"""
        # va = self.getVA(vals=vals)
        xvec = self.dataContainer('g')
        yvec = self.dataContainer('s')
        if usepos:
            pz = pg.y(self.dataContainer.sensorPositions())
            if squeeze:
                xvec = pz[xvec]
                yvec = pz[yvec]
            else:
                pz = pg.y(self.dataContainer.sensorPositions())
                raise Exception('Implement ME')
                # xvec = px[xvec]*1000 + pz[xvec]
                # xvec = px[yvec]*1000 + pz[yvec]

        plotVecMatrix(xvec, yvec, vals=t, squeeze=squeeze, ax=ax, name=name)

    def showVAold(self, vals=None, ax=None, usepos=True, name='va'):
        """show apparent velocity as image plot (old style)"""
        va = self.getVA(t=vals)
        data = self.dataContainer
        A = np.ones((data.sensorCount(), data.sensorCount())) * np.nan
        for i in range(data.size()):
            A[int(data('s')[i]), int(data('g')[i])] = va[i]

        if ax is None:
            fig, ax = plt.subplots()
            self.figs[name] = fig

        gci = ax.imshow(A, interpolation='nearest')
        ax.grid(True)
        if usepos:
            xt = np.linspace(0, data.sensorCount()-1, 7)
            xt.round()
            px = pg.abs(pg.y(self.dataContainer.sensorPositions()))
            ax.set_xticks(xt)
            ax.set_xticklabels([str(int(px[int(xti)])) for xti in xt])
            ax.set_yticks(xt)
            ax.set_yticklabels([str(int(px[int(xti)])) for xti in xt])

        plt.colorbar(gci, ax=ax)
        return va


if __name__ == '__main__':
    datafile = sys.argv[1]
    tomo = Tomography(datafile)
    print(tomo)
    tomo.createMesh(addpoints=[[6, -22], [0, -8]])
    tomo.showVA()
    tomo.estimateError(absoluteError=0.7, relativeError=0.001)
    tomo.run()
    tomo.showResult()
