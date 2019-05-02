#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Class for managing first arrival travel time inversions"""
import numpy as np

from matplotlib.collections import LineCollection

import pygimli as pg
from pygimli.frameworks import MeshModelling
from pygimli.manager import MeshMethodManager

from . raplot import drawTravelTimeData, drawVA, showVA
from . ratools import shotReceiverDistances


class TravelTimeDijkstraModelling(MeshModelling):
    def __init__(self, **kwargs):
        self.dijkstra = pg.TravelTimeDijkstraModelling()

        super(TravelTimeDijkstraModelling, self).__init__(**kwargs)

        self._refineSecNodes = 3
        self.jacobian = self.dijkstra.jacobian
        self.createJacobian = self.dijkstra.createJacobian

        self.setJacobian(self.dijkstra.jacobian())

    def regionManagerRef(self):
        # necessary because core dijkstra use its own RM
        return self.dijkstra.regionManagerRef()

    def refineFwdMesh(self):
        """Refine the current mesh for higher accuracy.

        This is called automatic when accesing self.mesh() so it ensures any
        effect of changing region properties (background, single).
        """
        pg.info("Creating refined mesh (secnodes: {0}) to "
                "solve forward task.".format(self._refineSecNodes))
        self._mesh = self._mesh.createMeshWithSecondaryNodes(self._refineSecNodes)
        pg.verbose(self._mesh)

    def setMeshPost(self, mesh):
        """
        """
        pg._r(mesh)
        self.dijkstra.setMesh(mesh)
        #self.dijkstra.setMesh(pg.Mesh(mesh))

    def setDataPost(self, data):
        """
        """
        # pg._r()
        self.dijkstra.setData(data)

    def createDefaultStartModel(self, dataVals):
        """
        """
        dists = shotReceiverDistances(self.data, full=True)

        aSlow = 1. / (dists / dataVals)

        sm = pg.Vector(self.regionManager().parameterCount(),
                       pg.median(aSlow))
        return sm

    def response(self, par):
        return self.dijkstra.response(par)

    def drawModel(self, ax, model, **kwargs):
        kwargs['label'] = pg.unit('vel')
        super(TravelTimeDijkstraModelling, self).drawModel(ax=ax,
                                                           model=model,
                                                           **kwargs)
        return ax

    def drawData(self, ax, data, err=None, **kwargs):
        kwargs['label'] = pg.unit('va')
        return showVA(self.data, vals=data, usePos=False,
                      ax=ax, **kwargs)


class TravelTimeManager(MeshMethodManager):
    """Manager for refraction seismics (traveltime tomography)

    TODO Document main members and use default MethodManager interface
    e.g., self.inv, self.fop, self.paraDomain, self.mesh, self.data
    """

    def __init__(self, **kwargs):
        """Init function with optional data load"""
        self._useFMM = False

        super(TravelTimeManager, self).__init__(**kwargs)

        self._dataToken = 't'
        self.inv.dataTrans = pg.RTransLog()

    def errorValues(self, data):
        """Return relative error values from a given DataContainer."""
        if not data.haveData('err'):
            pg.error('Datacontainer have no "err" values. Fallback set to 0.01')

        return pg.Vector(data('err') / data('t'))

    def setMesh(self, mesh, secNodes=0):
        """ """
        if secNodes > 0:
            self.fop._refineSecNodes = secNodes

        self.fop.setMesh(mesh)

    def createForwardOperator(self, **kwargs):
        """Create default forward operator for Traveltime modelling.

        Your want your Manager use a special forward operator you can add them
        here on default Dijkstra is used.
        """
        fop = TravelTimeDijkstraModelling(**kwargs)

        return fop

    def simulate(self, mesh, slowness, scheme, secNodes=2, **kwargs):
        """
        Simulate an Traveltime measurement.

        Perform the forward task for a given mesh,
        a slowness distribution (per cell) and return data
        (Traveltime) for a measurement scheme.
        This is a static method since it does not interfere with the Managers
        inversion approaches.

        Parameters
        ----------
        mesh : :gimliapi:`GIMLI::Mesh`
            Mesh to calculate for.

        slowness : array(mesh.cellCount()) | array(N, mesh.cellCount())
            slowness distribution for the given mesh cells can be:

            * a single array of len mesh.cellCount()
            * a matrix of N slowness distributions of len mesh.cellCount()
            * a res map as [[marker0, res0], [marker1, res1], ...]

        scheme : :gimliapi:`GIMLI::DataContainer`
            data measurement scheme

        Other Parameters
        ----------------
        noisify : add normal distributed noise based on scheme('err')
            IMPLEMENTME

        Returns
        -------
        t : array(N, data.size()) | DataContainer
            The resulting simulated travel time values.
            Either one column array or matrix in case of slowness matrix.
            A DataContainer is return if noisify set to True.

        """
        fop = self.createForwardOperator()

        fop.setData(scheme)
        self.setMesh(mesh, secNodes=secNodes)

        if len(slowness) == mesh.cellCount():
            if max(slowness) > 1.:
                print('Warning: slowness values larger than 1 (' +
                      str(max(slowness)) + ').. assuming that are velocity '
                      'values .. building reciprocity')
                t = fop.response(1./slowness)
            else:
                t = fop.response(slowness)
        else:
            print(mesh)
            print("slowness: ", slowness)
            raise BaseException("Simulate called with wrong slowness array.")

        ret = pg.DataContainer(scheme)
        ret.set('t', t)

        noiseLevel = kwargs.pop('noiseLevel', 0)

        if noiseLevel > 0:
            if not ret.allNonZero('err'):
                ret.set('t', t)
                ret.set('err', pg.physics.Refraction.estimateError(
                        ret, absoluteError=kwargs.pop('noiseAbs', 1e-4),
                        relativeError=noiseLevel))

            if self.verbose:
                print("Data error estimates (min:mpythonax) ",
                      min(ret('err')), ":", max(ret('err')))

            t += pg.randn(ret.size()) * ret('err')
            ret.set('t', t)

        if kwargs.pop('returnArray', False):
            return t

        return ret

    def invert(self, data=None, **kwargs):
        """Invert data.

        Parameters
        ----------
        data : pg.DataContainer()
            Data container with at least SensorIndieces 's g' and
            data values 't' (traveltime in ms) and 'err' (absolute error in ms)

        Other Parameters
        ----------------
        secNodes: int [2]
            Amount of secondary nodes used for ensure accuracy of the forward
            operator.
        """
        kwargs.setdefault('secNodes', 2)

        self.fop._refineSecNodes = kwargs.pop('secNodes', 2)

        if 'limits' in kwargs:
            if kwargs['limits'][0] > 1:
                tmp = kwargs['limits'][0]
                kwargs['limits'][0] = 1.0 / kwargs['limits'][1]
                kwargs['limits'][1] = 1.0 / tmp

        slowness = super(TravelTimeManager, self).invert(data, **kwargs)
        velocity = 1.0 / slowness
        self.fw.model = velocity
        return velocity

    def drawRayPaths(self, ax, model=None, **kwargs):
        """Draw the the ray paths for `model` or last model for
        which the last Jacobian was calculated.

        Parameters
        ----------
        model : array
            Velocity model for which to calculate and visualize ray paths (the
            default is model for last Jacobian calculation in self.velocity).
        ax : matplotlib.axes object
            To draw the model and the path into.
        **kwargs : type
            Additional arguments passed to LineCollection (alpha, linewidths,
            color, linestyles).

        Returns
        -------
        lc : matplotlib.LineCollection
        """
        if model is not None:
            self.fop.createJacobian(1/model)
        else:
            model = self.model

        _ = kwargs.setdefault("color", "w")
        _ = kwargs.setdefault("alpha", 0.5)
        _ = kwargs.setdefault("linewidths", 0.8)

        shots = self.fop.data.id("s")
        recei = self.fop.data.id("g")

        segs = []
        for s, g in zip(shots, recei):
            wi = self.fop.dijkstra.way(s, g)
            points = self.fop.dijkstra.mesh().positions(withSecNodes=True)[wi]
            segs.append(np.column_stack((pg.x(points), pg.y(points))))

        lc = LineCollection(segs, **kwargs)
        ax.add_collection(lc)

        return lc

    def showRayPaths(self, model=None, ax=None, **kwargs):
        """Show the model with ray paths for `model` or last model for
        which the last Jacobian was calculated.

        Parameters
        ----------
        model : array
            Velocity model for which to calculate and visualize ray paths (the
            default is model for last Jacobian calculation in self.velocity).
        ax : matplotlib.axes object
            To draw the model and the path into.
        **kwargs : type
            forward to drawRayPaths

        Returns
        -------
        ax : matplotlib.axes object
        cb : matplotlib.colorbar object (only if model is provided)

        Examples
        --------
        >>> # No reason to import matplotlib
        >>> import pygimli as pg
        >>> from pygimli.physics import TravelTimeManager
        >>> from pygimli.physics.traveltime import createRAData
        >>>
        >>> x, y = 8, 6
        >>> mesh = pg.createGrid(x, y)
        >>> data = createRAData([(0,0)] + [(x, i) for i in range(y)],
        ...                     shotdistance=y+1)
        >>> data.set("t", pg.RVector(data.size(), 1.0))
        >>> tt = TravelTimeManager()
        >>> tt.fop.setData(data)
        >>> tt.setMesh(mesh, secNodes=10)
        >>> ax, cb = tt.showRayPaths(showMesh=True, diam=0.1)
        """
        if model is None:
            if self.fop.jacobian().size() == 0:
                self.fop.mesh() # initialize any meshs .. just to be sure is 1
                model = pg.Vector(self.fop.regionManager().parameterCount(),
                                  1.0)
            else:
                model = self.model

        ax, cbar = self.showModel(ax=ax, model=model,
                                  showMesh=kwargs.pop('showMesh', None),
                                  diam=kwargs.pop('diam', None))
        self.drawRayPaths(ax, model=model, **kwargs)

        return ax, cbar


if __name__ == '__main__':
    pg.wait()
