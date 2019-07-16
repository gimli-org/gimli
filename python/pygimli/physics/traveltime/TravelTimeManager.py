#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Class for managing first arrival travel time inversions"""
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.collections import LineCollection

import pygimli as pg
from pygimli.frameworks import MeshModelling
from pygimli.frameworks import MeshMethodManager

from . raplot import drawTravelTimeData, drawVA, showVA
from . ratools import shotReceiverDistances


class TravelTimeDijkstraModelling(MeshModelling):
    def __init__(self, **kwargs):
        self.dijkstra = pg.core.TravelTimeDijkstraModelling()

        super(TravelTimeDijkstraModelling, self).__init__(**kwargs)

        self._refineSecNodes = 3
        self.jacobian = self.dijkstra.jacobian
        self.setThreadCount = self.dijkstra.setThreadCount
        #self.createJacobian = self.dijkstra.createJacobian
        self.setJacobian(self.dijkstra.jacobian())

    def regionManagerRef(self):
        # necessary because core dijkstra use its own RM
        return self.dijkstra.regionManagerRef()

    def createRefinedFwdMesh(self, mesh):
        """Refine the current mesh for higher accuracy.

        This is called automatic when accesing self.mesh() so it ensures any
        effect of changing region properties (background, single).
        """
        pg.info("Creating refined mesh (secnodes: {0}) to "
                "solve forward task.".format(self._refineSecNodes))
        m = mesh.createMeshWithSecondaryNodes(self._refineSecNodes)
        pg.verbose(m)
        return m

    def setMeshPost(self, mesh):
        """
        """
        # pg._r(mesh)
        self.dijkstra.setMesh(mesh)
        #self.dijkstra.setMesh(pg.Mesh(mesh))

    def setDataPost(self, data):
        """
        """
        # pg._r()
        self.dijkstra.setData(data)

    def createStartModel(self, dataVals):
        """
        """
        dists = shotReceiverDistances(self.data, full=True)
        aSlow = 1. / (dists / dataVals)

        # pg._r(self.regionManager().parameterCount())
        sm = pg.Vector(self.regionManager().parameterCount(),
                       pg.math.median(aSlow))
        return sm

    def createJacobian(self, par):
        if not self.mesh():
            pg.critical("no mesh")
        return self.dijkstra.createJacobian(par)

    def response(self, par):
        if not self.mesh():
            pg.critical("no mesh")
        return self.dijkstra.response(par)

    def drawModel(self, ax, model, **kwargs):
        kwargs['label'] = pg.unit('vel')
        ax, cBar = super(TravelTimeDijkstraModelling, self).drawModel(ax=ax,
                                                           model=model,
                                                           **kwargs)
        return ax, cBar

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

        self.inv.dataTrans = pg.trans.TransLog()

    def createForwardOperator(self, **kwargs):
        """Create default forward operator for Traveltime modelling.

        Your want your Manager use a special forward operator you can add them
        here on default Dijkstra is used.
        """
        fop = TravelTimeDijkstraModelling(**kwargs)
        return fop

    def dataCheck(self, data):
        """Return data from container"""
        if isinstance(data, pg.DataContainer):
            if not data.haveData('t'):
                pg.critical('Datacontainer have no "t" values.')
            return data['t']

        return data

    def errorCheck(self, err, dataVals):
        """Return relative error"""
        if isinstance(err, pg.DataContainer):
            if not err.haveData('err'):
                pg.error('Datacontainer have no "err" values. Fallback set to 0.01')
            return err['err'] / dataVals

        return err

    def setMesh(self, mesh, secNodes=0, ignoreRegionManager=False):
        """ """
        self.fop._refineSecNodes = secNodes
        if secNodes > 0:
            if ignoreRegionManager:
                mesh = self.fop.createRefinedFwdMesh(mesh)

        self.fop.setMesh(mesh, ignoreRegionManager=ignoreRegionManager)

    def simulate(self, slowness, scheme, mesh=None, secNodes=2,
                 noiseLevel=0.0, noiseAbs=0.0, **kwargs):
        """Simulate Traveltime measurements.

        Perform the forward task for a given mesh, a slowness distribution (per
        cell) and return data (traveltime) for a measurement scheme.

        Parameters
        ----------
        mesh : :gimliapi:`GIMLI::Mesh`
            Mesh to calculate for or use the last known mesh.
        slowness : array(mesh.cellCount()) | array(N, mesh.cellCount())
            Slowness distribution for the given mesh cells can be:

            * a single array of len mesh.cellCount()
            * a matrix of N slowness distributions of len mesh.cellCount()
            * a res map as [[marker0, res0], [marker1, res1], ...]
        scheme: :gimliapi:`GIMLI::DataContainer`
            Data measurement scheme needs 's' for shot and 'g' for geophone
            data token.
        secNodes: int [2]
            Number of refinement nodes to increase accuracy of the forward
            calculation.
        noiseLevel: float [0.0]
            Add relative noise to the simulated data. noiseLevel*100 in %
        noiseAbs: float [0.0]
            Add absolute noise to the simulated data in ms.

        Other Parameters
        ----------------
        returnArray: [False]
            Return only the calculated times.
        verbose: [self.verbose]
            Overwrite verbose level.
        **kwargs
            Additional kwargs ...

        Returns
        -------
        t : array(N, data.size()) | DataContainer
            The resulting simulated travel time values.
            Either one column array or matrix in case of slowness matrix.
        """
        verbose = kwargs.pop('verbose', self.verbose)

        fop = self.fop
        fop.data = scheme
        fop.verbose = verbose

        if mesh is not None:
            self.setMesh(mesh, secNodes=secNodes, ignoreRegionManager=True)

        if len(slowness) == self.fop.mesh().cellCount():
            if max(slowness) > 1.:
                pg.warn('slowness values larger than 1 ({0}), assuming velocity values .. building reciprocity.'.format(max(slowness)))
                t = fop.response(1./slowness)
            else:
                t = fop.response(slowness)
        else:
            print(self.fop.mesh())
            print("slowness: ", slowness)
            pg.critical("Simulate called with wrong slowness array.")

        ret = pg.DataContainer(scheme)
        ret.set('t', t)

        if noiseLevel > 0 or noiseAbs > 0:
            if not ret.allNonZero('err'):
                ret.set('t', t)
                err = noiseAbs + t * noiseLevel
                ret.set('err', err)

            pg.verbose("Absolute data error estimates (min:max) {0}:{1}".format(
                        min(ret('err')), max(ret('err'))))

            t += pg.math.randn(ret.size()) * ret('err')
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
        mesh = kwargs.pop('mesh', None)
        if mesh is not None:
            self.setMesh(mesh, secNodes=kwargs.pop('secNodes', 2))

        if 'limits' in kwargs:
            if kwargs['limits'][0] > 1:
                tmp = kwargs['limits'][0]
                kwargs['limits'][0] = 1.0 / kwargs['limits'][1]
                kwargs['limits'][1] = 1.0 / tmp
                pg.verbose('Switching velocity limits to slowness limits.',
                            kwargs['limits'])

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
        ...                     shotDistance=y+1)
        >>> data.set("t", pg.Vector(data.size(), 1.0))
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

    def rayCoverage(self):
        """ray coverage
        TODO little more
        """
        return self.fop.jacobian().transMult(np.ones(self.fop.data.size()))

    def standardizedCoverage(self):
        """standardized coverage vector (0|1) using neighbor info
        TODO little more
        """
        coverage = self.rayCoverage()
        C = self.fop.constraintsRef()
        return np.sign(np.absolute(C.transMult(C * coverage)))

    def showCoverage(self, ax=None, name='coverage', **kwargs):
        """shows the ray coverage in logscale"""
        if ax is None:
            fig, ax = plt.subplots()
            self.figs[name] = fig

        cov = self.rayCoverage()
        return pg.show(self.fop.paraDomain,
                       pg.log10(cov+min(cov[cov > 0])*.5), ax=ax,
                       coverage=self.standardizedCoverage(), **kwargs)
