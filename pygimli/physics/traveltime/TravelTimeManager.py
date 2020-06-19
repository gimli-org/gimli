#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Class for managing first arrival travel time inversions"""
import os
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.collections import LineCollection

import pygimli as pg
from pygimli.frameworks import MeshModelling
from pygimli.frameworks import MeshMethodManager

from pygimli.utils import getSavePath

from . raplot import showVA
from . ratools import shotReceiverDistances, createGradientModel2D


class TravelTimeDijkstraModelling(MeshModelling):
    """Forward modelling class for traveltime using Dijsktras method."""

    def __init__(self, **kwargs):
        self._core = pg.core.TravelTimeDijkstraModelling()

        super(TravelTimeDijkstraModelling, self).__init__(**kwargs)
        self._useGradient = None  # assumed to be [vTop, vBot] if set
        self._refineSecNodes = 3
        self.jacobian = self._core.jacobian
        self.setThreadCount = self._core.setThreadCount
        # self.createJacobian = self.dijkstra.createJacobian
        self.setJacobian(self._core.jacobian())

    @property
    def dijkstra(self):
        """Return current Dijkstra graph associated to mesh and model."""
        return self._core.dijkstra()

    def regionManagerRef(self):
        # necessary because core dijkstra use its own RM
        return self._core.regionManagerRef()

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
        self._core.setMesh(mesh)

    def setDataPost(self, data):
        """
        """
        self._core.setData(data)

    def createStartModel(self, dataVals):
        """
        """
        sm = None

        if self._useGradient is not None:
            [vTop, vBot] = self._useGradient  # something strange here!!!
            pg.info('Create gradient starting model. {0}: {1}'.format(vTop,
                                                                      vBot))
            sm = createGradientModel2D(self.data,
                                       self.paraDomain,
                                       vTop, vBot)
        else:
            dists = shotReceiverDistances(self.data, full=True)
            aSlow = 1. / (dists / dataVals)

            # pg._r(self.regionManager().parameterCount())
            sm = pg.Vector(self.regionManager().parameterCount(),
                           pg.math.median(aSlow))
            pg.info('Create constant starting model:', sm[0])

        return sm

    def createJacobian(self, par):
        """Create Jacobian (way matrix)."""
        if not self.mesh():
            pg.critical("no mesh")

        return self._core.createJacobian(par)

    def response(self, par):
        """Return forward response (simulated traveltimes)."""
        if not self.mesh():
            pg.critical("no mesh")
        return self._core.response(par)

    def way(self, s, g):
        """Return node indices for the way from the shot to the receiver.

        The index is based on the given data, mesh and last known model.
        """
        return self._core.way(s, g)

    def drawModel(self, ax, model, **kwargs):
        """Draw the model."""
        kwargs.setdefault('label', pg.unit('vel'))
        kwargs.setdefault('cMap', pg.utils.cMap('vel'))
        # kwargs['label'] = kwargs.pop('label', pg.unit('vel'))
        # kwargs['cMap'] = kwargs.pop('cMap', pg.utils.cMap('vel'))

        return super().drawModel(ax=ax, model=model,
                                 logScale=kwargs.pop('logScale', True),
                                 **kwargs)

    def drawData(self, ax, data=None, **kwargs):
        """Draw the data (as apparent velocity crossplot by default).

        Parameters
        ----------
        data: pg.DataContainer
        """
        kwargs.setdefault('label', pg.unit('va'))
        kwargs.setdefault('cMap', pg.utils.cMap('va'))

        if hasattr(data, '__iter__'):
            kwargs['vals'] = data
            data = self.data
        elif data is None:
            data = self.data

        if kwargs.pop('firstPicks', False):
            return pg.physics.traveltime.drawFirstPicks(ax, data)
        else:
            return showVA(data, usePos=False, ax=ax, **kwargs)


class TravelTimeManager(MeshMethodManager):
    """Manager for refraction seismics (traveltime tomography).

    TODO Document main members and use default MethodManager interface
    e.g., self.inv, self.fop, self.paraDomain, self.mesh, self.data
    """

    def __init__(self, data=None, **kwargs):
        """Create an instance of the Traveltime manager.

        Parameters
        ----------
        data: :gimliapi:`GIMLI::DataContainer` | str
            You can initialize the Manager with data or give them a dataset
            when calling the inversion.

        Other Parameters
        ----------------
        * useBert: bool [True]
            Use Bert forward operator instead of the reference implementation.
        * sr: bool [True]
            Calculate with singularity removal technique.
            Recommended but needs the primary potential.
            For flat earth cases the primary potential will be calculated
            analytically. For domains with topography the primary potential
            will be calculated numerical using a p2 refined mesh or
            you provide primary potentials with setPrimPot.
        """
        self._useFMM = False
        self.secNodes = 2  # default number of secondary nodes for inversion

        super(TravelTimeManager, self).__init__(data=data, **kwargs)

        self.inv.dataTrans = pg.trans.Trans()

    @property
    def velocity(self):
        # we can check here if there was an inversion run
        return self.fw.model

    def createForwardOperator(self, **kwargs):
        """Create default forward operator for Traveltime modelling.

        Your want your Manager use a special forward operator you can add them
        here on default Dijkstra is used.
        """
        fop = TravelTimeDijkstraModelling(**kwargs)
        return fop

    def load(self, fileName):
        self.data = pg.physics.traveltime.load(fileName)
        return self.data

    def createMesh(self, data=None, **kwargs):
        """Create default inversion mesh.

        Inversion mesh for traveltime inversion does not need boundary region.
        """
        d = data or self.data

        if d is None:
            pg.critical('Please provide a data file for mesh generation')

        return pg.meshtools.createParaMesh(data.sensors(),
                                           boundary=0, **kwargs)

    def checkData(self, data):
        """Return data from container."""
        if isinstance(data, pg.DataContainer):
            if not data.haveData('t'):
                pg.critical('DataContainer has no "t" values.')
            return data['t']

        return data

    def checkError(self, err, dataVals):
        """Return relative error."""
        if isinstance(err, pg.DataContainer):
            if not err.haveData('err'):
                pg.error('DataContainer has no "err" values. Fallback to 3%')
                return np.ones(err.size()) * 0.03
            return err['err'] / dataVals

        return err

    def applyMesh(self, mesh, secNodes=None, ignoreRegionManager=False):
        """Apply mesh, i.e. set mesh in the forward operator class."""
        if secNodes is None:
            secNodes = self.secNodes

        self.fop._refineSecNodes = secNodes
        if secNodes > 0:
            if ignoreRegionManager:
                mesh = self.fop.createRefinedFwdMesh(mesh)

        self.fop.setMesh(mesh, ignoreRegionManager=ignoreRegionManager)

    def simulate(self, mesh, scheme, slowness=None, vel=None, seed=None,
                 secNodes=2, noiseLevel=0.0, noiseAbs=0.0, **kwargs):
        """Simulate traveltime measurements.

        Perform the forward task for a given mesh, a slowness distribution (per
        cell) and return data (traveltime) for a measurement scheme.

        Parameters
        ----------
        mesh : :gimliapi:`GIMLI::Mesh`
            Mesh to calculate for or use the last known mesh.
        scheme: :gimliapi:`GIMLI::DataContainer`
            Data measurement scheme needs 's' for shot and 'g' for geophone
            data token.
        slowness : array(mesh.cellCount()) | array(N, mesh.cellCount())
            Slowness distribution for the given mesh cells can be:

            * a single array of len mesh.cellCount()
            * a matrix of N slowness distributions of len mesh.cellCount()
            * a res map as [[marker0, res0], [marker1, res1], ...]
        vel : array(mesh.cellCount()) | array(N, mesh.cellCount())
            Velocity distribution for the given mesh cells.
            Will overwrite given slowness.
        secNodes: int [2]
            Number of refinement nodes to increase accuracy of the forward
            calculation.
        noiseLevel: float [0.0]
            Add relative noise to the simulated data. noiseLevel*100 in %
        noiseAbs: float [0.0]
            Add absolute noise to the simulated data in ms.
        seed: int [None]
            Seed the random generator for the noise.

        Keyword Arguments
        -----------------
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
            self.applyMesh(mesh, secNodes=secNodes, ignoreRegionManager=True)

        if vel is not None:
            slowness = 1/vel

        if slowness is None:
            pg.critical("Need some slowness or velocity distribution for"
                        " simulation.")

        if len(slowness) == self.fop.mesh().cellCount():
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

            pg.verbose("Absolute error estimates (min:max) {0}:{1}".format(
                min(ret('err')), max(ret('err'))))

            t += pg.randn(ret.size(), seed=seed) * ret('err')
            ret.set('t', t)

        if kwargs.pop('returnArray', False) is True:
            return t

        return ret

    def invert(self, data=None, useGradient=True, vTop=500, vBottom=5000,
               secNodes=2, **kwargs):
        """Invert data.

        Parameters
        ----------
        data : pg.DataContainer()
            Data container with at least SensorIndieces 's g' and
            data values 't' (traveltime in ms) and 'err' (absolute error in ms)
        useGradient: bool [True]
            Use a gradient like starting model suited for standard flat
            earth cases. [Default]
            For cross tomography geometry you should set this to False for a
            non gradient startind model.
        vTop: float
            Top velocity for gradient stating model.
        vBottom: float
            Bottom velocity for gradient stating model.
        secNodes: int [2]
            Amount of secondary nodes used for ensure accuracy of the forward
            operator.

        Keyword Arguments
        -----------------
        ** kwargs:
            Inversion related arguments:
            See :py:mod:`pygimli.frameworks.MeshMethodManager.invert`
        """
        mesh = kwargs.pop('mesh', None)

        self.secNodes = secNodes

        if 'limits' in kwargs:
            if kwargs['limits'][0] > 1:
                tmp = kwargs['limits'][0]
                kwargs['limits'][0] = 1.0 / kwargs['limits'][1]
                kwargs['limits'][1] = 1.0 / tmp
                pg.verbose('Switching velocity limits to slowness limits.',
                           kwargs['limits'])

        if useGradient:
            self.fop._useGradient = [vTop, vBottom]
        else:
            self.fop._useGradient = None

        slowness = super().invert(data, mesh, **kwargs)
        velocity = 1.0 / slowness
        self.fw.model = velocity
        return velocity

    def drawRayPaths(self, ax, model=None, **kwargs):
        """Draw the the ray paths for model or last model.

        If model is not specifies, the last calculated Jacobian is used.

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
            if self.fop.jacobian().size() == 0 or model != self.model:
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
            wi = self.fop.way(s, g)
            points = self.fop._core.mesh().positions(withSecNodes=True)[wi]
            segs.append(np.column_stack((pg.x(points), pg.y(points))))

        lc = LineCollection(segs, **kwargs)
        ax.add_collection(lc)

        return lc

    def showRayPaths(self, model=None, ax=None, **kwargs):
        """Show the model with ray paths for given model.

        If not model specified, the last calculated Jacobian is taken.

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
        >>> tt.applyMesh(mesh, secNodes=10)
        >>> ax, cb = tt.showRayPaths(showMesh=True, diam=0.1)
        """
        if model is None:
            if self.fop.jacobian().size() == 0:
                self.fop.mesh()  # initialize any meshs .. just to be sure is 1
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
        """Ray coverage, i.e. summed raypath lengths."""
        return self.fop.jacobian().transMult(np.ones(self.fop.data.size()))

    def standardizedCoverage(self):
        """Standardized coverage vector (0|1) using neighbor info."""
        coverage = self.rayCoverage()
        C = self.fop.constraintsRef()
        return np.sign(np.absolute(C.transMult(C * coverage)))

    def showCoverage(self, ax=None, name='coverage', **kwargs):
        """Show the ray coverage in log-scale."""
        if ax is None:
            fig, ax = plt.subplots()

        cov = self.rayCoverage()
        return pg.show(self.fop.paraDomain,
                       pg.log10(cov+min(cov[cov > 0])*.5), ax=ax,
                       coverage=self.standardizedCoverage(), **kwargs)

    def saveResult(self, folder=None, size=(16, 10), verbose=False, **kwargs):
        """Save the results in a specified (or date-time derived) folder.

        Saved items are:
            * Resulting inversion model
            * Velocity vector
            * Coverage vector
            * Standardized coverage vector
            * Mesh (bms and vtk with results)

        Args
        ----
        path: str[None]
            Path to save into. If not set the name is automatically created
        size: (float, float) (16,10)
            Figure size.

        Keyword Args
        ------------
        Will be forwarded to showResults

        Returns
        -------
        str:
            Name of the result path.
        """
        subfolder = self.__class__.__name__
        path = getSavePath(folder, subfolder)

        if verbose:
            pg.info('Saving refraction data to: {}'.format(path))

        np.savetxt(os.path.join(path, 'velocity.vector'),
                   self.velocity)
        np.savetxt(os.path.join(path, 'velocity-cov.vector'),
                   self.rayCoverage())
        np.savetxt(os.path.join(path, 'velocity-scov.vector'),
                   self.standardizedCoverage())

        m = pg.Mesh(self.paraDomain)

        m['Velocity'] = self.paraModel(self.velocity)
        m['Coverage'] = self.rayCoverage()
        m['S_Coverage'] = self.standardizedCoverage()
        m.exportVTK(os.path.join(path, 'velocity'))
        m.saveBinaryV2(os.path.join(path, 'velocity-pd'))
        self.fop.mesh().save(os.path.join(path, 'velocity-mesh'))

        np.savetxt(os.path.join(path, 'chi.txt'), self.inv.chi2History)

        fig, ax = plt.subplots()
        self.showResult(ax=ax, cov=self.standardizedCoverage(), **kwargs)
        fig.set_size_inches(size)
        fig.savefig(os.path.join(path, 'velocity.pdf'), bbox_inches='tight')
        pg.plt.close(fig)
        return path
