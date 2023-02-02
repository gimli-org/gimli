#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Class for managing seismic refraction data and doing inversions"""

from math import pi
import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.collections import LineCollection


import pygimli as pg
import pygimli.meshtools as mt
from pygimli.viewer.mpl import drawModel, drawMesh, CellBrowser, createColorBar
from pygimli.utils.base import interperc, getSavePath
from pygimli.viewer.mpl.dataview import showVecMatrix

from pygimli.frameworks import MethodManager#, MethodManager0

# the explicit import with full name allow for:
# python ~/src/gimli/gimli/pygimli/physics/traveltime/refraction.py
from .utils import createGradientModel2D
from .plotting import drawTravelTimeData, drawFirstPicks  # , plotLines
from .importData import importGTT
from .modelling import FatrayDijkstraModelling


# class Refraction(MethodManager0):
#     """Manager for refraction seismics (traveltime tomography)

#     TODO Document main members and use default MethodeManager interface
#     e.g., self.inv, self.fop, self.paraDomain, self.mesh, self.data
#     """

#     def __init__(self, data=None, verbose=True, debug=False, fatray=False,
#                  frequency=1000., **kwargs):
#         """Init function with optional data load"""
#         pg.deprecated('Use TravelTimeManager insteas')
#         super().__init__(verbose=verbose, debug=debug, **kwargs)
#         self.figs = {}
#         self.axs = {}

#         self.doSave = kwargs.pop('doSave', False)
#         self.errIsAbsolute = True
#         self.method = None

#         # should be forwarded so it can be accessed from outside
#         self.mesh = None
#         self.poly = None
#         self.error = None
#         self.velocity = None
#         self.response = None
#         self.__dict__.update(**kwargs)
#         # self.start = []
#         self.pd = None

#         # CR!, check if this should be better a static member TG: no idea
#         self.dataToken_ = 't'
#         if fatray:
#             self.useFatray(True, frequency)
#         if isinstance(data, str):
#             self.loadData(data)
#         elif isinstance(data, pg.DataContainer):
#             self.setDataContainer(pg.DataContainer(data))
#             self.basename = kwargs.pop('name', 'new')
# #        if self.dataContainer is not None:
# #            self.createMesh()
# #        self.fop = self.createFOP(verbose=self.verbose)
# #        self.inv = self.createInv(self.fop,
# #                                  verbose=self.verbose, doSave=self.doSave)

# #    def __repr__(self):  # no need to overwrite with identical content
# #        """string representation of the class"""
# #        return self.__repr__()
# #
#     def __repr__(self):  # to be moved to Mesh/Data Method manager
#         """String representation of the class for the print function"""
#         out = type(self).__name__ + " object"
#         if hasattr(self, 'dataContainer'):
#             out += "\n" + self.dataContainer.__str__()
#         if hasattr(self, 'mesh'):
#             out += "\n" + self.mesh.__str__()
#         return out

#     def paraDomain(self):
#         """Return parameter domain mesh."""
#         return self.fop.regionManager().paraDomain()

#     def getModel(self):  # model collided with base method Manager attribute
#         """Return velocity vector."""
#         # (self.paraDomain.cellMarkers())
#         return self.velocity

#     def useFMM(self, fmm=True):
#         """Define whether to use Fast Marching Method (FMM).

#         Note that this method is more accurate but currently a lot slower!
#         """
#         self.fop = Refraction.createFOP(usefmm=fmm)

#     def useFatray(self, fatray=True, frequency=300.):
#         """Define whether to use Fatray jacobian computation."""
#         self.fop = Refraction.createFOP(fatray=fatray)
#         if fatray:
#             self.fop.frequency = frequency

#     @staticmethod
#     def createFOP(verbose=False, usefmm=False, fatray=False):
#         """Create default forward operator for Traveltime modelling.

#         usefmm forces Fast Marching Method, otherwise Dijkstra is used.
#         """
#         if usefmm:
#             from .FMModelling import TravelTimeFMM
#             fop = TravelTimeFMM(verbose=verbose)
#         else:
#             if fatray:
#                 fop = FatrayDijkstraModelling(verbose=verbose)
#             else:
#                 fop = pg.core.TravelTimeDijkstraModelling(verbose=verbose)

#         return fop

#     def createInv(self, fop, verbose=True, doSave=False):
#         """Create default inversion instance for Traveltime inversion."""
#         self.tD = pg.trans.Trans()
#         self.tM = pg.trans.TransLogLU()

#         inv = pg.Inversion(verbose, doSave)
#         inv.setTransData(self.tD)
#         inv.setTransModel(self.tM)
#         inv.setForwardOperator(fop)

#         return inv

#     def createApparentData(self, data):
#         """Create apparent slowness for given data."""
#         # hackish .. dislike!
#         self.setData(data)
#         return 1./(self.getOffset(data=data, full=True) / data('t'))

#     def dataVals(self, data):
#         """Return pure data values from a given DataContainer."""
#         return data('t')

#     def relErrorVals(self, data):
#         """Return pure data values from a given DataContainer."""
#         return data('err') / data('t')

#     def setData(self, data):
#         """Set data container (holding s and g indices and t floats)."""
#         if issubclass(type(data), pg.DataContainer):
#             self.setDataContainer(data)
#         else:
#             raise BaseException("Implement set data from type:", type(data))

#     def setDataContainer(self, data):
#         """Set data container from outside."""
#         self.dataContainer = data
#         self.checkData()
#         self.fop.setData(self.dataContainer)
#         self.inv.setData(self.dataContainer('t'))
#         if self.dataContainer.allNonZero('err'):
#             self.error = self.dataContainer('err')
#         else:
#             self.error = Refraction.estimateError(data)

#     def loadData(self, filename):
#         """Load data from file."""
#         if filename.endswith('.gtt'):
#             data = importGTT(filename)
#         else:
#             data = pg.DataContainer(filename, sensorTokens='s g')

#         self.basename = filename[:filename.rfind('.')]
#         self.setDataContainer(data)

#     def checkData(self):
#         """Check data

#         w.r.t. shot/geophone identity and zero/negative
#         traveltimes, plus check y/z sensor positions
#         """
#         oldsize = self.dataContainer.size()
#         self.dataContainer.markInvalid(pg.abs(self.dataContainer('s') -
#                                               self.dataContainer('g')) < 1)
#         self.dataContainer.markInvalid(self.dataContainer('t') <= 0.)
#         self.dataContainer.removeInvalid()
#         newsize = self.dataContainer.size()

#         if newsize < oldsize:
#             if self.verbose:
#                 print('Removed ' + str(oldsize - newsize) + ' values.')

#         maxyabs = max(pg.abs(pg.y(self.dataContainer.sensorPositions())))
#         maxzabs = max(pg.abs(pg.z(self.dataContainer.sensorPositions())))

#         if maxzabs > 0 and maxyabs == 0:
#             for i in range(self.dataContainer.sensorCount()):
#                 pos = self.dataContainer.sensorPosition(i).rotateX(-pi / 2)
#                 self.dataContainer.setSensorPosition(i, pos)

#         if self.verbose:
#             print(self.dataContainer)

#     def showData(self, data=None, response=None, ax=None, name='data'):
#         """Show data as travel time curves (optionally with response)

#         Parameters
#         ----------
#         data : pyGIMLi data Container [self.dataContainer]
#             data to show with points
#         response : array
#             response vector to draw with lines
#         ax : maxplotlib axes
#             axis to plot into, if not given, a new figure is created
#         """
#         if data is None:
#             data = self.dataContainer

#         if response is not None:
#             name = 'datafit'

#         if ax is None:
#             fig, ax = plt.subplots()
#             self.figs[name] = fig

#         self.axs[name] = ax
#         if response is None:
#             drawFirstPicks(ax, data)
#         else:
#             drawFirstPicks(ax, data, marker='+')
#             if response is True:
#                 response = self.response
#             drawFirstPicks(ax, data, np.asarray(response), marker='-')

#         return ax

#     def createMesh(self, depth=None, quality=34.3, paraDX=1, boundary=0,
#                    paraBoundary=0, secNodes=3, apply=True, **kwargs):
#         """Create (inversion) mesh using createParaDomain2D

#         Parameters
#         ----------
#         depth : float, optional
#             maximum depth, 0 (default) means maximum offset / 3.
#         paraDX : float
#             relative distance for refinement nodes between two sensors
#             e.g., 0 or 1 means no refinement
#             e.g., 0.5 means 1 additional node between two neighboring sensors
#             e.g., 0.33 means 2 additional equidistant nodes between two sensors
#         boundary : float, optional
#             boundary width to be appended for domain prolongation in absolute
#             para domain width.
#             values < 0 force the boundary to be 4 times para domain width.
#         paraBoundary : float, optional
#             margin for parameter domain in sensor distances (default 2)
#         quality : float, optional
#             mesh quality (smallest angle allowed)
#         apply : bool, optional
#             set mesh property of the underlying forward operator
#         secNodes : int (1)
#             Amount of secondary nodes to improve accuracy of the forward
#             solution.
#         **kwargs: Additional keyword arguments passed to
#             pygimli.meshtools.createParaMeshPLC

#         See also
#         --------
#         pygimli.meshtools.createParaMeshPLC
#         """
#         if self.dataContainer is None:
#             raise BaseException('Cannot create mesh without dataContainer.')

#         if depth is None:
#             depth = self.getDepth()

#         self.poly = mt.createParaMeshPLC(self.dataContainer.sensorPositions(),
#                                          paraDepth=depth, paraDX=paraDX,
#                                          paraBoundary=paraBoundary,
#                                          boundary=boundary, **kwargs)
#         mesh = mt.createMesh(self.poly, quality=quality, smooth=(1, 10))

#         if apply:
#             self.setMesh(mesh, secNodes=secNodes)

#         return mesh

#     def setMesh(self, mesh, refine=False, secNodes=1):
#         """Set mesh. To be removed from class once derived from MeshManager.

#         Parameters
#         ----------
#         secNodes : int (1)
#             Number of secondary nodes to improve accuracy of the forward
#             solution.
#         """
#         self.mesh = mesh
#         self.mesh.createNeighborInfos()
#         self.fop.setMesh(self.mesh)
#         self.fop.regionManager().setConstraintType(1)

#         if refine:
#             pg.warn("argument refine is deprecated .. use secnodes instead")
#             secNodes = 1

#         mesh = self.fop.regionManager().mesh().createMeshWithSecondaryNodes(
#                 secNodes)
#         self.fop.setMesh(mesh, ignoreRegionManager=True)

#         self.inv.setForwardOperator(self.fop)

#     def showMesh(self, ax=None, name='mesh'):
#         """show mesh in given ax or in a new figure"""
#         if ax is None:
#             fig, ax = plt.subplots()
#             self.figs[name] = fig

#         self.axs[name] = ax
#         drawMesh(ax, self.mesh)
# #        plt.show(block=False)
#         ax.set_aspect(1)

#         return ax

#     @staticmethod
#     def estimateError(data=None, absoluteError=0.001, relativeError=0.001):
#         """Estimate error composed of an absolute and a relative part

#         Parameters
#         ----------
#         absoluteError : float
#             absolute error of traveltimes (usually in s)
#         relativeError : float
#             relative error of traveltimes in 1 (e.g. 0.01 is 1%)

#         Returns
#         -------
#         err : array
#         """
#         # print(data)
#         # if not data.allNonZero('t'):
#         #     raise BaseException("We need travel time values (t) " +
#         #                         "in the data to estimate a data error.")

#         if relativeError >= 0.5:  # obviously in %
#             print("relativeError set to a value > 0.5 .. assuming this "
#                   "is a percentile error level dividing them by 100")
#             relativeError /= 100.

#         error = absoluteError + data('t') * relativeError
#         return error

#     def invert(self, data=None, t=None, err=None, mesh=None, **kwargs):
#         """Run actual inversion.

#         Values for result/response are stored in the class members
#         velocity/response

#         Parameters
#         ----------
#         useGradient : bool
#             Create gradient for starting model from vtop to vbottom.
#         vtop, vbottom : float
#             starting (gradient) model velocities on top/at bottom of the mesh
#         lam : float
#             regularization parameter describing the strength of smoothness
#         zWeight : float
#             relative weight for purely vertical boundaries
#         maxIter : int
#             Maximum number of iterations
#         startModel : array
#             Slowness starting model for the inversion
#         """
#         if 'verbose' in kwargs:
#             self.setVerbose(kwargs.pop('verbose'))

#         if data is not None:
#             # setDataContainer would be better
#             if t is not None:
#                 data.set('t', t)
#             self.setDataContainer(data)

#         if t is not None:
#             self.dataContainer.set('t', t)

#         if err is not None:
#             self.error = err

#         if mesh is not None:
#             self.setMesh(mesh)

#         if self.mesh is None:
#             self.createMesh(**kwargs)

#         startModel = kwargs.pop('startModel', None)
#         self.pd = self.fop.regionManager().paraDomain()

#         if startModel is None:
#             useGradient = kwargs.pop('useGradient', True)
#             if useGradient:
#                 startModel = createGradientModel2D(
#                     self.dataContainer, self.pd,
#                     kwargs.pop('vtop', 500.), kwargs.pop('vbottom', 5000.))
#             else:
#                 startModel = self.fop.createDefaultStartModel()
#         if isinstance(startModel, (float, int)):
#             startModel = pg.Vector(self.pd.cellCount(), startModel)

#         self.fop.setStartModel(startModel)

#         zWeight = kwargs.pop('zWeight', 0.2)
#         if 'zweight' in kwargs:
#             zWeight = kwargs.pop('zweight', 0.2)
#             print("zweight option will be removed soon. "
#                   "Please use zWeight instead.")
#         self.fop.regionManager().setZWeight(zWeight)

#         self.inv.setData(self.dataContainer('t'))
#         self.inv.setLambda(kwargs.pop('lam', 30.))

#         if 'threadCount' in kwargs:  # just for backward compatibility
#             self.fop.setThreadCount(kwargs.pop('threadCount'))
#         if 'max_iter' in kwargs:  # just for backward compatibility
#             self.inv.setMaxIter(kwargs.pop('max_iter'))
#         if 'maxIter' in kwargs:  # the better way
#             self.inv.setMaxIter(kwargs.pop('maxIter'))
#         if 'robustData' in kwargs:
#             self.inv.setRobustData(kwargs.pop('robustData'))
#         if 'blockyModel' in kwargs:
#             self.inv.setBlockyModel(kwargs.pop('blockyModel'))
#         if kwargs.pop('referenceModel', False):
#             self.inv.setReferenceModel(startModel)

#         if not hasattr(self.error, '__iter__'):
#             self.error = Refraction.estimateError(
#                 self.dataContainer, kwargs.pop('error', 0.003))  # abs err in s

#         self.inv.setAbsoluteError(self.error)
#         self.fop.jacobian().clear()

#         slowness = self.inv.run()
#         self.velocity = 1. / slowness
#         self.response = self.inv.response()

#         # use self.model() to access to this self.model = self.velocity

#         return self.velocity

#     @staticmethod
#     def simulate(mesh, slowness, scheme, verbose=False, **kwargs):
#         """Simulate a traveltime measurement.

#         Perform the forward task for a given mesh, a slowness distribution (per
#         cell) and return data (traveltime) for a measurement scheme. This is a
#         static method since it does not interfere with the managers inversion
#         approaches.

#         Parameters
#         ----------
#         mesh : :gimliapi:`GIMLI::Mesh`
#             Mesh to calculate for.

#         slowness : array(mesh.cellCount()) | array(N, mesh.cellCount())
#             slowness distribution for the given mesh cells can be:
#             . a single array of len mesh.cellCount()
#             . a matrix of N slowness distributions of len mesh.cellCount()
#             . a res map as [[marker0, res0], [marker1, res1], ...]

#         scheme : :gimliapi:`GIMLI::DataContainer`
#             data measurement scheme

#         verbose : boolean
#             Be verbose.

#         Other parameters
#         ----------------
#         noisify : boolean
#             add normally distributed noise based on scheme('err')

#         Returns
#         -------
#         t : array(N, data.size()) | DataContainer
#             The resulting simulated travel time values.
#             Either one column array or matrix in case of slowness matrix.
#             A DataContainer is return if noisify set to True.

#         """
#         fop = Refraction.createFOP(verbose=verbose)

#         fop.setData(scheme)
#         fop.setMesh(mesh, ignoreRegionManager=True)

#         if len(slowness) == mesh.cellCount():
#             if max(slowness) > 1.:
#                 print('Warning: slowness values larger than 1 (' +
#                       str(max(slowness)) + ').. assuming that are velocity '
#                       'values .. building reciprocity')
#                 t = fop.response(1./slowness)
#             else:
#                 t = fop.response(slowness)
#         else:
#             print(mesh)
#             print("slowness: ", slowness)
#             raise BaseException("Simulate called with wrong slowness array.")

#         ret = pg.DataContainer(scheme)
#         ret.set('t', t)

#         noiseLevel = kwargs.pop('noiseLevel', 0)
#         noiseAbs = kwargs.pop('noiseAbs', 0)

#         if noiseLevel > 0 or noiseAbs > 0:
#             if not ret.allNonZero('err'):
#                 ret.set('t', t)
#                 ret.set('err', pg.physics.Refraction.estimateError(
#                     ret, absoluteError=noiseAbs))

#             if verbose:
#                 print("Data error estimates (min:max) ",
#                       min(ret('err')), ":", max(ret('err')))

#             t += pg.math.randn(ret.size()) * ret('err')
#             ret.set('t', t)

#         if kwargs.pop('returnArray', False):
#             return t

#         return ret

#     @staticmethod
#     def drawApparentVelocities(ax, data, t=None, **kwargs):
#         """Plot data in for of apparent velocity image."""
#         tt = Refraction()
#         tt.setDataContainer(data)
#         tt.showVA(ax=ax, t=t, **kwargs)

#     @staticmethod
#     def drawTravelTimeData(ax, data, t=None):
#         """Plot travel time data as lines and points."""
#         drawTravelTimeData(ax, data, t)

#     def getOffset(self, data=None, full=False):
#         """Return vector of offsets (in m) between shot and receiver."""
#         if data is None:
#             data = self.dataContainer

#         if full:
#             pos = data.sensorPositions()
#             s, g = data('s'), data('g')
#             nd = data.size()
#             off = [pos[int(s[i])].distance(pos[int(g[i])]) for i in range(nd)]
#             return np.absolute(off)
#         else:
#             px = pg.x(data.sensorPositions())
#             gx = np.array([px[int(g)] for g in data("g")])
#             sx = np.array([px[int(s)] for s in data("s")])
#             return np.absolute(gx - sx)

#     def getMidpoint(self, data=None):
#         """Return vector of offsets (in m) between shot and receiver."""
#         if data is None:
#             data = self.dataContainer

#         px = pg.x(data.sensorPositions())
#         gx = np.array([px[int(g)] for g in data("g")])
#         sx = np.array([px[int(s)] for s in data("s")])
#         return (gx + sx) / 2

#     def showVA(self, data=None, t=None, name='va', pseudosection=False,
#                squeeze=True, full=True, ax=None, cmap=None, **kwargs):
#         """Show apparent velocity as image plot.

#         TODO showXXX commands need to return ax and cbar .. if there is one
#         """
#         if data is None:
#             data = self.dataContainer

#         if ax is None:
#             fig, ax = plt.subplots()
#             self.figs[name] = fig

#         self.axs[name] = ax
#         if t is None:
#             t = data('t')

#         px = pg.x(data.sensorPositions())
#         py = pg.y(data.sensorPositions())
#         if len(np.unique(py)) > len(np.unique(px)):  # probably crosshole
#             px = py

#         gx = np.array([px[int(g)] for g in data("g")])
#         sx = np.array([px[int(s)] for s in data("s")])
#         offset = self.getOffset(data=data, full=full)
#         kwargs.setdefault('vals', offset / t)

#         if pseudosection:
#             midpoint = (gx + sx) / 2
#             _, cb = showVecMatrix(midpoint, offset, squeeze=True, ax=ax,
#                                   label='Apparent slowness [s/m]', cmap=cmap,
#                                   **kwargs)
#         else:
#             _, cb = showVecMatrix(gx, sx, squeeze=squeeze, ax=ax,
#                                   label='Apparent velocity [m/s]', cmap=cmap,
#                                   **kwargs)
#         ax.figure.show()
#         return ax, cb

#     def getDepth(self):
#         """return a (a-priori guessed) depth of investigation"""
#         return max(self.getOffset()) / 3.0  # rule of thumb

#     def rayCoverage(self):
#         """return ray coverage"""
#         one = pg.Vector(self.dataContainer.size(), 1.)
#         return self.fop.jacobian().transMult(one)

#     def standardizedCoverage(self):
#         """return standardized coverage vector (0|1) using neighbor info"""
#         coverage = self.rayCoverage()
#         C = self.fop.constraintsRef()
#         return np.sign(np.absolute(C.transMult(C * coverage)))

#     def showRayPaths(self, model=None, ax=None, **kwargs):
#         """Show ray paths for `model` or last model for which the Jacobian was
#         calculated.

#         Parameters
#         ----------
#         model : array
#             Velocity model for which to calculate and visualize ray paths (the
#             default is model for last Jacobian calculation in self.velocity).
#         ax : matplotlib.axes
#             Axes for the plot (the default is None).
#         **kwargs : type
#             Additional arguments passed to LineCollection (alpha, linewidths,
#             color, linestyles).

#         Returns
#         -------
#         ax : matplotlib.axes object
#         cb : matplotlib.colorbar object (only if model is provided)
#         """
#         cbar = None
#         if model is None and self.velocity is None:
#             pg.info("No previous inversion result found and no model given.",
#                     "Using homogeneous slowness model.")
#             self.velocity = pg.Vector(self.mesh.cellCount(), 1.0)
#             self.fop.createJacobian(1./self.velocity)

#         if model is not None:
#             if self.velocity is not None:
#                 if not np.allclose(model, self.velocity):
#                     self.fop.createJacobian(1/model)

#             ax, cbar = self.showResult(ax=ax, val=model)
#             _ = kwargs.setdefault("color", "w")
#             _ = kwargs.setdefault("alpha", 0.5)
#             _ = kwargs.setdefault("linewidths", 0.8)
#         else:
#             ax = self.showMesh(ax=ax)

#         # Due to different numbering scheme of way matrix
#         _, shots = np.unique(self.dataContainer("s"), return_inverse=True)
#         _, receivers = np.unique(self.dataContainer("g"), return_inverse=True)

#         # Collecting way segments for all shot/receiver combinations
#         segs = []
#         for s, g in zip(shots, receivers):
#             wi = self.fop.way(s, g)
#             points = self.fop.mesh().positions(withSecNodes=True)[wi]
#             segs.append(np.column_stack((pg.x(points), pg.y(points))))

#         line_segments = LineCollection(segs, **kwargs)
#         ax.add_collection(line_segments)
#         return ax, cbar

#     def showCoverage(self, ax=None, name='coverage', **kwargs):
#         """shows the ray coverage in logscale"""
#         if ax is None:
#             fig, ax = plt.subplots()
#             self.figs[name] = fig

#         self.axs[name] = ax
#         cov = self.rayCoverage()
#         return pg.show(self.mesh, pg.log10(cov+min(cov[cov > 0])*.5), ax=ax,
#                        coverage=self.standardizedCoverage(), **kwargs)

#     def showModel(self, ax=None, vals=None, **kwargs):
#         """WRITEME"""
#         return self.showResult(ax=ax, val=vals, **kwargs)

#     def showResult(self, val=None, ax=None, cMin=None, cMax=None,
#                    logScale=False, rays=False, name='result', **kwargs):
#         """Show resulting velocity vector.

#         Parameters
#         ----------
#         val : result array [self.velocity]
#             field to show, usually the velocity vector
#         ax : matplotlib.axes
#             axes to plot into, if not give a new one-axis figure is created
#         cMin/cMax : float
#             minimum and maximum values for ranging colorbar
#         logScale : bool [False]
#             use logarithmic scale
#         rays : bool [False]
#             Show ray paths as well.

#         Other parameters
#         ----------------
#         useCoverage : bool [True]
#             use standardized (0 or 1) ray coverage as alpha-shading
#         label : str
#             label to write on colorbar
#         orientaton : str
#             color bar orientation
#         nLevs : int [7]
#             number of level values for color bar
#         **kwargs : keyword arguments passed to the show function

#         Returns
#         -------
#         ax : maxplotlib axes

#         cb : matplotlib color bar object
#         """
#         mesh = self.paraDomain()
#         if val is None:
#             val = self.velocity
#         if cMin is None or cMax is None:
#             cMin, cMax = interperc(val, 3)
#         coverage = 1
#         if kwargs.pop('useCoverage', True):
#             coverage = self.standardizedCoverage()
#         label = kwargs.pop("label", "Velocity (m/s)")
#         if ax is None:
#             fig, ax = plt.subplots()
#             self.figs[name] = fig
#             ax, cbar = pg.show(mesh, val, logScale=logScale, ax=ax,
#                                colorBar=True, cMin=cMin, cMax=cMax,
#                                coverage=coverage, label=label, hold=True,
#                                **kwargs)
#             self.figs[name] = plt.gcf()
#         else:
#             gci = drawModel(ax, mesh, val, logScale=logScale,
#                             colorBar=True, cMin=cMin, cMax=cMax,
#                             coverage=coverage, **kwargs)
#             labels = ['cMin', 'cMax', 'nLevs', 'orientation']
#             subkwargs = {key: kwargs[key] for key in labels if key in kwargs}
#             cbar = createColorBar(gci, label=label, **subkwargs)
#         if rays:
#             self.showRayPaths(ax=ax, alpha=0.5, color="w", linewidths=.8)

#         browser = CellBrowser(self.mesh, val, ax)
#         browser.connect()

#         self.axs[name] = ax
#         if 'lines' in kwargs:
#             plotLines(ax, kwargs['lines'])
#         return ax, cbar

#     def showResultAndFit(self, name='resultfit', **kwargs):
#         """show two vertical subplots with result and data (with response)"""
#         fig, axs = plt.subplots(nrows=2)
#         self.figs[name] = fig
#         self.showResult(ax=axs[0], **kwargs)
#         self.showData(ax=axs[1], response=self.response)

#     def saveFigures(self, name=None, ext='pdf'):
#         """save all existing figures to files"""
#         if name is None:
#             name = self.basename
#         if name is None or not any(name):
#             name = 'out'
#         for key in self.figs:
#             self.figs[key].savefig(name+'-'+key+'.'+ext, bbox_inches='tight')

#     def makeJacobianPDF(self, ind=None):
#         """Make multipage Jacobian PDF."""
#         from matplotlib.backends.backend_pdf import PdfPages
#         if ind is None:
#             ind = range(self.dataContainer.size())
#         with PdfPages(self.basename+'-jacobian.pdf') as pdf:
#             fig, ax = pg.plt.subplots()
#             for ii in ind:
#                 jj = self.fop.jacobian().row(ii)
#                 pg.show(self.mesh, jj, ax=ax, coverage=(jj > 0))
#                 fig.savefig(pdf, format='pdf')
#                 ax.cla()

#     def saveResult(self, folder=None, size=(16, 10), **kwargs):
#         """Save the results in the specified folder.

#         Saved items are:
#             Inverted profile
#             Velocity vector
#             Coverage vector
#             Standardized coverage vector
#             Mesh (bms and vtk with results)
#         """
#         # TODO: How to extract the chi2 etc. from each iteration???

#         subfolder = '/' + self.__class__.__name__
#         path = getSavePath(folder, subfolder)

#         if self.verbose:
#             print('Saving refraction data to: {}'.format(path))

#         np.savetxt(path + '/velocity.vector',
#                    self.velocity)
#         np.savetxt(path + '/velocity-cov.vector',
#                    self.rayCoverage())
#         np.savetxt(path + '/velocity-scov.vector',
#                    self.standardizedCoverage())

#         self.mesh.addExportData('Velocity', self.velocity)
#         self.mesh.addExportData('Coverage', self.rayCoverage())
#         self.mesh.addExportData('S_Coverage', self.standardizedCoverage())
#         self.mesh.exportVTK(path + 'velocity')
#         self.mesh.save(path + 'velocity-mesh')
#         self.pd.save(path + 'velocity-pd')

#         fig, ax = plt.subplots()
#         self.showResult(ax=ax, cov=self.standardizedCoverage(), **kwargs)
#         fig.set_size_inches(size)
#         fig.savefig(path + '/velocity.pdf', bbox_inches='tight')

#         return path, fig, ax


# def test_Refraction():
#     """Test Refraction manager stability some data/mesh set / data update"""
#     import os
#     datafile = os.path.dirname(__file__) + '/example_topo.sgt'

#     ra = Refraction(datafile, verbose=False, doSave=False)
#     ra.createMesh(depth=80)
#     ra.inv.setMaxIter(1)

#     ra.invert()
#     m1 = ra.model()
#     mesh = pg.Mesh(ra.mesh)

#     ra.setMesh(mesh)
#     ra.invert()
#     m2 = ra.model()

#     np.testing.assert_array_equal(m1, m2)

#     ra.setData(pg.DataContainer(datafile, 's g'))
#     m3 = ra.invert()

#     np.testing.assert_array_equal(m1, m3)


# class Tomography(Refraction):

#     """Traveltime tomography for tomographic (e.g. crosshole) measurements"""

#     def __init__(self, data=None, tcorr=0, name='new', **kwargs):
#         """Init function with optional data load

#         Parameters
#         ----------
#         data : pg.DataContainer or string

#         tcorr : float [0]
#             correct travel times by common shift

#         name : str [data if being string, otherwise 'new']
#             basename for saving Figures, results etc.

#         ndig : int [2]
#             number of digits to round positions (e.g. 2=cm), alternatively:

#         roundto : float [0]
#             unit spacing to round positions on
#         """
#         if isinstance(data, str):
#             name = data[:data.rfind('.')]
#             if data.lower()[-4:] == '.tom':
#                 data = readTOMfile(data, **kwargs)
#             else:
#                 data = pg.DataContainer(data, 's g')
#         if tcorr != 0:
#             data.set('t', data('t') + tcorr)
#         super(Tomography, self).__init__(data, name=name, **kwargs)

#     def createMesh(self, quality=34.6, maxarea=0.1, addpoints=None):
#         """Create (inversion) mesh by circumventing PLC"""
#         data = self.dataContainer
#         sx = list(pg.x(data.sensorPositions()))
#         sz = list(pg.y(data.sensorPositions()))

#         if addpoints is not None:
#             for po in addpoints:
#                 sx.append(po[0])
#                 sz.append(po[1])

#         iS = np.argsort(np.arctan2(sx-np.mean(sx), sz-np.mean(sz)))
#         plc = pg.Mesh(2)
#         nodes = [plc.createNode(sx[i], sz[i], 0) for i in iS]
#         for i in range(len(nodes)-1):
#             plc.createEdge(nodes[i], nodes[i+1])

#         plc.createEdge(nodes[-1], nodes[0])
#         tri = pg.core.TriangleWrapper(plc)
#         tri.setSwitches("-pzFq"+str(quality)+"a"+str(maxarea))
#         self.setMesh(tri.generate())

#     def offset(self):
#         """return shot-geophone distance"""
#         data = self.dataContainer
#         return np.array([data.sensorPosition(int(data('g')[i])).distance(
#           data.sensorPosition(int(data('s')[i]))) for i in range(data.size())])

#     def getVA(self, t=None):
#         """return apparent velocity"""
#         if t is None:
#             t = self.dataContainer('t')
#         return self.offset() / t

#     def createStartModel(self):
#         """create (gradient) starting model with vtop/vbottom bounds"""
#         va = self.getVA()
#         nModel = self.fop.regionManager().parameterCount()
#         return pg.Vector(nModel, 1./np.mean(va))

#     def showVA(self, t=None, ax=None, usepos=True, name='va', squeeze=True):
#         """show apparent velocity as image plot"""
#         # va = self.getVA(vals=vals)
#         xvec = self.dataContainer('g')
#         yvec = self.dataContainer('s')
#         if usepos:
#             pz = pg.y(self.dataContainer.sensorPositions())
#             if squeeze:
#                 xvec = pz[xvec]
#                 yvec = pz[yvec]
#             else:
#                 pz = pg.y(self.dataContainer.sensorPositions())
#                 raise Exception('Implement ME')
#                 # xvec = px[xvec]*1000 + pz[xvec]
#                 # xvec = px[yvec]*1000 + pz[yvec]

#         plotVecMatrix(xvec, yvec, vals=t, squeeze=squeeze, ax=ax, name=name)

#     def showVAold(self, vals=None, ax=None, usepos=True, name='va'):
#         """show apparent velocity as image plot (old style)"""
#         va = self.getVA(t=vals)
#         data = self.dataContainer
#         A = np.ones((data.sensorCount(), data.sensorCount())) * np.nan
#         for i in range(data.size()):
#             A[int(data('s')[i]), int(data('g')[i])] = va[i]

#         if ax is None:
#             fig, ax = plt.subplots()
#             self.figs[name] = fig

#         gci = ax.imshow(A, interpolation='nearest')
#         ax.grid(True)
#         if usepos:
#             xt = np.linspace(0, data.sensorCount()-1, 7)
#             xt.round()
#             px = pg.abs(pg.y(self.dataContainer.sensorPositions()))
#             ax.set_xticks(xt)
#             ax.set_xticklabels([str(int(px[int(xti)])) for xti in xt])
#             ax.set_yticks(xt)
#             ax.set_yticklabels([str(int(px[int(xti)])) for xti in xt])

#         plt.colorbar(gci, ax=ax)
#         return va


# def main():
#     """Main"""
#     parser = MethodManager.createArgParser(dataSuffix='sgt')
#     options = parser.parse_args()

#     ra = Refraction(verbose=not options.quiet, debug=pg.debug())

#     ra.loadData(options.dataFileName)

#     ra.showData()
#     ra.showVA()
#     ra.createMesh(depth=options.depth)
#     ra.showMesh()
#     ra.invert(lam=options.lam, max_iter=options.maxIter,
#               robustData=options.robustData, blockyModel=options.blockyModel)
#     ra.showResult()


# if __name__ == '__main__':
#     main()
#     pg.wait()
