#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Class for managing seismic refraction data and doing inversions"""

from math import pi
import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.meshtools import createParaMeshPLC, createMesh
from pygimli.mplviewer import drawModel, drawMesh, CellBrowser, createColorBar
from pygimli.utils.base import interperc, getSavePath
from pygimli.mplviewer.dataview import plotVecMatrix

from pygimli.manager import MethodManager

# the explicit import with full name allow for:
# python ~/src/gimli/gimli/python/pygimli/physics/traveltime/refraction.py
from pygimli.physics.traveltime.ratools import createGradientModel2D
from pygimli.physics.traveltime.raplot import plotFirstPicks, plotLines

from . raplot import drawTravelTimeData


class Refraction(MethodManager):
    """Class for managing refraction seismics data

        TODO Document main members and use default MethodeManager interface
        e.g., self.inv, self.fop, self.paraDomain, self.mesh, self.data
    """

    def __init__(self, data=None, verbose=False, debug=False, **kwargs):
        """Init function with optional data load"""
        MethodManager.__init__(self, verbose=verbose, debug=debug, **kwargs)
        self.figs = {}
        self.axs = {}

        self.doSave = kwargs.pop('doSave', False)
        self.errIsAbsolute = True

        # should be forwarded so it can be accessed from outside
        self.mesh = None
        self.poly = None
        self.error = None
        self.velocity = None
        self.response = None
        # self.start = []
        self.pd = None

        #CR!, check if this should be better a static member
        self.dataToken_ = 't'

        if isinstance(data, str):
            self.loadData(data)
        elif isinstance(data, pg.DataContainer):
            self.setDataContainer(data)
            self.basename = kwargs.pop('name', 'new')
#        if self.dataContainer is not None:
#            self.createMesh()
#        self.fop = self.createFOP(verbose=self.verbose)
#        self.inv = self.createInv(self.fop,
#                                  verbose=self.verbose, doSave=self.doSave)

#    def __str__(self):  # no need to overwrite with identical content
#        """string representation of the class"""
#        return self.__repr__()
#
    def __repr__(self):  # to be moved to Mesh/Data Method manager
        """string representation of the class for the print function"""
        out = type(self).__name__ + " object"
        if hasattr(self, 'dataContainer'):
            out += "\n" + self.dataContainer.__str__()
        if hasattr(self, 'mesh'):
            out += "\n" + self.mesh.__str__()
        return out

    def paraDomain(self):
        """base api """
        return self.fop.regionManager().paraDomain()

    def model(self):
        """base api """
        # (self.paraDomain.cellMarkers())
        return self.velocity

    @staticmethod
    def createFOP(verbose=False):
        """Create default forward operator for Traveltime modelling.
        base api

        Dijkstra, later FMM.
        """
#        if not hasattr(self, 'mesh'):  # self.mesh is None:
#        self.createMesh()
        fop = pg.TravelTimeDijkstraModelling(verbose=verbose)
        return fop

    def createInv(self, fop, verbose=True, doSave=False):
        """Create default inversion instance for Traveltime inversion.

        base api
        (typically done by run)
        """
        self.tD = pg.RTrans()
        self.tM = pg.RTransLogLU()

        inv = pg.RInversion(verbose, doSave)
        inv.setTransData(self.tD)
        inv.setTransModel(self.tM)
        inv.setForwardOperator(fop)

        return inv

    def createApparentData(self, data):
        """
        Create apparent slowness for given data.
        """
        # hackish .. dislike!
        self.setData(data)
        return 1./(self.getOffset(data=data, full=True) / data('t'))

    def dataVals(self, data):
        """Return pure data values from a given DataContainer. """
        return data('t')

    def relErrorVals(self, data):
        """Return pure data values from a given DataContainer. """
        return  data('err') / data('t')

    def setData(self, data):
        """Set data """
        if issubclass(type(data), pg.DataContainer):
            self.setDataContainer(data)
        else:
            raise BaseException("Implement set data from type:", type(data))

    def setDataContainer(self, data):
        """set data container from outside

        base api

        """
        self.dataContainer = data
        self.checkData()
        self.fop.setData(self.dataContainer)
        self.inv.setData(self.dataContainer('t'))
        if self.dataContainer.allNonZero('err'):
            self.error = self.dataContainer('err')
        else:
            self.error = Refraction.estimateError(data)

    def loadData(self, filename):
        """load data from file"""
        # check for file formats and import if necessary
        data = pg.DataContainer(filename, 's g')
        self.basename = filename[:filename.rfind('.')]
        self.setDataContainer(data)

    def checkData(self):
        """check data w.r.t. shot/geophone identity and zero/negative
        traveltimes, plus check y/z sensor positions """
        oldsize = self.dataContainer.size()
        self.dataContainer.markInvalid(pg.abs(self.dataContainer('s') -
                                              self.dataContainer('g')) < 1)
        self.dataContainer.markInvalid(self.dataContainer('t') <= 0.)
        self.dataContainer.removeInvalid()
        newsize = self.dataContainer.size()

        if newsize < oldsize:
            if self.verbose:
                print('Removed ' + str(oldsize - newsize) + ' values.')

        maxyabs = max(pg.abs(pg.y(self.dataContainer.sensorPositions())))
        maxzabs = max(pg.abs(pg.z(self.dataContainer.sensorPositions())))

        if maxzabs > 0 and maxyabs == 0:
            for i in range(self.dataContainer.sensorCount()):
                pos = self.dataContainer.sensorPosition(i).rotateX(-pi / 2)
                self.dataContainer.setSensorPosition(i, pos)

        if self.verbose:
            print(self.dataContainer)

    def showData(self, data=None, response=None, ax=None, name='data'):
        """Show data as travel time curves (optionally with response)

        Parameters
        ----------

        """
        if data is None:
            data = self.dataContainer

        if response is not None:
            name = 'datafit'

        if ax is None:
            fig, ax = plt.subplots()
            self.figs[name] = fig

        self.axs[name] = ax
        if response is None:
            drawFirstPicks(ax, data)
        else:
            drawFirstPicks(ax, data, marker='+')
            if response is True:
                response = self.response
            drawFirstPicks(ax, data, np.asarray(response), marker='-')

        # CR: don't use plt.show(..) at all TODO
        plt.show(block=False)
        return ax

    def createMesh(self, depth=None, quality=34.3, paraDX=0.5, boundary=0,
                   paraBoundary=5, apply=True, **kwargs):
        """Create (inversion) mesh using createParaDomain2D

        Parameters
        ----------

        apply : bool
            set Mesh property of the underlying forward operator

        """

        if self.dataContainer is None:
            raise BaseException('Cannot create mesh without dataContainer.')

        if depth is None:
            depth = self.getDepth()
        self.poly = createParaMeshPLC(self.dataContainer.sensorPositions(),
                                      paraDepth=depth, paraDX=paraDX,
                                      paraBoundary=paraBoundary,
                                      boundary=boundary)
        mesh = createMesh(self.poly, quality=quality, smooth=(1, 10))
        if apply:
            self.setMesh(mesh)
        return mesh

    def setMesh(self, mesh, refine=True):
        """
        base api
        """
        self.mesh = mesh
        self.mesh.createNeighbourInfos()
        self.fop.setMesh(self.mesh)
        self.fop.regionManager().setConstraintType(1)
        self.fop.createRefinedForwardMesh(refine)
        self.inv.setForwardOperator(self.fop)

    def showMesh(self, ax=None, name='mesh'):
        """show mesh in given ax or in a new figure"""
        if ax is None:
            fig, ax = plt.subplots()
            self.figs[name] = fig

        self.axs[name] = ax
        drawMesh(ax, self.mesh)
#        plt.show(block=False)
        ax.set_aspect(1)

        return ax

    @staticmethod
    def estimateError(data=None, absoluteError=0.001, relativeError=0.001):
        """Estimate error composed of an absolute and a relative part

        Parameters
        ----------
        absoluteError : float
            absolute error of traveltimes (usually in s)
        relativeError : float
            relative error of traveltimes in 1 (e.g. 0.01 is 1%)

        Returns
        -------
        err : array
        """
        # print(data)
        # if not data.allNonZero('t'):
        #     raise BaseException("We need travel time values (t) " +
        #                         "in the data to estimate a data error.")

        if relativeError >= 0.5:  # obviously in %
            print("relativeError set to a value > 0.5 .. assuming this "
                "is a percentile Error level dividing them by 100")
            relativeError /= 100.

        error = absoluteError + data('t') * relativeError
        return error

    def invert(self, data=None, t=None, err=None, mesh=None, **kwargs):
        """Run actual inversion (first creating inversion object if not there)

        result/response is stored in the class attribute velocity/response

        Parameters
        ----------
        vtop, vbottom : float
            starting (gradient) model velocities on top/at bottom of the mesh

        lam : float
            regularization parameter describing the strength of smoothness

        zWeight : float
            relative weight for purely vertical boundaries
        """
        if 'verbose' in kwargs:
            self.setVerbose(kwargs.pop('verbose'))

        if data is not None:
            # setDataContainer would be better
            if t is not None:
                data.set('t', t)
            self.setDataContainer(data)

        if t is not None:
            self.dataContainer.set('t', t)

        if err is not None:
            self.error = err

        if mesh is not None:
            self.setMesh(mesh)

        if self.mesh is None:
            self.createMesh(**kwargs)

        useGradient = kwargs.pop('useGradient', True)
        if useGradient:
            self.fop.setStartModel(createGradientModel2D(
                self.dataContainer, self.fop.regionManager().paraDomain(),
                kwargs.pop('vtop', 500.), kwargs.pop('vbottom', 5000.)))
        else:
            self.fop.setStartModel(self.fop.createDefaultStartModel())

        zWeight = kwargs.pop('zWeight', 0.2)
        if 'zweight' in kwargs:
            zWeight = kwargs.pop('zweight', 0.2)
            print("zweight option will be removed soon. "
                  "Please use zWeight instead.")
        self.fop.regionManager().setZWeight(zWeight)

        self.inv.setData(self.dataContainer('t'))
        self.inv.setLambda(kwargs.pop('lam', 30.))
        if 'max_iter' in kwargs:  # just for backward compatibility
            self.inv.setMaxIter(kwargs.pop('max_iter'))
        if 'maxIter' in kwargs:  # the better way
            self.inv.setMaxIter(kwargs.pop('maxIter'))
        if 'robustData' in kwargs:
            self.inv.setRobustData(kwargs.pop('robustData'))
        if 'blockyModel' in kwargs:
            self.inv.setBlockyModel(kwargs.pop('blockyModel'))

        if not hasattr(self.error, '__iter__'):
            self.error = Refraction.estimateError(
                self.dataContainer, kwargs.pop('error', 0.003))  # abs err in s

        self.inv.setAbsoluteError(self.error)
        self.fop.jacobian().clear()

        slowness = self.inv.run()
        self.velocity = 1. / slowness
        self.response = self.inv.response()

        # use self.model() to access to this self.model = self.velocity

        return self.velocity

    @staticmethod
    def simulate(mesh, slowness, scheme, verbose=False, **kwargs):
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

        **kwargs :
            * noisify : add normal distributed noise based on scheme('err')
                IMPLEMENTME

        Returns
        -------
        t : array(N, data.size()) | DataContainer
            The resulting simulated travel time values.
            Either one column array or matrix in case of slowness matrix.
            A DataContainer is return if noisify set to True.

        """

        fop = Refraction.createFOP(verbose=verbose)

        fop.setData(scheme)
        fop.setMesh(mesh, ignoreRegionManager=True)

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
                ret.set('err', pg.physics.Refraction.estimateError(ret,
                                  absoluteError=kwargs.pop('noiseAbs', 1e-4),
                                  relativeError=noiseLevel))

            if verbose:
                print("Data error estimates (min:max) ",
                      min(ret('err')), ":", max(ret('err')))

            t += pg.randn(ret.size()) * ret('err')
            ret.set('t', t)

        if kwargs.pop('returnArray', False):
            return t

        return ret

    @staticmethod
    def drawTravelTimeData(ax, data, t=None):
        """WRITEME"""
        drawTravelTimeData(ax, data, t)

    @staticmethod
    def drawApparentVelocities(ax, data, t=None, **kwargs):
        """WRITEME"""
        tt = Refraction()
        tt.setDataContainer(data)
        tt.showVA(ax=ax, t=t, **kwargs)

    def getOffset(self, data=None, full=False):
        """Return vector of offsets (in m) between shot and receiver."""

        if data is None:
            data = self.dataContainer

        if full:
            pos = data.sensorPositions()
            s, g = data('s'), data('g')
            nd = data.size()
            off = [pos[int(s[i])].distance(pos[int(g[i])]) for i in range(nd)]
            return np.absolute(off)
        else:
            px = pg.x(data.sensorPositions())
            gx = np.array([px[int(g)] for g in data("g")])
            sx = np.array([px[int(s)] for s in data("s")])
            return np.absolute(gx - sx)

    def getMidpoint(self, data=None):
        """Return vector of offsets (in m) between shot and receiver."""
        if data is None:
            data = self.dataContainer

        px = pg.x(data.sensorPositions())
        gx = np.array([px[int(g)] for g in data("g")])
        sx = np.array([px[int(s)] for s in data("s")])
        return (gx + sx) / 2

    def showVA(self, data=None, t=None, name='va', pseudosection=False,
               squeeze=True, full=True, ax=None):
        """Show apparent velocity as image plot.

        TODO showXXX commands need to return ax and cbar .. if there is one

        """
        if data is None:
            data = self.dataContainer

        if ax is None:
            fig, ax = plt.subplots()
            self.figs[name] = fig

        self.axs[name] = ax
        if t is None:
            t = data('t')

        px = pg.x(data.sensorPositions())
        gx = np.array([px[int(g)] for g in data("g")])
        sx = np.array([px[int(s)] for s in data("s")])
        offset = self.getOffset(data=data, full=full)
        va = offset / t

        if pseudosection:
            midpoint = (gx + sx) / 2
            plotVecMatrix(midpoint, offset, va, squeeze=True, ax=ax,
                          label='Apparent slowness [s/m]')
        else:
            plotVecMatrix(gx, sx, va, squeeze=squeeze, ax=ax,
                          label='Apparent velocity [m/s]')
        return va

    def getDepth(self):
        """return a (a-priori guessed) depth of investigation"""
        return max(self.getOffset()) / 3.0  # rule of thumb

    def rayCoverage(self):
        """return ray coverage"""
        one = pg.RVector(self.dataContainer.size(), 1.)
        return self.fop.jacobian().transMult(one)

    def standardizedCoverage(self):
        """return standardized coverage vector (0|1) using neighbor info"""
        coverage = self.rayCoverage()
        C = self.fop.constraintsRef()
        return np.sign(np.absolute(C.transMult(C * coverage)))

    def showCoverage(self, ax=None, name='coverage'):
        """shows the ray coverage in logscale"""
        if ax is None:
            fig, ax = plt.subplots()
            self.figs[name] = fig

        self.axs[name] = ax
        cov = self.rayCoverage()
        pg.show(self.mesh, pg.log10(cov+min(cov[cov > 0])*.5), ax=ax,
                coverage=self.standardizedCoverage())

    def showModel(self, ax=None, vals=None, **kwargs):
        """WRITEME"""
        self.showResult(ax=ax, val=vals, **kwargs)

    def showResult(self, val=None, ax=None, cMin=None, cMax=None,
                   logScale=False, name='result', **kwargs):
        """show resulting velocity vector"""
        mesh = self.paraDomain()
        if val is None:
            val = self.velocity
        if cMin is None or cMax is None:
            cMin, cMax = interperc(val, 3)
        if ax is None:
            fig, ax = plt.subplots()
            self.figs[name] = fig
            ax, cbar = pg.show(mesh, val, logScale=logScale, ax=ax,
                               colorBar=True, cMin=cMin, cMax=cMax,
                               coverage=self.standardizedCoverage(), **kwargs)
            self.figs[name] = plt.gcf()
        else:
            gci = drawModel(ax, mesh, val, logScale=logScale,
                            colorBar=True, cMin=cMin, cMax=cMax,
                            coverage=self.standardizedCoverage(), **kwargs)
            labels = ['cMin', 'cMax', 'nLevs', 'orientation', 'label']
            subkwargs = {key: kwargs[key] for key in labels if key in kwargs}
            cbar = createColorBar(gci, **subkwargs)
        browser = CellBrowser(self.mesh, val, ax)
        browser.connect()

        self.axs[name] = ax
        if 'lines' in kwargs:
            plotLines(ax, kwargs['lines'])
        return ax, cbar

    def showResultAndFit(self, name='resultfit', **kwargs):
        """show two vertical subplots with result and data (with response)"""
        fig, axs = plt.subplots(nrows=2)
        self.figs[name] = fig
        self.showResult(ax=axs[0], **kwargs)
        self.showData(ax=axs[1], response=self.response)

    def saveFigures(self, name=None, ext='pdf'):
        """save all existing figures to files"""
        if name is None:
            name = self.basename
        if name is None or not any(name):
            name = 'out'
        for key in self.figs:
            self.figs[key].savefig(name+'-'+key+'.'+ext, bbox_inches='tight')

    def saveResult(self, folder=None, size=(16, 10), **kwargs):
        """Save the results in the specified folder.

        Saved items are:
            Inverted profile
            Velocity vector
            Coverage vector
            Standardized coverage vector
            Mesh (bms and vtk with results)
        """
        # TODO: How to extract the chi2 etc. from each iteration???

        subfolder = '/' + self.__class__.__name__
        path = getSavePath(folder, subfolder)

        if self.verbose:
            print('Saving refraction data to: {}'.format(path))

        np.savetxt(path + '/velocity.vector',
                   self.velocity)
        np.savetxt(path + '/velocity-cov.vector',
                   self.rayCoverage())
        np.savetxt(path + '/velocity-scov.vector',
                   self.standardizedCoverage())

        self.mesh.addExportData('Velocity', self.velocity)
        self.mesh.addExportData('Coverage', self.rayCoverage())
        self.mesh.addExportData('S_Coverage', self.standardizedCoverage())
        self.mesh.exportVTK(path + 'velocity')
        self.mesh.save(path + 'velocity-mesh')
        self.pd.save(path + 'velocity-pd')

        fig, ax = plt.subplots()
        self.showResult(ax=ax, cov=self.standardizedCoverage(), **kwargs)
        fig.set_size_inches(size)
        fig.savefig(path + '/velocity.pdf', bbox_inches='tight')

        return path, fig, ax


def test_Refraction():
    """Test Refraction manager stability some data/mesh set / data update"""
    import os
    datafile = os.path.dirname(__file__) + '/example_topo.sgt'

    ra = Refraction(datafile, verbose=False, doSave=False)
    ra.createMesh(depth=80)
    ra.inv.setMaxIter(1)

    ra.invert()
    m1 = ra.model()
    mesh = pg.Mesh(ra.mesh)

    ra.setMesh(mesh)
    ra.invert()
    m2 = ra.model()

    np.testing.assert_array_equal(m1, m2)

    ra.setData(pg.DataContainer(datafile, 's g'))
    m3 = ra.invert()

    np.testing.assert_array_equal(m1, m3)


def main():
    """Main"""
    parser = MethodManager.createArgParser(dataSuffix='sgt')
    options = parser.parse_args()

    ra = Refraction(verbose=not options.quiet, debug=pg.debug())

    ra.loadData(options.dataFileName)

    ra.showData()
    ra.showVA()
    ra.createMesh(depth=options.depth)
    ra.showMesh()
    ra.invert(lam=options.lam, max_iter=options.maxIter,
              robustData=options.robustData, blockyModel=options.blockyModel)
    ra.showResult()

if __name__ == '__main__':
    main()
    pg.wait()
