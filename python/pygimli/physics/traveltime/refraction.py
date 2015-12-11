#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Class for managing seismic refraction data and doing inversions"""

import sys
from math import pi
import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
from pygimli.meshtools import createParaMeshPLC, createMesh
from pygimli.mplviewer import drawModel, drawMesh, CellBrowser, createColorbar
from pygimli.utils.base import interperc, getSavePath
from pygimli.mplviewer.dataview import plotVecMatrix

if __name__ == 'refraction':
    # keine Ahnung wie das eleganter geht
    # nur dafür:
    # python -c 'import refraction; refraction.test_Refraction()'

    from ratools import createGradientModel2D
    from raplot import plotFirstPicks, plotLines
else:
    from . ratools import createGradientModel2D
    from . raplot import plotFirstPicks, plotLines


class Refraction(object):  # to be derived from MethodManager(ND/Mesh)
    """ Class for managing refraction seismics data

        Document main members or use default MethodeManager interface
        e.g., self.inv, self.fop, self.paraDomain, self.mesh, self.data
    """

    def __init__(self, data=None, name='new', **kwargs):
        """Init function with optional data load"""
        self.figs = {}
        self.axs = {}

        self.verbose = kwargs.pop('verbose', True)
        self.doSave = kwargs.pop('doSave', False)

        # should be forwarded so it can be accessed from outside
        self.dataContainer = None
        self.mesh = None
        self.poly = None
        self.error = None
        self.velocity = None
        self.response = None
        self.start = []
        self.pd = None
        self.fop = self.createFop(verbose=self.verbose)
        self.inv = self.createInv(self.fop,
                                  verbose=self.verbose, doSave=self.doSave)

        if isinstance(data, str):
            self.loadData(data)
        elif isinstance(data, pg.DataContainer):
            self.setDataContainer(data)
            self.basename = name
        if self.dataContainer is not None:
            self.createMesh()

    def __str__(self):
        """string representation of the class"""
        return self.__repr__()

    def __repr__(self):
        """string representation of the class for the print function"""
        out = type(self).__name__ + " object"
        if hasattr(self, 'dataContainer'):
            out += "\n" + self.dataContainer.__str__()
        if hasattr(self, 'mesh'):
            out += "\n" + self.mesh.__str__()
        return out

    def paraDomain(self):
        """ base api """
        return self.fop.regionManager().paraDomain()

    def model(self):
        """ base api """
        return self.velocity

    def createFop(self, verbose=False):
        """Create default forward operator for Traveltime modelling.
        base api

        Dijkstra, later FMM.
        """
        if not hasattr(self, 'mesh'):  # self.mesh is None:
            self.makeMesh()
        fop = pg.TravelTimeDijkstraModelling(verbose=True)
        return fop

    def createInv(self, fop, verbose=True, doSave=False):
        """Create default inversion instance for Traveltime inversion.

        base api
        (typically done by run)
        """
        self.tD = pg.RTrans()
        self.tM = pg.RTransLog()

        inv = pg.RInversion(verbose, doSave)
        inv.setTransData(self.tD)
        inv.setTransModel(self.tM)
        inv.setForwardOperator(fop)

        return inv

    def setDataContainer(self, data):
        """ set data container from outside

        base api

        """
        self.dataContainer = data
        self.checkData()
        self.fop.setData(self.dataContainer)
        self.inv.setData(self.dataContainer('t'))
        if self.dataContainer.allNonZero('err'):
            self.error = self.dataContainer('err')
        else:
            self.estimateError()

    def loadData(self, filename):
        """load data from file"""
        # check for file formats and import if necessary
        data = pg.DataContainer(filename, 's g')
        self.basename = filename[:filename.rfind('.')]
        self.setDataContainer(data)

    def checkData(self):
        """ check data w.r.t. shot/geophone identity and zero/negative
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

    def showData(self, ax=None, response=None, name='data'):
        """show data as travel time curves (optionally with response)"""
        if response is not None:
            name = 'datafit'
        if ax is None:
            fig, ax = plt.subplots()
            self.figs[name] = fig

        self.axs[name] = ax
        if response is None:
            plotFirstPicks(ax, self.dataContainer)
        else:
            plotFirstPicks(ax, self.dataContainer, marker='+')
            if response is True:
                response = self.response
            plotFirstPicks(ax, self.dataContainer, np.asarray(response),
                           marker='-')

        plt.show(block=False)

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

    def createMesh(self, depth=None, quality=34.3, paraDX=0.5, boundary=0,
                   paraBoundary=5):
        """Create (inversion) mesh using createParaDomain2D"""
        if depth is None:
            depth = self.getDepth()
        self.poly = createParaMeshPLC(self.dataContainer.sensorPositions(),
                                      paraDepth=depth, paraDX=paraDX,
                                      paraBoundary=paraBoundary,
                                      boundary=boundary)
        mesh = createMesh(self.poly, quality=quality, smooth=(1, 10))
        self.setMesh(mesh)

    def showMesh(self, ax=None):
        """show mesh in given axes or in a new figure"""
        if ax is None:
            fig, ax = plt.subplots()
            self.figs['mesh'] = fig

        self.axs['mesh'] = ax
        drawMesh(ax, self.mesh)
        plt.show(block=False)

    def estimateError(self, absoluteError=0.001, relativeError=0.001):
        """estimate error composed of an absolute and a relative part"""
        if relativeError >= 0.5:  # obviously in %
            relativeError /= 100.
        self.error = absoluteError + self.dataContainer('t') * relativeError
        self.inv.setAbsoluteError(self.error)

    def createStartModel(self, vtop=500., vbottom=5000.):
        """create (gradient) starting model with vtop/vbottom bounds"""
        self.start = createGradientModel2D(self.dataContainer, self.mesh,
                                           vtop, vbottom)

    def run(self, **kwargs):
        """run actual inversion (first creating inversion object if not there)

        result/response is stored in the class attribute velocity/response

        Parameters
        ----------
        vtop, vbottom : float
            starting (gradient) model velocities on top/at bottom of the mesh

        lam : float
            regularization parameter describing the strength of smoothness

        zweight : float
            relative weight for purely vertical boundaries
        """
#        self.start = createGradientModel2D(self.dataContainer, self.mesh,
#                                           vtop, vbottom)
        if self.fop.regionManager().parameterCount() != len(self.start):
            self.createStartModel(kwargs.pop('vtop', 500.),
                                  kwargs.pop('vbottom', 5000.))
        self.fop.setStartModel(self.start)
        self.fop.regionManager().setZWeight(kwargs.pop('zweight', 0.2))
        self.inv.setData(self.dataContainer('t'))
        self.inv.setLambda(kwargs.pop('lam', 30.))
        self.inv.setMaxIter(kwargs.pop('max_iter', 20))
        self.inv.setRobustData(kwargs.pop('setRobust', False))
        if not hasattr(self.error, '__iter__'):
            self.estimateError(kwargs.pop('error', 0.003))  # abs. error in ms
        self.inv.setAbsoluteError(self.error)
        self.inv.setModel(self.start)  # somehow doubled with 3 lines above
        self.fop.jacobian().clear()
        slowness = self.inv.run()
        self.velocity = 1. / slowness
        self.response = self.inv.response()

    def getOffset(self):
        """return vector of offsets (in m) between shot and receiver"""
        px = pg.x(self.dataContainer.sensorPositions())
        gx = np.array([px[int(g)] for g in self.dataContainer("g")])
        sx = np.array([px[int(s)] for s in self.dataContainer("s")])
        return np.absolute(gx - sx)

    def getMidpoint(self):
        """return vector of offsets (in m) between shot and receiver"""
        px = pg.x(self.dataContainer.sensorPositions())
        gx = np.array([px[int(g)] for g in self.dataContainer("g")])
        sx = np.array([px[int(s)] for s in self.dataContainer("s")])
        return (gx + sx) / 2

    def showVA(self, ax=None, t=None, name='va', pseudosection=False,
               squeeze=True):
        """show apparent velocity as image plot"""
        if ax is None:
            fig, ax = plt.subplots()
            self.figs[name] = fig

        self.axs[name] = ax
        if t is None:
            t = self.dataContainer('t')
        px = pg.x(self.dataContainer.sensorPositions())
        gx = np.array([px[int(g)] for g in self.dataContainer("g")])
        sx = np.array([px[int(s)] for s in self.dataContainer("s")])
        offset = np.absolute(gx - sx)
        va = offset / t
        if pseudosection:
            midpoint = (gx + sx) / 2
            plotVecMatrix(midpoint, offset, va, squeeze=True)
        else:
            plotVecMatrix(gx, sx, va, squeeze=squeeze)
#        va = showVA(ax, self.dataContainer)
#        plt.show(block=False)
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

    def showCoverage(self, ax=None):
        """shows the ray coverage in logscale"""
        cov = self.rayCoverage()
        pg.show(self.mesh, pg.log10(cov+min(cov[cov > 0])*.5), axes=ax,
                coverage=self.standardizedCoverage())

    def showResult(self, val=None, ax=None, cMin=None, cMax=None,
                   logScale=False, name='result', **kwargs):
        """show resulting velocity vector"""
        if val is None:
            val = self.velocity
        if cMin is None or cMax is None:
            cMin, cMax = interperc(val, 3)
        if ax is None:
            ax, cbar = pg.show(self.mesh, val, logScale=logScale,
                               colorBar=True, cMin=cMin, cMax=cMax,
                               coverage=self.standardizedCoverage(), **kwargs)
            self.figs[name] = plt.gcf()
        else:
            gci = drawModel(ax, self.mesh, val, logScale=logScale,
                            colorBar=True, cMin=cMin, cMax=cMax,
                            coverage=self.standardizedCoverage(), **kwargs)
            labels = ['cMin', 'cMax', 'nLevs', 'orientation', 'label']
            subkwargs = {key: kwargs[key] for key in labels if key in kwargs}
            cbar = createColorbar(gci, **subkwargs)
        browser = CellBrowser(self.mesh, val, ax)
        browser.connect()

        self.axs[name] = ax
        if 'lines' in kwargs:
            plotLines(ax, kwargs['lines'])
        return ax, cbar

    def showResultAndFit(self, **kwargs):
        """show two vertical subplots with result and data (with response)"""
        fig, ax = plt.subplots(nrows=2)
        self.figs['resultfit'] = fig
        self.showResult(ax=ax[0], **kwargs)
        self.showData(ax=ax[1], response=self.response)

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
#        TODO: How to extract the chi2 etc. from each iteration???

        subfolder = '/' + self.__class__.__name__
        path = getSavePath(folder, subfolder)

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
        fig.savefig(path + '/velocity.pdf')

        return path, fig, ax


def test_Refraction():

    import os
    datafile = os.path.dirname(__file__) + '/example_topo.sgt'

    ra = Refraction(datafile, verbose=False, doSave=False)
    ra.createMesh(depth=80)
    ra.inv.setMaxIter(1)

    ra.start()
    m1 = ra.model()
    mesh = pg.Mesh(ra.mesh)

    ra.setMesh(mesh)
    ra.start()
    m2 = ra.model()

    np.testing.assert_array_equal(m1, m2)

    ra.setData(pg.DataContainer(datafile, 's g'))
    m3 = ra.start()

    np.testing.assert_array_equal(m1, m3)

def main(argv):
    """
    """
    implementme()


if __name__ == '__main__':
    datafile = sys.argv[1]
    ra = Refraction(datafile)
    print(ra)
    ra.showData()
    ra.showVA()
    ra.createMesh(depth=100)
    ra.showMesh()
    ra.run()
    ax, cbar = ra.showResult()
#    pg.showNow()
