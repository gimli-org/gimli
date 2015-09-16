#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Class for managing seismic refraction data and doing inversions"""

import sys
from math import pi
import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
from pygimli.meshtools import createParaDomain2D, createMesh
from pygimli.mplviewer import drawModel, drawMesh, CellBrowser, createColorbar
from pygimli.utils.base import interperc
from .ratools import createGradientModel2D, getSavePath
from .raplot import plotFirstPicks, showVA, plotLines


class Refraction(object):

    """ Class for managing refraction seismics data"""

    def __init__(self, data=None, verbose=True, **kwargs):
        """Init function with optional data load"""
        self.figs = {}
        self.axs = {}
        if isinstance(data, str):
            self.load(data)
        elif isinstance(data, pg.DataContainer):
            self.setData(data)
            self.basename = 'new'

    def __str__(self):
        """string representation of the class"""
        return self.__repr__()

    def __repr__(self):
        """string representation of the class for the print function"""
        out = "Refraction object"
        if hasattr(self, 'data'):
            out += "\n" + self.data.__str__()
        if hasattr(self, 'mesh'):
            out += "\n" + self.mesh.__str__()
        return out

    def setData(self, data):
        """ set data container from outside """
        self.data = data
        self.checkData()

    def load(self, filename):
        """load data from file"""
        # check for file formats and import if necessary
        self.data = pg.DataContainer(filename, 's g')
        self.basename = filename[:filename.rfind('.')]
        self.checkData()

    def checkData(self):
        """ check data w.r.t. shot/geophone identity and zero/negative
        traveltimes, plus check y/z sensor positions """
        oldsize = self.data.size()
        self.data.markInvalid(pg.abs(self.data('s') - self.data('g')) < 1)
        self.data.markInvalid(self.data('t') <= 0.)
        self.data.removeInvalid()
        newsize = self.data.size()
        if newsize < oldsize:
            print('Removed ' + str(oldsize - newsize) + ' values.')
        maxyabs = max(pg.abs(pg.y(self.data.sensorPositions())))
        maxzabs = max(pg.abs(pg.z(self.data.sensorPositions())))
        if maxzabs > 0 and maxyabs == 0:
            for i in range(self.data.sensorCount()):
                pos = self.data.sensorPosition(i).rotateX(-pi / 2)
                self.data.setSensorPosition(i, pos)

        print(self.data)

    def showData(self, ax=None, response=None):
        """show data in form of travel time curves (optional with response)"""
        if ax is None:
            fig, ax = plt.subplots()
            self.figs['data'] = fig

        self.axs['data'] = ax
        if response is None:
            plotFirstPicks(ax, self.data)
        else:
            plotFirstPicks(ax, self.data, marker='x')
            if response is True:
                response = self.response
            plotFirstPicks(ax, self.data, np.asarray(response), marker='-')

        plt.show(block=False)

    def showVA(self, ax=None):
        """show apparent velocity as image plot"""
        if ax is None:
            fig, ax = plt.subplots()
            self.figs['va'] = fig

        self.axs['va'] = ax
        showVA(ax, self.data)
        plt.show(block=False)

    def getOffset(self):
        """return vector of offsets (in m) between shot and receiver"""
        px = pg.x(self.data.sensorPositions())
        gx = np.array([px[int(g)] for g in self.data("g")])
        sx = np.array([px[int(s)] for s in self.data("s")])
        return np.absolute(gx - sx)

    def makeMesh(self, depth=None, quality=34.3, paraDX=0.5, boundary=0,
                 paraBoundary=5):
        """create (inversion) mesh using createParaDomain2D"""
        if depth is None:
            depth = max(self.getOffset()) / 3.
        self.poly = createParaDomain2D(self.data.sensorPositions(),
                                       paraDepth=depth, paraDX=paraDX,
                                       paraBoundary=paraBoundary,
                                       boundary=boundary)
        self.mesh = createMesh(self.poly, quality=quality, smooth=(1, 10))
        self.mesh.createNeighbourInfos()

    def showMesh(self, ax=None):
        """show mesh in given axes or in a new figure"""
        if ax is None:
            fig, ax = plt.subplots()
            self.figs['mesh'] = fig

        self.axs['mesh'] = ax
        drawMesh(ax, self.mesh)
        plt.show(block=False)

    def createFOP(self, refine=True):  # Dijkstra, later FMM
        """create forward operator working on refined mesh"""
        if not hasattr(self, 'mesh'):  # self.mesh is None:
            self.makeMesh()
        self.f = pg.TravelTimeDijkstraModelling(self.mesh, self.data, True)
        self.f.regionManager().setConstraintType(1)
        self.f.createRefinedForwardMesh(refine)
        self.pd = self.f.regionManager().paraDomain()  # no idea why needed

    def estimateError(self, absoluteError=0.001, relativeError=0.001):
        """estimate error composed of an absolute and a relative part"""
        if relativeError >= 0.5:  # obviously in %
            relativeError /= 100.
        self.error = absoluteError / self.data('t') + relativeError
        if hasattr(self, 'INV'):  # self.INV is not None:
            self.INV.setRelativeError(self.error)

    def createInv(self, verbose=True, dosave=False):
        """create inversion instance (typically done by run)"""
        if not hasattr(self, 'f'):
            self.createFOP()
        if not hasattr(self, 'error'):
            self.estimateError()

        self.tD = pg.RTrans()
        self.tM = pg.RTransLog()
        self.INV = pg.RInversion(self.data('t'), self.f, True, False)
        self.INV.setTransData(self.tD)
        self.INV.setTransModel(self.tM)
        self.INV.setRelativeError(self.error)

    def run(self, vtop=500., vbottom=5000., zweight=0.2, lam=30.):
        """run actual inversion (first creating inversion object if not there)

        result/response is stored in the class attribute velocity/response

        Parameters
        ----------
        vtop, vbottom : float
            starting model velocities on top/bottom of the mesh

        lam : float
            regularization parameter describing the strength of smoothness

        zweight : float
            relative weight for purely vertical boundaries
        """
        if not hasattr(self, 'INV'):  # self.f is None:
            self.createInv()

        self.start = createGradientModel2D(self.data, self.mesh, vtop, vbottom)
        self.f.setStartModel(self.start)
        self.f.regionManager().setZWeight(zweight)
        self.INV.setLambda(lam)
        self.INV.setModel(self.start)
        slowness = self.INV.run()
        self.velocity = 1. / slowness
        self.response = self.INV.response()

    def rayCoverage(self):
        """return ray coverage"""
        return self.f.jacobian().transMult(pg.RVector(self.data.size(), 1.))

    def standardizedCoverage(self):
        """return standardized coverage vector (0|1) using neighbor info"""
        coverage = self.rayCoverage()
        C = self.f.constraintsRef()
        return np.sign(np.absolute(C.transMult(C * coverage)))

    def showResult(self, ax=None, cMin=None, cMax=None, logScale=False,
                   **kwargs):
        """show resulting velocity vector"""
        if cMin is None or cMax is None:
            cMin, cMax = interperc(self.velocity, 3)
        if ax is None:
            ax, cbar = pg.show(self.mesh, self.velocity, logScale=logScale,
                               colorBar=True, cMin=cMin, cMax=cMax, **kwargs)
            self.figs['result'] = plt.gcf()
        else:
            gci = drawModel(ax, self.mesh, self.velocity, logScale=logScale,
                            colorBar=True, cMin=cMin, cMax=cMax, **kwargs)
            createColorbar(gci, **kwargs)
        browser = CellBrowser(self.mesh, self.velocity, ax)
        browser.connect()
        plt.show()  # block=False)

        self.axs['result'] = ax
        if 'lines' in kwargs:
            plotLines(ax, kwargs['lines'])

    def saveFigures(self, name=None, ext='pdf'):
        """save all existing figures to files"""
        if name is None:
            name = self.basename
        if name is None or not any(name):
            name = 'out'
        for key in self.figs:
            self.figs[key].savefig(name+'-'+key+'.'+ext, bbox_inches='tight')

    def saveResult(self, folder=None, size=(16, 10), **kwargs):
        """
        Saves the results in the specified folder.

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


if __name__ == '__main__':
    datafile = sys.argv[1]
    ra = Refraction(datafile)
    print(ra)
    ra.showData()
    ra.showVA()
    ra.makeMesh(depth=100)
    ra.showMesh()
    ra.run()
    ax, cbar = ra.showResult()
#    pg.showNow()
