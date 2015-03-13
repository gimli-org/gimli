#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Class for managing seismic refraction data and doing inversions"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
from pygimli.meshtools import createParaDomain2D, createMesh
from pygimli.mplviewer import drawModel, drawMesh, CellBrowser, createColorbar
from pygimli.utils.base import interperc
from math import pi
import time
import os


def plotLines(ax, line_filename, step=1):
    xz = np.loadtxt(line_filename)
    n_points = xz.shape[0]
    if step == 2:
        for i in range(0, n_points, step):
            x = xz[i:i + step, 0]
            z = xz[i:i + step, 1]
            ax.plot(x, z, 'k-')
    if step == 1:
        ax.plot(xz[:, 0], xz[:, 1], 'k-')


def createResultFolder(subfolder):
    now = time.localtime()
    results = str(now.tm_year) + str(now.tm_mon).zfill(2) + \
        str(now.tm_mday).zfill(2) + '-' + str(now.tm_hour).zfill(2) + '.' + \
        str(now.tm_min).zfill(2)

    return createfolders(['./', results, subfolder])


def createfolders(foldername_list):
    """

    """
    path = ''

    for s in foldername_list:
        path = path + s + '/'

    try:
        os.makedirs(path)
    except OSError as e:
        if os.path.exists(path):
            print('Path "{}" already exists.'.format(path))
        else:
            print('Unable to create path "{}".'.format(path))
            raise(e)

    return path


def getSavePath(folder=None, subfolder=''):
    if folder is None:
        path = createResultFolder(subfolder)
    else:
        path = createfolders([folder, subfolder])
    return path


def plotFirstPicks(ax, data, tt=None, plotva=False, marker='x-'):
    """ plot first arrivals as lines """
    px = pg.x(data.sensorPositions())
    gx = np.array([px[int(g)] for g in data("g")])
    sx = np.array([px[int(s)] for s in data("s")])
    if tt is None:
        tt = np.array(data("t"))
    if plotva:
        tt = np.absolute(gx - sx) / tt

    uns = np.unique(sx)
    cols = 'brgcmyk'
    for i, si in enumerate(uns):
        ti = tt[sx == si]
        gi = gx[sx == si]
        ii = gi.argsort()
        ax.plot(gi[ii], ti[ii], marker, color=cols[i % 7])
        ax.plot(si, 0., 's', color=cols[i % 7], markersize=8)

    ax.grid(True)


def showVA(ax, data, usepos=True):
    """ show apparent velocity as image plot """
    px = pg.x(data.sensorPositions())
    gx = np.asarray([px[int(g)] for g in data("g")])
    sx = np.asarray([px[int(s)] for s in data("s")])
    va = np.absolute(gx - sx) / data('t')
    A = np.ones((data.sensorCount(), data.sensorCount())) * np.nan
    for i in range(data.size()):
        A[int(data('s')[i]), int(data('g')[i])] = va[i]

    gci = ax.imshow(A, interpolation='nearest')
    ax.grid(True)
    xt = np.arange(0, data.sensorCount(), 50)
    if usepos:
        ax.set_xticks(xt)
        ax.set_xticklabels([str(int(px[xti])) for xti in xt])
        ax.set_yticks(xt)
        ax.set_yticklabels([str(int(px[xti])) for xti in xt])

    plt.colorbar(gci, ax=ax)
    return va


def createGradientModel2D(data, mesh, VTop, VBot):
    """
    Create 2D velocity gradient model.

    Creates a smooth, linear, starting model that takes the slope
    of the topography into account. This is done by fitting a straight line
    and using the distance to that as the depth value.
    Known as "The Marcus method"


    Parameters
    ----------
    data : pygimli DataContainer
        The topography list is in here.
    mesh : pygimli.Mesh
        The parametric mesh used for the inversion
    VTop : float
        The velocity at the surface of the mesh
    VBot : float
        The velocity at the bottom of the mesh

    Returns
    -------
    model : pygimli Vector, length M
        A numpy array with slowness values that can be used to start
        the inversion.
    """

    p = np.polyfit(pg.x(data.sensorPositions()), pg.y(data.sensorPositions()),
                   deg=1)  # slope-intercept form
    n = np.asarray([-p[0], 1.0])  # normal vector
    nLen = np.sqrt(np.dot(n, n))

    x = pg.x(mesh.cellCenters())
    z = pg.y(mesh.cellCenters())
    pos = np.column_stack((x, z))
    d = np.array([np.abs(np.dot(pos[i, :], n) - p[1]) / nLen
                  for i in range(pos.shape[0])])

    return np.interp(d, [min(d), max(d)], [1.0 / VTop, 1.0 / VBot])


class Refraction(object):

    """ Class for managing a refraction seismics"""

    def __init__(self, data=None, verbose=True, **kwargs):
        """ Init function with optional data load """
        self.figs = {}
        self.axs = {}
        if isinstance(data, str):
            self.load(data)
        elif isinstance(data, pg.DataContainer):
            self.setData(data)
            self.basename = 'new'

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        """ string representation of the class """
        out = "Refraction object"
        if hasattr(self, 'data'):
            out += "\n" + self.data.__str__()
        if hasattr(self, 'mesh'):
            out += "\n" + self.mesh.__str__()
        return out

    def setData(self, data):
        self.data = data
        self.checkData()

    def load(self, filename):
        """" load data from file """
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
        """ show data in form of travel time curves """
        if ax is None:
            fig, ax = plt.subplots()
            self.figs['data'] = fig

        self.axs['data'] = ax
        if response is None:
            plotFirstPicks(ax, self.data)
        else:
            plotFirstPicks(ax, self.data, marker='x')
            plotFirstPicks(ax, self.data, np.asarray(response), marker='-')

        plt.show(block=False)

    def showVA(self, ax=None):
        """ show apparent velocity as image plot """
        if ax is None:
            fig, ax = plt.subplots()
            self.figs['va'] = fig

        self.axs['va'] = ax
        showVA(ax, self.data)
        plt.show(block=False)

    def getOffset(self):
        """ return vector of offsets (in m) between shot and receiver """
        px = pg.x(self.data.sensorPositions())
        gx = np.array([px[int(g)] for g in self.data("g")])
        sx = np.array([px[int(s)] for s in self.data("s")])
        return np.absolute(gx - sx)

    def makeMesh(self, depth=None, quality=34.3, paraDX=0.5, boundary=0,
                 paraBoundary=5):
        """ create (inversion) """
        if depth is None:
            depth = max(self.getOffset()) / 3.
        self.poly = createParaDomain2D(self.data.sensorPositions(),
                                       paraDepth=depth, paraDX=paraDX,
                                       paraBoundary=paraBoundary,
                                       boundary=boundary)
        self.mesh = createMesh(self.poly, quality=quality, smooth=(1, 10))
        self.mesh.createNeighbourInfos()

    def showMesh(self, ax=None):
        """ show mesh in given axes or in a new figure """
        if ax is None:
            fig, ax = plt.subplots()
            self.figs['mesh'] = fig

        self.axs['mesh'] = ax
        drawMesh(ax, self.mesh)
        plt.show(block=False)

    def createFOP(self, refine=True):  # Dijkstra, later FMM
        """ create forward operator working on refined mesh """
        if not hasattr(self, 'mesh'):  # self.mesh is None:
            self.makeMesh()
        self.f = pg.TravelTimeDijkstraModelling(self.mesh, self.data, True)
        self.f.regionManager().setConstraintType(1)
        self.f.createRefinedForwardMesh(refine)
        self.pd = self.f.regionManager().paraDomain()  # no idea why needed

    def estimateError(self, absoluteError=0.001, relativeError=0.001):
        """ estimate error composed of an absolute and a relative part """
        if relativeError > 1:  # obviously in %
            relativeError /= 100.
        self.error = absoluteError / self.data('t') + relativeError
        if hasattr(self, 'INV'):  # self.INV is not None:
            self.INV.setRelativeError(self.error)

    def createInv(self, verbose=True, dosave=False):
        """ create inversion instance """
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
        """ run inversion """
        if not hasattr(self, 'INV'):  # self.f is None:
            self.createInv()

        self.start = createGradientModel2D(self.data, self.mesh, vtop, vbottom)
        self.f.setStartModel(self.start)
        self.f.regionManager().setZWeight(zweight)
        self.INV.setLambda(lam)
        self.INV.setModel(self.start)
        slowness = self.INV.run()
        self.velocity = 1. / slowness

    def rayCoverage(self):
        """ return ray coverage """
        return self.f.jacobian().transMult(pg.RVector(self.data.size(), 1.))

    def standardizedCoverage(self):
        """ return standardized coverage vector (0|1) using neighbor info """
        coverage = self.rayCoverage()
        C = self.f.constraintsRef()
        return np.sign(np.absolute(C.transMult(C * coverage)))

    def showResult(self, ax=None, cMin=None, cMax=None, logScale=False,
                   **kwargs):
        """ show resulting velocity vector """
        if cMin is None or cMax is None:
            cMin, cMax = interperc(self.velocity, 3)
        if ax is None:
            ax, cbar = pg.show(self.mesh, self.velocity, logScale=logScale,
                               colorBar=True, cMin=cMin, cMax=cMax, **kwargs)
            fig, ax = plt.subplots()
            self.figs['result'] = fig
        else:
            gci = drawModel(ax, self.mesh, self.velocity, logScale=logScale,
                            colorBar=True, cMin=cMin, cMax=cMax, **kwargs)
            cbar = createColorbar(gci, **kwargs)
            browser = CellBrowser(self.mesh, self.velocity, ax)
            browser.connect()
            plt.show()  # block=False)

        self.axs['result'] = ax
        if 'lines' in kwargs:
            plotLines(ax, kwargs['lines'])

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
    pg.showLater(True)
    ra.showData()
    ra.showVA()
    ra.makeMesh(depth=100)
    ra.showMesh()
    ra.run()
    ax, cbar = ra.showResult()
#    pg.showNow()
