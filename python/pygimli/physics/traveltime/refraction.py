#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
    Class for managing seismic refraction data and doing inversions
'''

# general modules to import according to standards
import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
from pygimli.meshtools import createParaDomain2D, createMesh
from pygimli.mplviewer import drawModel, drawMesh, CellBrowser, createColorbar


def plotFirstPicks(ax, data, plotva=False):
    ''' plot first arrivals as lines '''
    px = pg.x(data.sensorPositions())
    gx = np.array([px[int(g)] for g in data("g")])
    sx = np.array([px[int(s)] for s in data("s")])
    tt = np.array(data("t"))
    if plotva:
        tt = np.absolute(gx-sx) / tt

    uns = np.unique(sx)
    cols = 'brgcmyk'
    for i, si in enumerate(uns):
        ti = tt[sx == si]
        gi = gx[sx == si]
        offset = gi - si
        print(si, len(ti), min(offset), max(offset))
        ii = gi.argsort()
        ax.plot(gi[ii], ti[ii], 'x-', color=cols[i % 7])
        ax.plot(si, 0., 's', color=cols[i % 7], markersize=8)

    ax.grid(True)


def showVA(ax, data):
    ''' show apparent velocity as image plot '''
    px = pg.x(data.sensorPositions())
    gx = np.asarray([px[int(g)] for g in data("g")])
    sx = np.asarray([px[int(s)] for s in data("s")])
    va = np.absolute(gx-sx) / data('t')
    A = np.ones((data.sensorCount(), data.sensorCount())) * np.nan
    for i in range(data.size()):
        A[int(data('s')[i]), int(data('g')[i])] = va[i]

    ax.imshow(A, interpolation='nearest')
    ax.grid(True)
    xt = np.arange(0, data.sensorCount(), 50)
    ax.set_xticks(xt)
    ax.set_xticklabels([str(px[xti]) for xti in xt])
    return va


def createGradientModel2D(data, mesh, VTop, VBot):
    """
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
    d = np.array([np.abs(np.dot(pos[i, :], n) - p[1])/nLen
                  for i in range(pos.shape[0])])

    return np.interp(d, [min(d), max(d)], [1.0/VTop, 1.0/VBot])


# Data handling class (logics)
class Refraction():
    '''
    class for managing a refraction seismics
    '''
    def __init__(self, filename=None, verbose=True, **kwargs):
        ''' init function with optional data load '''
        if filename is not None:
            self.load(filename)

    def __repr__(self):  # for print function
        out = ""
        return out

    def load(self, filename):
        ''' load data from file '''
        # check for file formats and import if necessary
        self.data = pg.DataContainer(filename, 's g')
        print(self.data)

    def showData(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()

        plotFirstPicks(ax, self.data)
        plt.show(block=False)

    def showVA(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()

        showVA(ax, self.data)
        plt.show(block=False)

    def getOffset(self):
        px = pg.x(self.data.sensorPositions())
        gx = np.array([px[int(g)] for g in self.data("g")])
        sx = np.array([px[int(s)] for s in self.data("s")])
        return np.absolute(gx-sx)

    def makeMesh(self, depth=None, quality=34.3):
        if depth is None:
            depth = max(self.getOffset())/3.
        self.poly = createParaDomain2D(self.data.sensorPositions(),
                                       paradepth=depth,
                                       paraboundary=0, boundary=0)
        self.mesh = createMesh(self.poly, quality=quality, smooth=(1, 10))
        self.mesh.createNeighbourInfos()

    def showMesh(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()

        drawMesh(ax, self.mesh)
        plt.show(block=False)

    def createFOP(self):  # type Dijkstra
        if not hasattr(self,'mesh'): # self.mesh is None:
            self.makeMesh()
        self.f = pg.TravelTimeDijkstraModelling(self.mesh, self.data, True)
        self.f.regionManager().setConstraintType(1)
        self.f.createRefinedForwardMesh(True)
        self.pd = self.f.regionManager().paraDomain() # no idea why needed

    def estimateError(self, absoluteError=0.001, relativeError=0.001):
        if relativeError > 1: # obviously in %
            relativeError /= 100.
        self.error = absoluteError / self.data('t') + relativeError
        if hasattr(self, 'INV'):  # self.INV is not None:
            self.INV.setRelativeError(self.error)

    def createInv(self):
        if not hasattr(self,'f'):
            self.createFOP()
        if not hasattr(self, 'error'):
            self.estimateError()

        self.tD = pg.RTrans()
        self.tM = pg.RTransLog()
        self.INV = pg.RInversion(self.data('t'), self.f, True, True)
        self.INV.setTransData(self.tD)
        self.INV.setTransModel(self.tM)
        self.INV.setRelativeError(self.error)

    def run(self, vtop=500., vbottom=5000., zweight=0.3, lam=30.):
        if not hasattr(self,'INV'):  # self.f is None:
            self.createInv()

        self.start = createGradientModel2D(ra.data, ra.mesh, vtop, vbottom)
        self.f.setStartModel(self.start)
        self.f.regionManager().setZWeight(zweight)
        self.INV.setLambda(lam)
        self.INV.setModel(self.start)
        slowness = self.INV.run()
        self.velocity = 1. / slowness

    def showResult(self, ax=None, **kwargs):
        if ax is None:
            ax, cbar = pg.show(self.mesh, self.velocity)
            # fig, ax = plt.subplots()
        else:
            gci = drawModel(ax, self.mesh, self.velocity, **kwargs)
            cbar = createColorbar(gci, *args, **kwargs)
            browser = CellBrowser(self.mesh, self.velocity, ax)
            browser.connect()
            plt.show() # block=False)


# # # # # # # # MAIN # # # # # # # #
if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser("usage: %prog [options] ",
                          version="%prog: " + pg.__version__)
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
                      help="be verbose", default=False)
    (options, args) = parser.parse_args()

    if options.verbose:
        __verbose__ = True

    if len(args) == 0:
        parser.print_help()
        print("Use with a data file name.")
        raise SystemExit
    else:
        datafile = args[0]

    ra = Refraction(datafile)
    print(ra)
    pg.showLater()
    ra.showData()
    ra.showVA()
    ra.makeMesh()
    ra.mesh.save('mesh')
    print(ra.mesh)
    ra.showMesh()
    
    ra.createFOP()
    
    slo = createGradientModel2D(ra.data, ra.mesh, 500., 5000.)
    
    #print('#'*100)
    resp = ra.f(slo)
    #print('#'*100)
    #ra.velocity = 1. / slo
    #ra.createInv()
    #ra.run()
    #ra.showResult()

    pg.showNow()
