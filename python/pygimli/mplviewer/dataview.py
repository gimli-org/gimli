#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Some data related viewer."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Rectangle, Wedge

import pygimli as pg

from .utils import updateAxes as updateAxes_


def generateMatrix(xvec, yvec, vals, **kwargs):
    """Generate a data matrix from x/y and value vectors.

    Parameters
    ----------
    xvec, yvec, vals : iterables (list, np.array, pg.Vector) of same length

    full: bool [False]
        generate a fully symmetric matrix containing all unique xvec+yvec
        otherwise A is squeezed to the individual unique vectors

    Returns
    -------
    A : np.ndarray(2d)
        matrix containing the values sorted according to unique xvec/yvec
    xmap/ymap : dict {key: num}
        dictionaries for accessing matrix position (row/col number from x/y[i])
    """
    if kwargs.pop('full', False):
        xymap = {xy: ii
                 for ii, xy in enumerate(np.unique(np.hstack((xvec, yvec))))}
        xmap = xymap
        ymap = xymap
    else:
        xu, yu = np.unique(xvec), np.unique(yvec)
        if kwargs.pop('fillx', False):
            print('filling x', len(xu))
            dx = np.median(np.diff(xu)).round(1)
            xu = np.arange(0, xu[-1] - xu[0] + dx * 0.5, dx) + xu[0]
            print(len(xu))
        if kwargs.pop('filly', False):
            dy = np.median(np.diff(yu)).round(1)
            yu = np.arange(0, yu[-1] - yu[0] + dy * 0.5, dy) + yu[0]
        xmap = {xx: ii for ii, xx in enumerate(xu)}
        ymap = {yy: ii for ii, yy in enumerate(yu)}

    A = np.zeros((len(ymap), len(xmap)))
    inot = []
    nshow = min([len(xvec), len(yvec), len(vals)])

    for i in range(nshow):
        xi, yi = xvec[i], yvec[i]
        if A[ymap[yi], xmap[xi]]:
            inot.append(i)
        A[ymap[yi], xmap[xi]] = vals[i]

    if len(inot) > 0:
        print(len(inot), "data of", nshow, "not shown")
        if len(inot) < 30:
            print(inot)

    return A, xmap, ymap

def showValMapPatches(vals, xVec=None, yVec=None, dx=1, dy=None, **kwargs):
    """ """
    ax, _ = pg.show(ax=kwargs.pop('ax', None))

    gci, ymap = drawValMapPatches(ax, vals, xVec=xVec, yVec=yVec, dx=dx, dy=dy,
                            **kwargs)
    
    cbar = None
    if not kwargs.pop('colorBar', False):
        cbar = pg.mplviewer.createColorBar(gci, **kwargs)

    return ax, cbar, ymap

def drawValMapPatches(ax, vals, xVec=None, yVec=None, dx=1, dy=None, **kwargs):
    """ """
    recs = []

    circular = kwargs.pop('circular', False)
    if circular:
        recs = [None] * len(xVec)
        if dy is None:  # map y values to unique
            ymap = {xy: ii for ii, xy in enumerate(np.unique(yVec))}

            xyMap = {}
            for i, y in enumerate(yVec):
                if y not in xyMap:
                    xyMap[y] = []
                xyMap[y].append(i)

            # maxR = max(ymap.values())  # what's that for? not used
            dR = 1 / (len(ymap.values())+1)
            # dOff = np.pi / 2  # what's that for? not used

            for y, xIds in xyMap.items():
                r = 1. - dR*(ymap[y]+1)
                # ax.plot(r * np.cos(xvec[xIds]),
                #         r * np.sin(xvec[xIds]), 'o')

                # print(y, ymap[y])
                for i in xIds:
                    phi = xVec[i]
                    # x = r * np.cos(phi)  # what's that for? not used
                    y = r * np.sin(phi)

                    dPhi = (xVec[1] - xVec[0])

                    recs[i] = Wedge((0., 0.), r + dR/1.5,
                                    (phi - dPhi)*360/(2*np.pi),
                                    (phi + dPhi)*360/(2*np.pi),
                                    width=dR,
                                    zorder=1+r)
                    # if i < 5:
                    #     ax.text(x, y, str(i))
                    # pg.wait()
        else:
            raise("Implementme")
    else:
        if dy is None:  # map y values to unique
            ymap = {xy: ii for ii, xy in enumerate(np.unique(yVec))}
            for i in range(len(vals)):
                recs.append(Rectangle((xVec[i] - dx / 2, ymap[yVec[i]] - 0.5),
                                      dx, 1))
        else:
            for i in range(len(vals)):
                recs.append(Rectangle((xVec[i] - dx / 2, yVec[i] - dy / 2),
                                      dx, dy))
        ax.set_xlim(min(xVec) - dx / 2, max(xVec) + dx / 2)
        ax.set_ylim(len(ymap) - 0.5, -0.5)

    pp = PatchCollection(recs)
    pp.set_edgecolor(None)
    pp.set_linewidths(0.0)
    pp.set_array(vals)
    
    gci = ax.add_collection(pp)

    if circular:
        pp.set_edgecolor('black')
        pp.set_linewidths(0.1)

    return gci, ymap


def patchValMap(vals, xvec=None, yvec=None, ax=None, cMin=None, cMax=None,
                logScale=None, label=None, dx=1, dy=None, **kwargs):
    """Plot previously generated (generateVecMatrix) y map (category).

    Parameters
    ----------
    vals : iterable
        Data values to show.
    xvec : dict {i:num}
        dict (must match vals.shape[0])
    ymap : iterable
        vector for x axis (must match vals.shape[0])
    ax : mpl.axis
        axis to plot, if not given a new figure is created
    cMin/cMax : float
        minimum/maximum color values
    logScale : bool
        logarithmic colour scale [min(vals)>0]
    label : string
        colorbar label
    ** kwargs:
        * circular : bool
            Plot in polar coordinates.
    """
    if cMin is None:
        cMin = np.min(vals)
    if cMax is None:
        cMax = np.max(vals)

    if logScale is None:
        logScale = (cMin > 0.0)

    norm = None
    if logScale and cMin > 0:
        norm = LogNorm(vmin=cMin, vmax=cMax)
    else:
        norm = Normalize(vmin=cMin, vmax=cMax)

    if ax is None:
        ax = plt.subplots()[1]

    recs = []

    circular = kwargs.pop('circular', False)
    if circular:
        recs = [None] * len(xvec)
        if dy is None:  # map y values to unique
            ymap = {xy: ii for ii, xy in enumerate(np.unique(yvec))}

            xyMap = {}
            for i, y in enumerate(yvec):
                if y not in xyMap:
                    xyMap[y] = []
                xyMap[y].append(i)

            # maxR = max(ymap.values())  # what's that for? not used
            dR = 1 / (len(ymap.values())+1)
            # dOff = np.pi / 2  # what's that for? not used

            for y, xIds in xyMap.items():
                r = 1. - dR*(ymap[y]+1)
                # ax.plot(r * np.cos(xvec[xIds]),
                #         r * np.sin(xvec[xIds]), 'o')

                # print(y, ymap[y])
                for i in xIds:
                    phi = xvec[i]
                    # x = r * np.cos(phi)  # what's that for? not used
                    y = r * np.sin(phi)

                    dPhi = (xvec[1] - xvec[0])

                    recs[i] = Wedge((0., 0.), r + dR/1.5,
                                    (phi - dPhi)*360/(2*np.pi),
                                    (phi + dPhi)*360/(2*np.pi),
                                    width=dR,
                                    zorder=1+r)
                    # if i < 5:
                    #     ax.text(x, y, str(i))
                    # pg.wait()
        else:
            raise("Implementme")
    else:
        if dy is None:  # map y values to unique
            ymap = {xy: ii for ii, xy in enumerate(np.unique(yvec))}
            for i in range(len(vals)):
                recs.append(Rectangle((xvec[i] - dx / 2, ymap[yvec[i]] - 0.5),
                                      dx, 1))
        else:
            for i in range(len(vals)):
                recs.append(Rectangle((xvec[i] - dx / 2, yvec[i] - dy / 2),
                                      dx, dy))
        ax.set_xlim(min(xvec) - dx / 2, max(xvec) + dx / 2)
        ax.set_ylim(len(ymap) - 0.5, -0.5)

    pp = PatchCollection(recs)
    # ax.clear()
    col = ax.add_collection(pp)
    pp.set_edgecolor(None)
    pp.set_linewidths(0.0)

    if circular:
        pp.set_edgecolor('black')
        pp.set_linewidths(0.1)

    cmap = pg.mplviewer.cmapFromName(**kwargs)
    if kwargs.pop('markOutside', False):
        cmap.set_bad('grey')
        cmap.set_under('darkgrey')
        cmap.set_over('lightgrey')
        cmap.set_bad('black')
    pp.set_cmap(cmap)

    pp.set_norm(norm)
    pp.set_array(vals)
    pp.set_clim(cMin, cMax)

    updateAxes_(ax)
    cbar = kwargs.pop('colorBar', True)

    ori = kwargs.pop('orientation', 'horizontal')
    if cbar in ['horizontal', 'vertical']:
        ori = cbar
        cbar = True

    if cbar is True:  # not for cbar=1, which is really confusing!
        cbar = pg.mplviewer.createColorBar(col, cMin=cMin, cMax=cMax,
                                           nLevs=5, label=label,
                                           orientation=ori)

    elif cbar is not False:
        # .. cbar is an already existing cbar .. so we update its values
        pg.mplviewer.updateColorBar(cbar, cMin=cMin, cMax=cMax,
                                    nLevs=5, label=label)

    updateAxes_(ax)
    return ax, cbar, ymap


def patchMatrix(mat, xmap=None, ymap=None, ax=None, cMin=None, cMax=None,
                logScale=None, label=None, dx=1, **kwargs):
    """Plot previously generated (generateVecMatrix) matrix.

    Parameters
    ----------
    mat : numpy.array2d
        matrix to show
    xmap : dict {i:num}
        dict (must match A.shape[0])
    ymap : iterable
        vector for x axis (must match A.shape[0])
    ax : mpl.axis
        axis to plot, if not given a new figure is created
    cMin/cMax : float
        minimum/maximum color values
    logScale : bool
        logarithmic colour scale [min(A)>0]
    label : string
        colorbar label
    dx : float
        width of the matrix elements (by default 1)
    """
    mat = np.ma.masked_where(mat == 0.0, mat, False)
    if cMin is None:
        cMin = np.min(mat)
    if cMax is None:
        cMax = np.max(mat)
    if logScale is None:
        logScale = (cMin > 0.0)
    if logScale:
        norm = LogNorm(vmin=cMin, vmax=cMax)
    else:
        norm = Normalize(vmin=cMin, vmax=cMax)

    if 'ax' is None:
        ax = plt.subplots()[1]

    iy, ix = np.nonzero(mat)  # != 0)
    recs = []
    vals = []
    for i, _ in enumerate(ix):
        recs.append(Rectangle((ix[i] - dx / 2, iy[i] - 0.5), dx, 1))
        vals.append(mat[iy[i], ix[i]])

    pp = PatchCollection(recs)
    col = ax.add_collection(pp)
    pp.set_edgecolor(None)
    pp.set_linewidths(0.0)
    if 'cmap' in kwargs:
        pp.set_cmap(kwargs.pop('cmap'))
    if 'cMap' in kwargs:
        pp.set_cmap(kwargs.pop('cMap'))
    pp.set_norm(norm)
    pp.set_array(np.array(vals))
    pp.set_clim(cMin, cMax)
    xval = [k for k in xmap.keys()]
    ax.set_xlim(min(xval) - dx / 2, max(xval) + dx / 2)
    ax.set_ylim(len(ymap) + 0.5, -0.5)

    updateAxes_(ax)
    cbar = None
    if kwargs.pop('colorBar', True):
        ori = kwargs.pop('orientation', 'horizontal')
        cbar = pg.mplviewer.createColorBar(col, cMin=cMin, cMax=cMax, nLevs=5,
                                           label=label, orientation=ori)
    return ax, cbar


def plotMatrix(mat, *args, **kwargs):
    """Naming conventions. Use drawMatrix or showMatrix"""
    pg.deprecated("use drawMatrix or showMatrix")
    return showMatrix(*args, **kwargs)

def showMatrix(mat, xmap=None, ymap=None, **kwargs):
    """Show value map as matrix.

    Returns
    -------
    ax : matplotlib axes object
        axes object
    cb : matplotlib colorbar
        colorbar object
    """
    ax, _ = pg.show(ax=kwargs.pop('ax', None))

    # pg._r(ax)

    gci = drawMatrix(ax, mat, **kwargs)
    # pg._y(**kwargs)

    cbar = None
    if kwargs.pop('colorBar', True):
        ori = kwargs.pop('orientation', 'horizontal')
        cMin = kwargs.pop('cMin', None)
        cMax = kwargs.pop('cMax', None)
        label = kwargs.pop('label', None)
        cbar = pg.mplviewer.createColorBar(gci, cMin=cMin, cMax=cMax, nLevs=5,
                                           label=label, orientation=ori)

    return ax, cbar

def drawMatrix(ax, mat, xmap=None, ymap=None, cMin=None, cMax=None,
               logScale=None, label=None, **kwargs):
    """Draw previously generated (generateVecMatrix) matrix.

    Parameters
    ----------
    ax : mpl.axis
        axis to plot, if not given a new figure is created
    mat : numpy.array2d
        matrix to show
    xmap : dict {i:num}
        dict (must match A.shape[0])
    ymap : iterable
        vector for x axis (must match A.shape[0])
    cMin/cMax : float
        minimum/maximum color values
    logScale : bool
        logarithmic colour scale [min(A)>0]
    label : string
        colorbar label
    """
    if xmap is None:
        xmap = {i: i for i in range(mat.shape[0])}
    if ymap is None:
        ymap = {i: i for i in range(mat.shape[1])}
    if isinstance(mat, np.ma.MaskedArray):
        mat_ = mat
    else:
        mat_ = np.ma.masked_where(mat == 0.0, mat, False)

    if cMin is None:
        cMin = np.min(mat_)
    if cMax is None:
        cMax = np.max(mat_)
    if logScale is None:
        logScale = (cMin > 0.0)
    if logScale:
        norm = LogNorm(vmin=cMin, vmax=cMax)
    else:
        norm = Normalize(vmin=cMin, vmax=cMax)

    gci = ax.imshow(mat_, norm=norm, interpolation='nearest')
    if 'cmap' in kwargs:
        pg.deprecated('use cMap') #190422
        gci.set_cmap(kwargs.pop('cmap'))
    if 'cMap' in kwargs:
        gci.set_cmap(kwargs.pop('cMap'))
    ax.set_aspect(kwargs.pop('aspect', 1))

    ax.grid(True)
    xt = np.unique(ax.get_xticks().clip(0, len(xmap) - 1))
    yt = np.unique(ax.get_xticks().clip(0, len(ymap) - 1))

    if kwargs.pop('showally', False):
        yt = np.arange(len(ymap))
    else:
        yt = np.round(np.linspace(0, len(ymap) - 1, 5))

    xx = np.sort([k for k in xmap])
    ax.set_xticks(xt)
    ax.set_xticklabels(['{:g}'.format(round(xx[int(ti)], 2)) for ti in xt])

    yy = np.unique([k for k in ymap])
    ax.set_yticks(yt)
    ax.set_yticklabels(['{:g}'.format(round(yy[int(ti)], 2)) for ti in yt])
    return gci


def plotVecMatrix(xvec, yvec, vals, full=False, **kwargs):
    """DEPRECATED for nameing
    """
    pg.deprecated("use drawVecMatrix or showVecMatrix")
    return showVecMatrix(xvec, yvec, vals, full, **kwargs)


def showVecMatrix(xvec, yvec, vals, full=False, **kwargs):
    """Plot three vectors as matrix.

    Parameters
    ----------
    xvec, yvec : iterable (e.g. list, np.array, pg.Vector) of identical length
        vectors defining the indices into the matrix
    vals : iterable of same length as xvec/yvec
        vector containing the values to show
    full: bool [False]
        use a fully symmetric matrix containing all unique xvec+yvec
        otherwise A is squeezed to the individual unique xvec/yvec values

    **kwargs: forwarded to plotMatrix
        *   ax : mpl.axis
            Axis to plot, if not given a new figure is created
        * cMin/cMax : float
            Minimum/maximum color values
        * logScale : bool
            Lgarithmic colour scale [min(A)>0]
        * label : string
            Colorbar label

    Returns
    -------
    ax : matplotlib axes object
        axes object
    cb : matplotlib colorbar
        colorbar object
    """
    A, xmap, ymap = generateMatrix(xvec, yvec, vals, full=full)
    return showMatrix(A, xmap=xmap, ymap=ymap, **kwargs)

def drawVecMatrix(ax, xvec, yvec, vals, full=False, **kwargs):
    A, xmap, ymap = generateMatrix(xvec, yvec, vals, full=full)
    return drawMatrix(ax, A, xmap=xmap, ymap=ymap, **kwargs)


def plotDataContainerAsMatrix(data, x=None, y=None, v=None, **kwargs):
    """Plot data container as matrix.

    for each x, y and v token strings or vectors should be given
    """
    if isinstance(x, str):
        x = data(x)
    if isinstance(y, str):
        y = data(y)
    if isinstance(v, str):
        v = data(v)
    if x is None or y is None or v is None:
        raise Exception("Vectors or strings must be given")
    if len(x) != len(y) or len(x) != len(v):
        raise Exception("lengths x/y/v not matching: {:d}!={:d}!={:d}".format(
            len(x), len(y), len(v)))
    return showVecMatrix(x, y, v, **kwargs)


def drawSensorAsMarker(ax, data):
    """Draw Sensor marker, these marker are pickable."""
    elecsX = []
    elecsY = []

    for i in range(len(data.sensorPositions())):
        elecsX.append(data.sensorPositions()[i][0])
        elecsY.append(data.sensorPositions()[i][1])

    electrodeMarker, = ax.plot(elecsX, elecsY, 'x', color='black', picker=5.)

    ax.set_xlim([data.sensorPositions()[0][0] - 1.,
                 data.sensorPositions()[data.sensorCount() - 1][0] + 1.])

    return electrodeMarker
