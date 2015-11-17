#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Some data related viewer
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import pygimli as pg


def generateMatrix(xvec, yvec, vals, squeeze=False):
    """ return data matrix from x/y and value vector """
    if squeeze:
        xmap = {xx: ii for ii, xx in enumerate(np.unique(xvec))}
        ymap = {yy: ii for ii, yy in enumerate(np.unique(yvec))}
    else:
        xymap = {xy: ii for ii, xy in enumerate(np.unique(np.hstack((xvec,
                                                                     yvec))))}
        xmap = xymap
        ymap = xymap
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


def plotMatrix(A, xmap=None, ymap=None, ax=None, cMin=None, cMax=None,
               logScale=None, label=None, **kwargs):
    """ plot previously generated matrix """
    if xmap is None:
        xmap = {i: i for i in range(A.shape[0])}
    if ymap is None:
        ymap = {i: i for i in range(A.shape[1])}
    mat = np.ma.masked_where(A == 0.0, A, False)
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
    if ax is None:
        fig, ax = plt.subplots()

    im = ax.imshow(mat, norm=norm, interpolation='nearest')
    pg.mplviewer.createColorbar(im, cMin=cMin, cMax=cMax, nLevs=5, label=label)
    ax.grid(True)
    xt = np.unique(ax.get_xticks().clip(0, len(xmap)-1))
    yt = np.unique(ax.get_xticks().clip(0, len(ymap)-1))
    xx = [k for k in xmap]
    ax.set_xticks(xt)
    ax.set_xticklabels(['{:g}'.format(round(xx[int(ti)], 2)) for ti in xt])
    yy = [k for k in ymap]
    ax.set_yticks(yt)
    ax.set_yticklabels(['{:g}'.format(round(yy[int(ti)], 2)) for ti in yt])
    return ax


def plotVecMatrix(xvec, yvec, vals, squeeze=True, **kwargs):
    """ plot three vectors as matrix """
    A, xmap, ymap = generateMatrix(xvec, yvec, vals, squeeze)
    return plotMatrix(A, xmap, ymap, **kwargs)


def plotDataContainerAsMatrix(data, x=None, y=None, v=None, **kwargs):
    """ plot data container as matrix

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
    return plotVecMatrix(x, y, v, **kwargs)


def drawSensorAsMarker(ax, data):
    """
        Draw Sensor marker, these marker are pickable
    """
    elecsX = []
    elecsY = []

    for i in range(len(data.sensorPositions())):
        elecsX.append(data.sensorPositions()[i][0])
        elecsY.append(data.sensorPositions()[i][1])

    electrodeMarker, =  ax.plot(elecsX, elecsY, 'x', color='black', picker=5.)

    ax.set_xlim([data.sensorPositions()[0][0] - 1.,
                 data.sensorPositions()[data.sensorCount() - 1][0] + 1.])

    return electrodeMarker


def drawTravelTimeData(a, data):
    """
        Draw first arrival traveltime data into mpl axes a.
        data of type \ref DataContainer must contain sensorIdx 's' and 'g'
        and thus being numbered internally [0..n)
    """

    x = pg.x(data.sensorPositions())
#    z = pg.z(data.sensorPositions())

    shots = pg.unique(pg.sort(data('s')))
    geoph = pg.unique(pg.sort(data('g')))

    startOffsetIDX = 0

    if min(min(shots), min(geoph)) == 1:
        startOffsetIDX = 1

    a.set_xlim([min(x), max(x)])
    a.set_ylim([max(data('t')), -0.002])
    a.figure.show()

    for shot in shots:
        gIdx = pg.find(data('s') == shot)
        sensorIdx = [int(i__ - startOffsetIDX) for i__ in data('g')[gIdx]]
        a.plot(x[sensorIdx], data('t')[gIdx], 'x-')

    yPixel = a.transData.inverted().transform_point((1, 1))[1] - \
        a.transData.inverted().transform_point((0, 0))[1]
    xPixel = a.transData.inverted().transform_point((1, 1))[0] - \
        a.transData.inverted().transform_point((0, 0))[0]

    # draw shot points
    a.plot(x[[int(i__ - startOffsetIDX) for i__ in shots]],
           np.zeros(len(shots)) + 8. * yPixel, 'gv', markersize=8)

    # draw geophone points
    a.plot(x[[int(i__ - startOffsetIDX) for i__ in geoph]],
           np.zeros(len(geoph)) + 3. * yPixel, 'r^', markersize=8)

    a.grid()
    a.set_ylim([max(data('t')), +16. * yPixel])
    a.set_xlim([min(x) - 5. * xPixel, max(x) + 5. * xPixel])

    a.set_xlabel('x-Coordinate [m]')
    a.set_ylabel('Traveltime [ms]')
# def drawTravelTimeData(...)
