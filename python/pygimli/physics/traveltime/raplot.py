#!/usr/bin/env python
# -*- coding: utf-8 -*-
""""WRITEME"""

import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg


def drawFirstPicks(ax, data, tt=None, plotva=False, marker='x-'):
    """Naming convention. drawFOO(ax, ... )"""
    return plotFirstPicks(ax=ax, data=data, tt=tt,
                          plotva=plotva, marker=marker)


def drawVA(ax, data, usePos=True):
    """Naming convention. drawFOO(ax, ... )"""
    return showVA(ax=ax, data=data, usepos=usePos)


def drawTravelTimeData(ax, data, t=None):
    """
        Draw first arrival traveltime data into mpl ax a.
        data of type \ref DataContainer must contain sensorIdx 's' and 'g'
        and thus being numbered internally [0..n)
    """
    x = pg.x(data.sensorPositions())
    # z = pg.z(data.sensorPositions())

    shots = pg.unique(pg.sort(data('s')))
    geoph = pg.unique(pg.sort(data('g')))

    startOffsetIDX = 0

    if min(min(shots), min(geoph)) == 1:
        startOffsetIDX = 1

    tShow = data('t')
    if t is not None:
        tShow = t

    ax.set_xlim([min(x), max(x)])
    ax.set_ylim([max(tShow), -0.002])
    ax.figure.show()

    for shot in shots:
        gIdx = pg.find(data('s') == shot)
        sensorIdx = [int(i__ - startOffsetIDX) for i__ in data('g')[gIdx]]
        ax.plot(x[sensorIdx], tShow[gIdx], 'x-')

    yPixel = ax.transData.inverted().transform_point((1, 1))[1] - \
        ax.transData.inverted().transform_point((0, 0))[1]
    xPixel = ax.transData.inverted().transform_point((1, 1))[0] - \
        ax.transData.inverted().transform_point((0, 0))[0]

    # draw shot points
    ax.plot(x[[int(i__ - startOffsetIDX) for i__ in shots]],
              np.zeros(len(shots)) + 8. * yPixel, 'gv', markersize=8)

    # draw geophone points
    ax.plot(x[[int(i__ - startOffsetIDX) for i__ in geoph]],
              np.zeros(len(geoph)) + 3. * yPixel, 'r^', markersize=8)

    ax.grid()
    ax.set_ylim([max(tShow), +16. * yPixel])
    ax.set_xlim([min(x) - 5. * xPixel, max(x) + 5. * xPixel])

    ax.set_xlabel('x-Coordinate [m]')
    ax.set_ylabel('Traveltime [ms]')
# def drawTravelTimeData(...)


def plotFirstPicks(ax, data, tt=None, plotva=False, marker='x-'):
    """plot first arrivals as lines"""
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
    """show apparent velocity as image plot"""
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
