import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg


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
