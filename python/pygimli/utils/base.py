# -*- coding: utf-8 -*-
"""
pygimli base functions
"""
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

import pygimli as pg
from pygimli.mplviewer.colorbar import setMappableData


def gmat2numpy(mat):
    """convert pygimli matrix into numpy.array.

    TODO implement correct rval
    """
    nmat = np.zeros((len(mat), len(mat[0])))
    for i, row in enumerate(mat):
        nmat[i] = row
    return nmat


def numpy2gmat(nmat):
    """convert numpy.array into pygimli RMatrix.

    TODO implement correct rval
    """
    gmat = pg.RMatrix()
    for arr in nmat:
        gmat.push_back(arr)  # pg.asvector(arr))
    return gmat


def rndig(a, ndig=3):
    """round float using a number of counting digits."""
    if np.abs(a) < 1e-4:
        return a
    else:
        return np.around(a, ndig - int(np.ceil(np.log10(np.abs(a) + 1e-4))))


def num2str(a, fmtstr='%g'):
    """List of strings (deprecated, for backward-compatibility) """
    return [fmtstr % rndig(ai) for ai in a]


def inthist(a, vals, bins=None, islog=False):
    """What is this good for?"""
    if bins is None:
        bins = np.min((np.round(len(a) / 20), 10))

    if islog:
        hists, edges = np.histogram(np.log(a), bins=bins)
    else:
        hists, edges = np.histogram(a, bins=bins)

    cums = np.cumsum(np.hstack((0., hists))) / np.sum(hists) * 100.
    out = np.interp(vals, cums, edges)
    if islog:
        return np.exp(out)
    else:
        return out


def interperc(a, trimval=3.0, islog=False, bins=None):
    """What is this good for?"""
    return inthist(
        a, np.array([trimval, 100. - trimval]), bins=bins, islog=islog)


def interpExtrap(x, xp, yp):
    """np.interp function with linear extrapolation"""
    y = np.interp(x, xp, yp)
    y = np.where(x < xp[0], yp[0]+(x-xp[0])*(yp[0]-yp[1])/(xp[0]-xp[1]), y)
    return np.where(x > xp[-1], yp[-1]+(x-xp[-1])*(yp[-1]-yp[-2]) /
                    (xp[-1]-xp[-2]), y)


def showmymatrix(mat, x, y, dx=2, dy=1, xlab=None, ylab=None, cbar=None):
    """What is this good for?"""
    plt.imshow(mat, interpolation='nearest')
    plt.xticks(np.arange(0, len(x), dx), ["%g" % rndig(xi, 2) for xi in x])
    plt.yticks(np.arange(0, len(y), dy), ["%g" % rndig(yi, 2) for yi in y])
    plt.ylim((len(y) - 0.5, -0.5))
    if xlab is not None:
        plt.xlabel(xlab)
    if ylab is not None:
        plt.ylabel(ylab)
    plt.axis('auto')
    if cbar is not None:
        plt.colorbar(orientation=cbar)
    return


def draw1dmodel(x, thk=None, xlab=None, zlab="z in m", islog=True, fs=14, z0=0,
                **kwargs):
    """draw 1d block model defined by value and thickness vectors."""
#    if xlab is None:
#        xlab = "$\\rho$ in $\\Omega$m"

    if thk is None:  # gimli blockmodel (thk+x together) given
        nl = int(np.floor((len(x) - 1) / 2.)) + 1
        thk = np.asarray(x)[:nl - 1]
        x = np.asarray(x)[nl - 1:nl * 2 - 1]

    z1 = np.concatenate(([0], np.cumsum(thk))) + z0
    z = np.concatenate((z1, [z1[-1] * 1.2]))
    nl = len(x)  # x.size()
    px = np.zeros((nl * 2, 1))
    pz = np.zeros((nl * 2, 1))
    for i in range(nl):
        px[2 * i] = x[i]
        px[2 * i + 1] = x[i]
        pz[2 * i + 1] = z[i + 1]
        if i < nl - 1:
            pz[2 * i + 2] = z[i + 1]

#    plt.cla()
    li = []
    if islog:
        li = plt.semilogx(px, pz, **kwargs)
    else:
        li = plt.plot(px, pz, **kwargs)

    plt.gca().xaxis.set_label_position('top')

    locs = plt.xticks()[0]
    if len(locs) < 2:
        locs = np.hstack((min(x), locs, max(x)))
    elif len(locs) < 5:
        locs[0] = max(locs[0], min(x))
        locs[-1] = min(locs[-1], max(x))

    a = []
    for l in locs:
        a.append('%g' % rndig(l))

    plt.xticks(locs, a, fontsize=fs)
    plt.yticks(fontsize=fs)

    plt.xlim((np.min(x) * 0.9, np.max(x) * 1.1))
    plt.ylim((max(z1) * 1.15, 0.))
    if xlab is not None:
        plt.xlabel(xlab, fontsize=fs)
    if zlab is not None:
        plt.ylabel(zlab, fontsize=fs)
    plt.grid(which='both')
    # plt.show()
    return li


def draw1dmodelErr(x, xL, xU=None, thk=None, xcol='g', ycol='r', **kwargs):
    """TODO"""
    if thk is None:
        nlay = (len(x) + 1) / 2
        thk = np.array(x)[:nlay - 1]
        x = np.asarray(x)[nlay - 1:nlay * 2 - 1]
        thkL = np.array(xL)[:nlay - 1]
        thkU = np.array(xU)[:nlay - 1]
        xL = np.asarray(xL)[nlay - 1:nlay * 2 - 1]
        xU = np.asarray(xU)[nlay - 1:nlay * 2 - 1]

#    thk0 = np.hstack((thk, 0.))
#    thkL0 = np.hstack((thkL, 0.))
#    thkU0 = np.hstack((thkU, 0.))
    zm = np.hstack((np.cumsum(thk) - thk / 2, np.sum(thk) * 1.2))  # midpoint
    zc = np.cumsum(thk)  # cumulative
    draw1dmodel(x, thk, **kwargs)
    plt.xlim(min(xL) * 0.95, max(xU) * 1.05)
    plt.ylim(zm[-1] * 1.1, 0.)
    plt.errorbar(
        x, zm, fmt='.', xerr=np.vstack(
            (x - xL, xU - x)), ecolor=xcol, **kwargs)
    plt.errorbar((x[:-1] + x[1:]) / 2, zc, fmt='.',
                 yerr=np.vstack((thk - thkL, thkU - thk)), ecolor=ycol)


def draw1dmodelLU(x, xL, xU, thk=None, **kwargs):
    """draw 1d model with lower and upper bounds"""
    draw1dmodel(x, thk, color='red', **kwargs)
    for i in range(len(x)):
        x1 = np.array(x)
        x1[i] = xL[i]
        draw1dmodel(x1, thk, color='blue')
        x1[i] = xU[i]
        draw1dmodel(x1, thk, color='blue')

    li = draw1dmodel(x, thk, color='red', **kwargs)
    plt.xlim((min(xL) * 0.9, max(xU) * 1.1))
    return li


def showStitchedModels(models, ax=None, x=None, cmin=None, cmax=None,
                       islog=True, title=None, cmap='jet'):
    """show several 1d block models as (stitched) section"""
    if x is None:
        x = np.arange(len(models))

    nlay = int(np.floor((len(models[0]) + 1) / 2.))

    fig = None
    if ax is None:
        fig, ax = plt.subplots()

    dxmed2 = np.median(np.diff(x)) / 2.
    vals = np.zeros((len(models), nlay))
    patches = []
    maxz = 0.
    for i, imod in enumerate(models):
        if isinstance(imod, pg.RVector):
            vals[i, :] = imod(nlay - 1, 2 * nlay - 1)
            thk = np.asarray(imod(0, nlay - 1))
        else:
            vals[i, :] = imod[nlay - 1:2 * nlay - 1]
            thk = imod[:nlay - 1]

        thk = np.hstack((thk, thk[-1]*3))
        z = np.hstack((0., np.cumsum(thk)))
        maxz = max(maxz, z[-1])

        for j in range(nlay):
            rect = Rectangle((x[i] - dxmed2, z[j]), dxmed2 * 2, thk[j])
            patches.append(rect)

    p = PatchCollection(patches, cmap=cmap, linewidths=0)

    if cmin is not None:
        p.set_clim(cmin, cmax)

#    p.set_array( np.log10( vals.ravel() ) )
    setMappableData(p, vals.ravel(), logScale=islog)
    ax.add_collection(p)

    ax.set_ylim((maxz, 0.))
    ax.set_xlim((min(x) - dxmed2, max(x) + dxmed2))
    if title is not None:
        ax.set_title(title)

    pg.mplviewer.createColorbar(p, cMin=cmin, cMax=cmax, nLevs=5)

#    cb = plt.colorbar(p, orientation='horizontal',aspect=50,pad=0.1)
#    xt = [10, 20, 50, 100, 200, 500]
#    cb.set_ticks( xt, [str(xti) for xti in xt] )

    plt.draw()
    return fig, ax


def showStitchedModelsOld(models, x=None, cmin=None, cmax=None,
                          islog=True, title=None):
    """show several 1d block models as (stitched) section"""
    if x is None:
        x = np.arange(len(models))

    nlay = int(np.floor((len(models[0]) - 1) / 2.)) + 1
    if cmin is None or cmax is None:
        cmin = 1e9
        cmax = 1e-9
        for model in models:
            res = np.asarray(model)[nlay - 1:nlay * 2 - 1]
            cmin = min(cmin, min(res))
            cmax = max(cmax, max(res))

        print("cmin=", cmin, " cmax=", cmax)

    dx = np.diff(x)
    dx = np.hstack((dx, dx[-1]))
    x1 = x - dx / 2
    ax = plt.gcf().add_subplot(111)
    ax.cla()
    mapsize = 64
    # cmap = jetmap(mapsize)
    plt.plot(x, np.zeros(len(x)), 'k.')
    maxz = 0.
    for mod in models:
        mod1 = np.asarray(mod)
        res = mod1[nlay - 1:]
        if islog:
            res = np.log(res)
            cmi = np.log(cmin)
            cma = np.log(cmax)
        else:
            cmi = cmin
            cma = cmax

        thk = mod1[:nlay - 1]
        thk = np.hstack((thk, thk[-1]))
        z = np.hstack((0., np.cumsum(thk)))
        maxz = max(maxz, z[-1])
        nres = (res - cmi) / (cma - cmi)
        cind = np.around(nres * mapsize)
        cind[cind >= mapsize] = mapsize - 1
        cind[cind < 0] = 0
        # for j in range(len(thk)):
        #   fc = cmap[cind[j], :]
        #   rect = Rectangle((x1[i], z[j]), dx[i], thk[j], fc=fc)
        #   plt.gca().add_patch(rect)

    ax.set_ylim((maxz, 0.))
    ax.set_xlim((x1[0], x1[-1] + dx[-1]))
    if title is not None:
        plt.title(title)

    plt.draw()
    return


def showfdemsounding(freq, inphase, quadrat, response=None, npl=2):
    """ show FDEM sounding as real(inphase) and imaginary (quadrature)
        fields normalized by the (purely real) free air solution """
    nf = len(freq)
    fig = plt.figure(1)
    fig.clf()
    ax1 = fig.add_subplot(1, npl, npl - 1)
    plt.semilogy(inphase, freq, 'x-')
    if response is not None:
        plt.semilogy(np.asarray(response)[:nf], freq, 'x-')

    plt.grid(which='both')
    ax2 = fig.add_subplot(1, npl, npl)
    plt.semilogy(quadrat, freq, 'x-')
    if response is not None:
        plt.semilogy(np.asarray(response)[nf:], freq, 'x-')
    plt.grid(which='both')
    fig.show()
    ax = [ax1, ax2]
    if npl > 2:
        ax3 = fig.add_subplot(1, npl, 1)
        ax.append(ax3)

    return ax


def insertUnitAtNextLastTick(ax, unit, xlabel=True, position=-2):
    """replaces the last-but-one tick label by unit symbol"""
    if xlabel:
        labels = ax.get_xticks().tolist()
        labels[position] = unit
        ax.set_xticklabels(labels)
    else:
        labels = ax.get_yticks().tolist()
        labels[position] = unit
        ax.set_yticklabels(labels)


def saveResult(fname, data, rrms=None, chi2=None, mode='w'):
    """save rms/chi2 results into filename"""
    with open(fname, mode) as f:
        np.savetxt(f, data)
        if rrms is not None:
            f.write('\nrrms:{}\n'.format(rrms))
        if chi2 is not None:
            f.write('\nchi2:{}\n'.format(chi2))


def getSavePath(folder=None, subfolder='', now=None):
    """TODO"""
    if folder is None:
        path = createResultFolder(subfolder, now)
    else:
        path = createfolders([folder, subfolder])
    return path


def createResultFolder(subfolder, now=None):
    """create a result Folder"""
    result = createDateTimeString(now)
    return createfolders(['./', result, subfolder])


def createDateTimeString(now=None):
    """returns datetime as string (e.g. for saving results)"""
    if now is None:
        now = time.localtime()
    return str(now.tm_year) + str(now.tm_mon).zfill(2) + \
        str(now.tm_mday).zfill(2) + '-' + \
        str(now.tm_hour).zfill(2) + '.' + \
        str(now.tm_min).zfill(2)


def createfolders(foldername_list):
    """
    Creates the folder structure specified by the list.
    """

    path = ''

    for s in foldername_list:
        if s != '/':
            path = path + s + '/'

    try:
        os.makedirs(path)
    except OSError as e:
        if os.path.exists(path):
            print('Path "{}" already exists.'.format(path))
        else:
            print('Unable to create path "{}".'.format(path))
            raise e

    return path
