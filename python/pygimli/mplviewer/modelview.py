# -*- coding: utf-8 -*-
"""Model viewer functions."""

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm

import pygimli as pg

from pygimli.mplviewer.colorbar import setMappableData
from pygimli.utils import rndig


def showmymatrix(mat, x, y, dx=2, dy=1, xlab=None, ylab=None, cbar=None):
    """What is this good for?."""
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


def draw1dmodelErr(x, xL, xU=None, thk=None, xcol='g', ycol='r', **kwargs):
    """TODO."""
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
    """Draw 1d model with lower and upper bounds."""
    raise BaseException("IMPLEMENTME")
    # draw1dmodel(x, thk, color='red', **kwargs)
    # for i in range(len(x)):
    #     x1 = np.array(x)
    #     x1[i] = xL[i]
    #     draw1dmodel(x1, thk, color='blue')
    #     x1[i] = xU[i]
    #     draw1dmodel(x1, thk, color='blue')
    #
    # li = draw1dmodel(x, thk, color='red', **kwargs)
    # plt.xlim((min(xL) * 0.9, max(xU) * 1.1))
    # return li


def showStitchedModels(models, ax=None, x=None, cmin=None, cmax=None,
                       islog=True, title=None, cmap='jet'):
    """Show several 1d block models as (stitched) section."""
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


def showStitchedModels_Redundant(mods, ax=None,
                                 cmin=None, cmax=None, **kwargs):
    """Show several 1d block models as (stitched) section."""
    x = kwargs.pop('x', np.arange(len(mods)))
    topo = kwargs.pop('topo', x*0)

    nlay = int(np.floor((len(mods[0]) - 1) / 2.)) + 1
    if cmin is None or cmax is None:
        cmin = 1e9
        cmax = 1e-9
        for model in mods:
            res = np.asarray(model)[nlay - 1:nlay * 2 - 1]
            cmin = min(cmin, min(res))
            cmax = max(cmax, max(res))

    if kwargs.pop('sameSize', True):  # all having the same width
        dx = np.ones_like(x)*np.median(np.diff(x))
    else:
        dx = np.diff(x) * 1.05
        dx = np.hstack((dx, dx[-1]))

    x1 = x - dx / 2
    if ax is None:
        fig, ax = plt.subplots()
    else:
        ax = ax
        fig = ax.figure

#    ax.plot(x, x * 0., 'k.')
    zm = kwargs.pop('zm', None)
    maxz = 0.
    if zm is not None:
        maxz = zm
    recs = []
    RES = []
    for i, mod in enumerate(mods):
        mod1 = np.asarray(mod)
        res = mod1[nlay - 1:]
        RES.extend(res)

        thk = mod1[:nlay - 1]
        thk = np.hstack((thk, thk[-1]))
        z = np.hstack((0., np.cumsum(thk)))
        if zm is not None:
            thk[-1] = zm - z[-2]
            z[-1] = zm
        else:
            maxz = max(maxz, z[-1])

        for j, _ in enumerate(thk):
            recs.append(Rectangle((x1[i], topo[i]-z[j]), dx[i], -thk[j]))

    pp = PatchCollection(recs, edgecolors=kwargs.pop('edgecolors', 'none'))
    pp.set_edgecolor(kwargs.pop('edgecolors', 'none'))
    pp.set_linewidths(0.0)
    ax.add_collection(pp)
    if 'cmap' in kwargs:
        pp.set_cmap(kwargs['cmap'])

    print(cmin, cmax)
    norm = LogNorm(cmin, cmax)
    pp.set_norm(norm)
    pp.set_array(np.array(RES))
#    pp.set_clim(cmin, cmax)
    ax.set_ylim((-maxz, max(topo)))
    ax.set_xlim((x1[0], x1[-1] + dx[-1]))

    cbar = None
    if kwargs.pop('colorBar', True):
        cbar = plt.colorbar(pp, ax=ax, norm=norm, orientation='horizontal',
                            aspect=60)  # , ticks=[1, 3, 10, 30, 100, 300])
        if 'ticks' in kwargs:
            cbar.set_ticks(kwargs['ticks'])
#        cbar.autoscale_None()
    if ax is None:  # newly created fig+ax
        return fig, ax
    else:  # already given, better give back color bar
        return cbar


def showStitchedModelsOld(models, x=None, cmin=None, cmax=None,
                          islog=True, title=None):
    """Show several 1d block models as (stitched) section."""
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


def insertUnitAtNextLastTick(ax, unit, xlabel=True, position=-2):
    """Replace the last-but-one tick label by unit symbol."""
    if xlabel:
        labels = ax.get_xticks().tolist()
        labels[position] = unit
        ax.set_xticklabels(labels)
    else:
        labels = ax.get_yticks().tolist()
        labels[position] = unit
        ax.set_yticklabels(labels)


def drawModel1D(ax, thickness, values, plotfunction='plot',
                xlabel=r'Resistivity $[\Omega$ m$]$', **kwargs):
    """Draw 1d block model into axis ax.

    Draw 1d block model into axis ax defined by values and thickness vectors
    using plotfunction.
    """
    nLayers = len(thickness) + 1
    px = np.zeros(nLayers * 2)
    pz = np.zeros(nLayers * 2)
    z1 = np.cumsum(thickness)

    for i in range(nLayers):
        px[2 * i] = values[i]
        px[2 * i + 1] = values[i]

        if i == nLayers - 1:
            pz[2 * i + 1] = z1[i - 1] * 1.2
        else:
            pz[2 * i + 1] = z1[i]
            pz[2 * i + 2] = z1[i]

    if plotfunction == 'loglog' or plotfunction == 'semilogy':
        pz[0] = thickness[0] * 0.8

    try:
        plot = getattr(ax, plotfunction)
        plot(px, pz, **kwargs)
    except BaseException as e:
        print(e)

    ax.set_ylabel('Depth [m]')
    ax.set_xlabel(xlabel)
    ax.set_ylim(pz[-1], pz[0])
    ax.grid()
    return ax

# def draw1dmodel(... )


def draw1dmodel(x, thk=None, xlab=None, zlab="z in m", islog=True, z0=0):
    """DEPRECATED."""
    print("STYLE_WARNING!!!!!!! don't use this call. Use show1dmodel instead.")
    show1dmodel(x, thk, xlab, zlab, islog, z0)


def show1dmodel(x, thk=None, xlab=None, zlab="z in m", islog=True, z0=0):
    """Show 1d block model defined by value and thickness vectors."""
    if xlab is None:
        xlab = "$\\rho$ in $\\Omega$m"

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
    if islog:
        plt.semilogx(px, pz)
    else:
        plt.plot(px, pz)

    plt.ion()
    plt.grid(which='both')
    plt.xlim((np.min(x) * 0.9, np.max(x) * 1.1))
    plt.ylim((max(z1) * 1.15, 0.))
    plt.xlabel(xlab)
    plt.ylabel(zlab)
    plt.show()
    return


def draw1dmodel__Redundant(x, thk=None, xlab=None, zlab="z in m", islog=True,
                           fs=14, z0=0, **kwargs):
    """Draw 1d block model defined by value and thickness vectors."""
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


def showfdemsounding(freq, inphase, quadrat, response=None, npl=2):
    """Show FDEM sounding as real(inphase) and imaginary (quadrature) fields.

    Show FDEM sounding as real(inphase) and imaginary (quadrature) fields
        normalized by the (purely real) free air solution.
    """
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
