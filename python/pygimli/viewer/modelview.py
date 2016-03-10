# -*- coding: utf-8 -*-
"""
    pygimli model viewer functions.
"""

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm


def drawModel1D(ax, thickness, values, plotfunction='plot',
                xlabel='Resistivity $[\Omega$ m$]$', *args, **kwargs):
    """Draw 1d block model into axis ax defined by values and thickness vectors
    using plotfunction."""

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
        plot(px, pz, *args, **kwargs)
    except Exception as e:
        print(e)

    ax.set_ylabel('Depth [m]')
    ax.set_xlabel(xlabel)
    ax.set_ylim(pz[-1], pz[0])
    ax.grid()
    return ax

# def draw1dmodel(... )


def draw1dmodel(x, thk=None, xlab=None, zlab="z in m", islog=True, z0=0):
    """
        DEPRECATED
    """
    print("STYLE_WARNING!!!!!!! don't use this call. Use show1dmodel instead.")
    show1dmodel(x, thk, xlab, zlab, islog, z0)


def show1dmodel(x, thk=None, xlab=None, zlab="z in m", islog=True, z0=0):
    """draw 1d block model defined by value and thickness vectors."""
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


def showStitchedModels(mods, axes=None, cmin=None, cmax=None, **kwargs):
    """
        Show several 1d block models as (stitched) section.
    """
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
    if axes is None:
        fig, ax = plt.subplots()
    else:
        ax = axes
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

        for j in range(len(thk)):
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
    if axes is None:  # newly created fig+ax
        return fig, ax
    else:  # already given, better give back color bar
        return cbar


def showfdemsounding(freq, inphase, quadrat, response=None, npl=2):
    """
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
