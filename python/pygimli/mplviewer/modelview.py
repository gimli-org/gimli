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


def drawModel1D(ax, thickness=None, values=None, model=None, depths=None,
                plot='plot',
                xlabel=r'Resistivity $[\Omega$m$]$', zlabel='Depth [m]',
                z0=0,
                **kwargs):
    """Draw 1d block model into axis ax.

    Draw 1d block model into axis ax defined by values and thickness vectors
    using plot function.
    For log y cases, z0 should be set > 0 so that the default becomes 1.

    Parameters
    ----------
    ax : mpl axes
        Matplotlib Axes object to plot into.

    values : iterable [float]
        [N] Values for each layer plus lower background.

    thickness : iterable [float]
        [N-1] thickness for each layer. Either thickness or depths must be set.

    depths : iterable [float]
        [N-1] Values for layer depths (positive z-coordinates).
        Either thickness or depths must be set.

    model : iterable [float]
        Shortcut to use default model definition.
        thks = model[0:nLay]
        values = model[nLay:]

    plot : string
        mpl plot funktion
        'plot', 'semilogx', 'semilogy', 'loglog'

    xlabel : str
        Label for x axis.

    ylabel : str
        Label for y axis.

    z0 : float
        Starting depth [m]

    **kwargs : dict()
        Forwarded to the plot routine

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import pygimli as pg
    >>> # plt.style.use('ggplot')
    >>> thk = [1, 4, 4]
    >>> res = np.array([10., 5, 15, 50])
    >>> fig, ax = plt.subplots()
    >>> pg.mplviewer.drawModel1D(ax, values=res*5, depths=np.cumsum(thk),
    ...                          plot='semilogx', color='blue')
    >>> pg.mplviewer.drawModel1D(ax, values=res, thickness=thk, z0=1,
    ...                          plot='semilogx', color='red')
    >>> pg.wait()
    """
    if model is not None:
        nLayers = (len(model)-1)//2
        thickness = model[:nLayers]
        values = model[nLayers:]

    if thickness is None and depths is None:
        raise Exception("Either thickness or depths must be given.")

    nLayers = len(values)
    px = np.zeros(nLayers * 2)
    pz = np.zeros(nLayers * 2)

    if thickness is not None:
        z1 = np.cumsum(thickness) + z0
    else:
        z1 = depths

    for i in range(nLayers):
        px[2 * i] = values[i]
        px[2 * i + 1] = values[i]

        if i == nLayers - 1:
            pz[2 * i + 1] = z1[i - 1] * 1.2
        else:
            pz[2 * i + 1] = z1[i]
            pz[2 * i + 2] = z1[i]

    if plot == 'loglog' or plot == 'semilogy':
        if z0 == 0:
            pz[0] = z1[0] / 2.
        else:
            pz[0] = z0

    try:
        plot = getattr(ax, plot)
        plot(px, pz+z0, **kwargs)
    except BaseException as e:
        print(e)

    ax.set_ylabel(zlabel)
    ax.set_xlabel(xlabel)
    # assume positive depths pointing upward
    ax.set_ylim(pz[-1], pz[0])
    ax.grid(True)


def showStitchedModels(models, x=None,
                       logScale=True, title=None, ax=None,
                       zMin=0, zMax=0, zLog=True,
                       cMin=None, cMax=None, cMap='jet',
                       useMesh=False, **kwargs):
    """Show several 1d block models as (stitched) section.

    Parameters
    ----------
    useMesh: bool (False)
        use pygimli mesh instead of patch graph. Experimental.

    """

    nLayers = int(np.floor((len(models[0]) + 1) / 2.))

    fig = None
    if ax is None:
        fig, ax = plt.subplots()

    noXVals = False
    if x is None:
        noXVals = True
        x = np.arange(len(models)) * 1.0

    dxmed2 = np.median(np.diff(x)) / 2.

    vals = np.zeros((len(models), nLayers))

    mesh = None
    if useMesh:
        mesh = pg.Mesh(2)
    patches = []

    zMaxLimit = 0
    z = None

    for i, imod in enumerate(models):
        if isinstance(imod, pg.RVector):
            vals[i, :] = imod(nLayers - 1, 2 * nLayers - 1)
            thk = np.asarray(imod(0, nLayers - 1))
        else:
            vals[i, :] = imod[nLayers - 1:2 * nLayers - 1]
            thk = imod[:nLayers - 1]

        if zMax > 0:
            z = np.hstack((0., np.cumsum(thk)))
            z = np.hstack((z, zMax))
        else:
            thk = np.hstack((thk, thk[-1]*3))
            z = np.hstack((0., np.cumsum(thk)))

        zMaxLimit = max(zMaxLimit, z[-1])

        if mesh is not None:
            for j in range(nLayers):
                #n1 = mesh.createNode([x[i],   z[j]])
                #n2 = mesh.createNode([x[i+1], z[j]])
                #n3 = mesh.createNode([x[i+1], z[j+1]])
                #n4 = mesh.createNode([x[i],   z[j+1]])

                n1 = mesh.createNode([x[i] - dxmed2, -z[j]])
                n2 = mesh.createNode([x[i] + dxmed2, -z[j]])
                n3 = mesh.createNode([x[i] + dxmed2, -z[j+1]])
                n4 = mesh.createNode([x[i] - dxmed2, -z[j+1]])

                mesh.createQuadrangle(n1, n2, n3, n4)
        else:
            for j in range(nLayers):
                rect = Rectangle((x[i] - dxmed2, z[j]),
                                dxmed2 * 2, z[j+1]-z[j])
                patches.append(rect)


    if mesh is not None:
        pg.show(mesh, vals.ravel(), ax=ax, cMin=cMin, cMax=cMax, cMap=cMap,
                label='Parameters')
    else:

        if 'cmap' in kwargs:
            print("DeprecationWarning: please use 'cMap' instead of 'cmap'")
            cMap = kwargs.pop('cmap', cMap)

        p = PatchCollection(patches, cmap=cMap, linewidths=0)

        if cMin is not None and cMax is not None:
            p.set_clim(cMin, cMax)

        if 'islog' in kwargs:
            print("DeprecationWarning: pleas use 'logScale' instead of 'islog'")
            logScale = kwargs.pop('islog', logScale)

        setMappableData(p, vals.ravel(), logScale=logScale)
        ax.add_collection(p)

        ax.set_ylim((zMaxLimit, zMin))

        if zLog:
            ax.set_yscale("log", nonposy='clip')

        ax.set_xlim((min(x) - dxmed2, max(x) + dxmed2))

        if title is not None:
            ax.set_title(title)

        if cMin is not None and cMax is not None:
            pg.mplviewer.createColorBar(p, cMin=cMin, cMax=cMax, nLevs=5)

    #    cb = plt.colorbar(p, orientation='horizontal',aspect=50,pad=0.1)
    #    xt = [10, 20, 50, 100, 200, 500]
    #    cb.set_ticks( xt, [str(xti) for xti in xt] )

    if noXVals:
        #ax.figure.set_size((20,20))
        #ax.set_aspect(0.2)
        ax.set_aspect(abs(max(x))/max(abs(z))/3, adjustable='box')
        #plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        #plt.tight_layout()
        #ax.figure.set_size_inches(200,20)
        #ax.axis('scaled')
        #plt.axis('tight')

    return ax

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


def draw1dmodel(x, thk=None, xlab=None, zlab="z in m", islog=True, z0=0):
    """DEPRECATED."""
    print("STYLE_WARNING!!!!!!! don't use this call. "
          "Use show1dmodel or drawModel1D instead.")
    show1dmodel(x, thk, xlab, zlab, islog, z0)


def show1dmodel(x, thk=None, xlab=None, zlab="z in m", islog=True, z0=0,
                **kwargs):
    """Show 1d block model defined by value and thickness vectors."""
    print("STYLE_WARNING!!!!!!! don't use this call. "
          "WHO use this anymore??.")

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
        plt.semilogx(px, pz, **kwargs)
    else:
        plt.plot(px, pz, **kwargs)

    plt.ion()
    plt.grid(which='both')
    plt.xlim((np.min(x) * 0.9, np.max(x) * 1.1))
    plt.ylim((max(z1) * 1.15, 0.))
    plt.xlabel(xlab)
    plt.ylabel(zlab)
    plt.show()
    return


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
