#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Define special colorbar behavior."""


import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

import pygimli as pg
from .. mplviewer import saveFigure, updateAxes
from .. utils import prettyFloat
from .. core.logger import renameKwarg

def autolevel(z, nLevs, logScale=None, zmin=None, zmax=None):
    """Create nLevs bins for the data array z based on matplotlib ticker.

    Examples
    --------
    >>> import numpy as np
    >>> from pygimli.mplviewer import autolevel
    >>> x = np.linspace(1, 10, 100)
    >>> autolevel(x, 3)
    array([ 1.,  4.,  7., 10.])
    >>> autolevel(x, 3, logScale=True)
    array([  0.1,   1. ,  10. , 100. ])
    """
    locator = None
    if logScale and min(z) > 0:
        locator = ticker.LogLocator()
    else:
        # print('MaxNLocator(nBins=nLevs + 1)', nLevs)
        locator = ticker.LinearLocator(numticks=nLevs+1)
        # locator = ticker.MaxNLocator(nBins=nLevs + 1)
        # locator = ticker.MaxNLocator(nBins='auto')

    if zmin is None:
        zmin = min(z)

    if zmax is None:
        zmax = max(z)

    return locator.tick_values(zmin, zmax)


def cmapFromName(cmapname='jet', ncols=256, bad=None, **kwargs):
    """Get a colormap either from name or from keyworld list.

    See http://matplotlib.org/examples/color/colormaps_reference.html

    Parameters
    ----------
    cmapname : str
        Name for the colormap.

    ncols : int
        Amount of colors.

    bad : [r,g,b,a]
        Default color for bad values [nan, inf] [white]

    ** kwargs :
        cMap : str
            Name for the colormap
        cmap : str
            colormap name (old)
    Returns
    -------
    cMap:
        matplotlib Colormap
    """

    if not bad:
        bad = [1.0, 1.0, 1.0, 0.0]

    renameKwarg('cmap', 'cMap', kwargs)

    if 'cmap' in kwargs:
        cmapname = kwargs.pop('cmap', cmapname)
    elif 'cMap' in kwargs:
        cmapname = kwargs.pop('cMap', cmapname)

    cMap = None
    if cmapname is None:
        cmapname = 'jet'

    if cmapname == 'b2r':
        pg.warn("Don't use manual b2r cMap, use MPL internal 'RdBu' instead.")
        cMap = "RdBu_r"
    else:
        try:
            cMap = mpl.cm.get_cmap(cmapname, ncols)
        except BaseException as e:
            pg.warn("Could not retrieve colormap ", cmapname, e)

    cMap.set_bad(bad)
    return cMap


def findAndMaskBestClim(dataIn, cMin=None, cMax=None, dropColLimitsPerc=5,
                        logScale=True):
    """TODO Documentme."""
    data = np.asarray(dataIn)

    if min(data) < 0:
        logScale = False
    if logScale:
        data = np.log10(data)

    xHist = np.histogram(data, bins=100)[1]

    if not cMin:
        cMin = xHist[dropColLimitsPerc]
        if logScale:
            cMin = pow(10.0, cMin)

    if not cMax:
        cMax = xHist[100 - dropColLimitsPerc]
        if logScale:
            cMax = pow(10.0, cMax)

    if logScale:
        data = pow(10.0, data)

    data[np.where(data < cMin)] = cMin
    data[np.where(data > cMax)] = cMax

    return data, cMin, cMax


def findColorBar(ax):
    """Find the colorbar of an axes.

    Find the colorbar that is associated with given axes or return None.
    """
    pg.critical('Who use this?') #190530
    for i, ai in enumerate(ax.figure.axes):
        print(i, ai)

    for c in ax.collections:
        if isinstance(c, mpl.cm.ScalarMappable):
            if c.colorbar is not None:
                print("cbar:,", c.colorbar)
                return c.colorbar

    raise BaseException("Implement me")

    # print(ax.colorbar)
    # print(ax.images)
    # return None


def updateColorBar(cbar, gci=None, cMin=None, cMax=None, cMap=None,
                   logScale=None, nLevs=5, label=None, **kwargs):
    """Update colorbar values.

    Update limits and label of a given colorbar.

    Parameters
    ----------
    cbar: matplotlib colorbar

    gci : matplotlib graphical instance

    cMin: float

    cMax: float

    cLog: bool

    cMap: matplotlib colormap

    nLevs: int

    label: str

    """
    # print("update cbar:", cMin, cMax, label)
    if gci is not None:
        if min(gci.get_array()) < 1e12:
            norm = mpl.colors.Normalize(vmin=min(gci.get_array()),
                                        vmax=min(gci.get_array()))
            cbar.set_norm(norm)
            gci.set_norm(norm)
        cbar.on_mappable_changed(gci)

    if cMap is not None:
        if isinstance(cMap, str):
            cMap = cmapFromName(cMap, ncols=256, bad=[1.0, 1.0, 1.0, 0.0])

        cbar.mappable.set_cmap(cMap)

    needLevelUpdate = False

    if cMin is not None or cMax is not None or nLevs is not None:
        needLevelUpdate = True

    if logScale is not None:
        needLevelUpdate = True

        if cMin is None:
            cMin = cbar.get_clim()[0]
        if cMax is None:
            cMax = cbar.get_clim()[1]

        if logScale:
            if cMin < 1e-12:
                cMin = min(filter(lambda _x: _x > 0.0,
                                  cbar.mappable.get_array()))

            norm = mpl.colors.LogNorm(vmin=cMin, vmax=cMax)
        else:
            norm = mpl.colors.Normalize(vmin=cMin, vmax=cMax)

        cbar.set_norm(norm)
        cbar.mappable.set_norm(norm)

    if needLevelUpdate:
        setCbarLevels(cbar, cMin, cMax, nLevs)

    if label is not None:
        cbar.set_label(label)

    return cbar


def createColorBar(gci, orientation='horizontal', size=0.2, pad=None,
                   **kwargs):
    """Create a Colorbar.

    Shortcut to create a matplotlib colorbar within the ax for a given
    patchset. The colorbar is stored in the axes object as __cBar__
    to avoid duplicates.

    Parameters
    ----------
    gci: matplotlib graphical instance

    orientation: string

    size: float

    pad: float

    **kwargs :
        Forwarded to updateColorBar
    """
    cbarTarget = plt
    cax = None
    divider = None
    #    if hasattr(patches, 'figure'):
    #       cbarTarget = patches.figure

    ax = None
    if hasattr(gci, 'ax'):
        ax = gci.ax
    if hasattr(gci, 'axes'):
        ax = gci.axes
    elif hasattr(gci, 'get_axes'):
        ax = gci.get_axes()

    cbar = None
    if hasattr(ax, '__cBar__'):
        cbar = ax.__cBar__
        updateColorBar(cbar, gci, **kwargs)
    else:
        divider = make_axes_locatable(ax)

        if divider:
            if orientation == 'horizontal':
                if pad is None:
                    pad = 0.5
                cax = divider.append_axes("bottom", size=size, pad=pad)
            else:
                if pad is None:
                    pad = 0.1
                cax = divider.append_axes("right", size=size, pad=pad)

        cbar = cbarTarget.colorbar(gci, cax=cax, orientation=orientation)
        #store the cbar into the axes to reuse it on the next call
        ax.__cBar__ = cbar

        updateColorBar(cbar, **kwargs)

    return cbar


def createColorBarOnly(cMin=1, cMax=100, logScale=False, cMap=None, nLevs=5,
                       label=None, orientation='horizontal', savefig=None,
                       ax=None, **kwargs):
    """Create figure with a colorbar.

    Create figure with a colorbar.

    Parameters
    ----------
    **kwargs:
        Forwarded to mpl.colorbar.ColorbarBase.

    Returns
    -------
    fig:
        The created figure.

    Examples
    --------
    >>> # import pygimli as pg
    >>> # from pygimli.mplviewer import createColorBarOnly
    >>> # createColorBarOnly(cMin=0.2, cMax=5, logScale=False,
    >>> #                   cMap='b2r',
    >>> #                   nLevs=7,
    >>> #                   label=r'Ratio',
    >>> #                   orientation='horizontal')
    >>> # pg.wait()
    """
    if ax is None:
        fig = plt.figure()
        if orientation is 'horizontal':
            ax = fig.add_axes([0.035, 0.6, 0.93, 0.05])
        else:
            ax = fig.add_axes([0.30, 0.02, 0.22, 0.96])

    norm = None
    if cMin > 0 and logScale is True:
        norm = mpl.colors.LogNorm(vmin=cMin, vmax=cMax)
    else:
        norm = plt.Normalize(vmin=cMin, vmax=cMax)

    cmap = cmapFromName(cMap)

    cbar = mpl.colorbar.ColorbarBase(ax, norm=norm, cmap=cmap,
                                     orientation=orientation, **kwargs)

    #        cbar.labelpad = -20
    #        cbar.ax.yaxis.set_label_position('left')
    updateColorBar(cbar, cMin=cMin, cMax=cMax, nLevs=nLevs, label=label,
                   **kwargs)

    if savefig is not None:
        saveFigure(fig, savefig)

    return ax


def setCbarLevels(cbar, cMin=None, cMax=None, nLevs=5):
    """Set colorbar levels given a number of levels and min/max values."""
    if cMin is None:
        cMin = cbar.get_clim()[0]
    if cMax is None:
        cMax = cbar.get_clim()[1]

    if cMin == cMax:
        cMin *= 0.999
        cMax *= 1.001

    norm = None
    if hasattr(cbar, 'mappable'):
        norm = cbar.mappable.norm
    elif hasattr(cbar, 'norm'):
        norm = cbar.norm

    if isinstance(norm, mpl.colors.LogNorm):
        cbarLevels = np.logspace(np.log10(cMin), np.log10(cMax), nLevs)
    else:
        cbarLevels = np.linspace(cMin, cMax, nLevs)

    # FIXME: [10.1, 10.2, 10.3] mapped to [10 10 10]

    cbarLevelsString = []
    if np.all(np.array(cbarLevels) < 1e-2):
        pg.debug("All values smaller than 1e-4, avoiding additional rounding.")
        roundValue = False
    else:
        roundValue = True

    for i in cbarLevels:
        cbarLevelsString.append(prettyFloat(i, roundValue))
        # print(i, prettyFloat(i))

    if hasattr(cbar, 'mappable'):
        cbar.mappable.set_clim(vmin=cMin, vmax=cMax)
        cbar.set_clim(cMin, cMax)

    cbar.set_ticks(cbarLevels)
    cbar.set_ticklabels(cbarLevelsString)
    cbar.draw_all()

    # necessary since mpl 3.0
    cbar.ax.minorticks_off()


def setMappableValues(mappable, dataIn):
    """Change the data values for a given mapable."""
    pg.critical('remove me')
    data = dataIn
    if not isinstance(data, np.ma.core.MaskedArray):
        data = np.array(dataIn)

    # set bad value color to white
    if mappable.get_cmap() is not None:
        mappable.get_cmap().set_bad([1.0, 1.0, 1.0, 0.0])

    mappable.set_array(data)


def setMappableData(mappable, dataIn, cMin=None, cMax=None, logScale=None,
                    **kwargs):
    """Change the data values for a given mappable.
    """
    data = dataIn
    if not isinstance(data, np.ma.core.MaskedArray):
        data = np.array(dataIn)

    # set bad value color to white
    if mappable.get_cmap() is not None:
        mappable.get_cmap().set_bad([1.0, 1.0, 1.0, 0.0])

    if not cMin:
        cMin = data.min()
    if not cMax:
        cMax = data.max()

    oldLog = None
    if cMin <= 0.0:
        oldLog = isinstance(mappable.norm, mpl.colors.LogNorm)
        if oldLog is True or logScale is True:
            if cMax > 0:
                cMin = min(data[data > 0.0])
                data = np.ma.masked_array(data, data <= 0.0)
            else:
                # if all data are negative switch to lin scale
                return setMappableData(mappable, dataIn, cMin, cMax,
                                       logScale=False)
    if logScale is True:
        mappable.set_norm(mpl.colors.LogNorm())
    elif logScale is False:
        mappable.set_norm(mpl.colors.Normalize(vmin=cMin, vmax=cMax))

    #pg._g(oldLog, logScale, cMin, cMax, mappable.norm)
    mappable.set_array(data)
    mappable.set_clim(cMin, cMax)

    if mappable.colorbar is not None:
        updateColorBar(mappable.colorbar, cMin=cMin, cMax=cMax)

def addCoverageAlpha(patches, coverage, dropThreshold=0.4):
    """Add alpha values to the colors of a polygon collection.

    Parameters
    ----------

    patches : 2D mpl mappable

    coverage : array
        coverage values. Maximum coverage mean no opaqueness.

    dropThreshold : float
        relative minimum coverage
    """
    patches.set_antialiaseds(True)

    # generate individual color values here
    patches.update_scalarmappable()

    cols = patches.get_facecolor()

    C = np.asarray(coverage)
    #    print(np.min(C), np.max(C))

    if (np.min(C) < 0.) | (np.max(C) > 1.) | (np.max(C) < 0.5):

        nn, hh = np.histogram(C, 50)
        nnn = nn.cumsum(axis=0) / float(len(C))

        #        print("min-max nnn ", min(nnn), max(nnn))
        mi = hh[min(np.where(nnn > 0.02)[0])]

        if min(nnn) > dropThreshold:
            ma = max(C)
        else:
            ma = hh[max(np.where(nnn < dropThreshold)[0])]

#            mi = hh[min(np.where(nnn > 0.2)[0])]
#            ma = hh[max(np.where(nnn < 0.7)[0])]

        C = (C - mi) / (ma - mi)
        C[np.where(C < 0.)] = 0.0
        C[np.where(C > 1.)] = 1.0

#    else:
#        print('taking the values directly')

    # add alpha value to the color values
    cols[:, 3] = C

    patches._facecolors = cols

    # delete patch data to avoid automatically rewrite of _facecolors
    patches._A = None

    if hasattr(patches, 'ax'):
        updateAxes(patches.ax)
    elif hasattr(patches, 'get_axes'):
        updateAxes(patches.get_axes())
