# -*- coding: utf-8 -*-
"""
    Define special colorbar behavior.
"""
import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable

cdict = {'red': ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)),
         'green': ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)),
         'blue': ((0.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0))}

blueRedCMap = mpl.colors.LinearSegmentedColormap('my_colormap', cdict, 256)


def autolevel(z, N, logscale=None):
    """Create N levels for the data array z based on matplotlib ticker.

    Examples
    --------
    >>> import numpy as np
    >>> from pygimli.mplviewer import autolevel
    >>> x = np.linspace(1,10,100)
    >>> autolevel(x, 3)
    array([  0. ,   2.5,   5. ,   7.5,  10. ])
    >>> autolevel(x, 3, logscale=True)
    array([   0.1,    1. ,   10. ,  100. ])
    """
    locator = None
    if logscale:
        locator = ticker.LogLocator()
    else:
        locator = ticker.MaxNLocator(N + 1)

    zmin = min(z)
    zmax = max(z)

    return locator.tick_values(zmin, zmax)


def cmapFromName(cmapname, ncols=256, bad=None):
    """
        Do we need this?
    """
    if not bad:
        bad = [1.0, 1.0, 1.0, 0.0]

    cmap = mpl.cm.get_cmap('jet', ncols)

    if cmapname is not None:
        if cmapname == 'b2r':
            cmap = mpl.colors.LinearSegmentedColormap('my_colormap',
                                                      cdict, ncols)
        else:
            try:
                cmap = mpl.cm.get_cmap(cmapname, ncols)
            except Exception as e:
                print("Could not retrieve colormap ", cmapname, e)

    cmap.set_bad(bad)
    return cmap


def findAndMaskBestClim(dataIn, cMin=None, cMax=None,
                        dropColLimitsPerc=5, logScale=True):
    """What is this?"""
    data = np.asarray(dataIn)

    # if type( dataIn ) == g.RVector:
    # data = np.asarray( dataIn )
    # elif type( dataIn ) == list:
    # data = np.array( dataIn )
    # else:
    # data = array( dataIn )

    if (min(data) < 0):
        logScale = False
    if (logScale):
        data = np.log10(data)

    Nhist, xHist = np.histogram(data, bins=100)

    if not cMin:
        cMin = xHist[dropColLimitsPerc]
        if logScale:
            cMin = pow(10.0, cMin)

    if not cMax:
        cMax = xHist[100 - dropColLimitsPerc]
        if logScale:
            cMax = pow(10.0, cMax)

    if (logScale):
        data = pow(10.0, data)

    data[np.where(data < cMin)] = cMin
    data[np.where(data > cMax)] = cMax

    return data, cMin, cMax


def createColorbar(patches, cMin=None, cMax=None, nLevs=5,
                   label=None, orientation='horizontal', **kwargs):
    """
    Create a Colobar.

    Shortcut to create a matplotlib colorbar within the axes for a given
    patchset.

    Parameters
    ----------
    **kwargs :
        * size : with or height of the colobar
        * pad : padding distance from axes
    """
    cbarTarget = plt
    cax = None
    divider = None
#    if hasattr(patches, 'figure'):
#       cbarTarget = patches.figure

    if hasattr(patches, 'ax'):
        divider = make_axes_locatable(patches.ax)
    elif hasattr(patches, 'get_axes'):
        divider = make_axes_locatable(patches.get_axes())

    if divider:
        if orientation == 'horizontal':
            size = kwargs.pop('size', 0.2)
            pad = kwargs.pop('pad', 0.5)
            cax = divider.append_axes("bottom", size=size, pad=pad)
        else:
            size = kwargs.pop('size', 0.2)
            pad = kwargs.pop('pad', 0.1)
            cax = divider.append_axes("right", size=size, pad=pad)

    cbar = cbarTarget.colorbar(patches, cax=cax,
                               orientation=orientation)

    setCbarLevels(cbar, cMin, cMax, nLevs)

    if label is not None:
        cbar.set_label(label)

    return cbar


def setCbarLevels(cbar, cMin=None, cMax=None, nLevs=5):
    """What's that"""

    # print "setCbarLevels", cMin, cMax

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
        cMin = norm.vmin
        cMax = norm.vmax

    if isinstance(norm, mpl.colors.LogNorm):
        cbarLevels = np.logspace(np.log10(cMin), np.log10(cMax), nLevs)
    else:
        cbarLevels = np.linspace(cMin, cMax, nLevs)

    # print cbarLevels
    cbarLevelsString = []
    for i in cbarLevels:
        if abs(i) == 0.0:
            cbarLevelsString.append("0")
        elif abs(i) > 1e4 or abs(i) <= 1e-4:
            cbarLevelsString.append("%.1e" % i)
        elif abs(i) < 1e-3:
            cbarLevelsString.append("%.5f" % i)
        elif abs(i) < 1e-2:
            cbarLevelsString.append("%.4f" % i)
        elif abs(i) < 1e-1:
            cbarLevelsString.append("%.3f" % i)
        elif abs(i) < 1e0:
            cbarLevelsString.append("%.2f" % i)
        elif abs(i) < 1e1:
            cbarLevelsString.append("%.1f" % i)
        else:
            cbarLevelsString.append("%.0f" % i)

    if hasattr(cbar, 'mappable'):
        cbar.mappable.set_clim(cMin, cMax)

    # print ticks, cbarLevels, cbarLevelsString
    # if cbar.orientation == 'vertical':
    cbar.set_ticks(cbarLevels)
    cbar.set_ticklabels(cbarLevelsString)

    # print cbar._ticker()

    #    else:
    #        cbar.ax.set_xticks( ticks )
    #        cbar.ax.set_xticklabels( cbarLevelsString )

    cbar.draw_all()


def setMappableData(mappable, dataIn, cMin=None, cMax=None, logScale=False):
    """
        Change the data values for a given mappable.
    """

    data = dataIn

    if not isinstance(data, np.ma.core.MaskedArray):
        data = np.array(dataIn)

    if logScale and data.min() <= 0:
        data = np.ma.masked_array(data, data <= 0.0)

    # set bad value color to white
    if mappable.get_cmap() is not None:
        mappable.get_cmap().set_bad([1.0, 1.0, 1.0, 0.0])

    if not cMin:
        cMin = data.min()
    if not cMax:
        cMax = data.max()

    if cMin > 0.0 and logScale:
        mappable.set_norm(mpl.colors.LogNorm())
    else:
        mappable.set_norm(mpl.colors.Normalize())

#    print("set mappable data, log: ", logScale, "cmin: ", cMin, "cmax: ", cMax)
    mappable.set_array(data)
    # mappable.set_level(10)
    mappable.set_clim(cMin, cMax)


def addCoverageAlpha(patches, coverage, dropThreshold=0.4):
    """
    Add alpha values to the colors of a polygon collection.

    Parameters
    ----------

    patches : 2D mpl mappable

    coverage : array
        coverage values. Maximum coverage mean no opaqueness.

    dropThreshold : float
        relative minimum coverage
    """

    patches.set_antialiaseds(True)
    patches.set_linewidth(0.000)

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
    # patches._edgecolor = 'None'

    # delete patch data to avoid automatically rewrite of _facecolors
    patches._A = None
