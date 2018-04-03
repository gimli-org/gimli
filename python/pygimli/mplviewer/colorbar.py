#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Define special colorbar behavior."""

from distutils.version import LooseVersion

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

from pygimli.mplviewer import saveFigure, updateAxes


cdict = {'red': ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)),
         'green': ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)),
         'blue': ((0.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0))}

blueRedCMap = mpl.colors.LinearSegmentedColormap('my_colormap', cdict, 256)

def autolevel(z, nLevs, logScale=None, zmin=None, zmax=None):
    """Create nLevs bins for the data array z based on matplotlib ticker.

    Examples
    --------
    >>> import numpy as np
    >>> from pygimli.mplviewer import autolevel
    >>> x = np.linspace(1, 10, 100)
    >>> autolevel(x, 3)
    array([  1.,   4.,   7.,  10.])
    >>> autolevel(x, 3, logScale=True)
    array([   0.1,    1. ,   10. ,  100. ])
    """
    locator = None
    if logScale and min(z) > 0:
        locator = ticker.LogLocator()
    else:
        #print('MaxNLocator(nBins=nLevs + 1)', nLevs)
        locator = ticker.LinearLocator(numticks=nLevs+1)
        #locator = ticker.MaxNLocator(nBins=nLevs + 1)
        #locator = ticker.MaxNLocator(nBins='auto')

    if zmin is None:
        zmin = round(min(z), 2)

    if zmax is None:
        zmax = round(max(z), 2)

    #print("autolevel", zmin, zmax)
    #print(locator.tick_values(zmin, zmax))
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
    cmap:
        matplotlib Colormap
    """

    if not bad:
        bad = [1.0, 1.0, 1.0, 0.0]

    if 'cmap' in kwargs:
        cmapname = kwargs.pop('cmap', cmapname)
    elif 'cMap' in kwargs:
        cmapname = kwargs.pop('cMap', cmapname)

    cmap = None
    if cmapname is None:
        cmapname = 'jet'

    if cmapname == 'b2r':
        cmap = mpl.colors.LinearSegmentedColormap('my_colormap', cdict, ncols)
    elif cmapname == 'viridis' and \
            LooseVersion(mpl.__version__) < LooseVersion('1.5.0'):

        print("Mpl:", mpl.__version__, " using HB viridis")
        cmap = LinearSegmentedColormap.from_list('viridis', viridis_data[::-1])
    elif cmapname == 'viridis_r':
        print("Using HB viridis_r")
        cmap = LinearSegmentedColormap.from_list('viridis', viridis_data)
    else:
        try:
            cmap = mpl.cm.get_cmap(cmapname, ncols)
        except BaseException as e:
            print("Could not retrieve colormap ", cmapname, e)

    cmap.set_bad(bad)
    return cmap


def findAndMaskBestClim(dataIn,
                        cMin=None,
                        cMax=None,
                        dropColLimitsPerc=5,
                        logScale=True):
    """TODO Documentme."""
    data = np.asarray(dataIn)

    # if type( dataIn ) == g.RVector:
    # data = np.asarray( dataIn )
    # elif type( dataIn ) == list:
    # data = np.array( dataIn )
    # else:
    # data = array( dataIn )

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
                   logScale=None, nLevs=5, label=None):
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
        pass
        # check the following first
        # cbar.on_mappable_changed(gci)

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
                cMin = min(filter(lambda _x: _x > 0.0, cbar.mappable.get_array()))

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


def createColorBar(gci, orientation='horizontal', size=0.2, pad=None, **kwargs):
    """Create a Colorbar.

    Shortcut to create a matplotlib colorbar within the ax for a given
    patchset.

    Parameters
    ----------

    gci : matplotlib graphical instance

    orientation : string

    size : float

    pad : float


    **kwargs :
        Forwarded to updateColorBar
    """
    cbarTarget = plt
    cax = None
    divider = None
    #    if hasattr(patches, 'figure'):
    #       cbarTarget = patches.figure

    if hasattr(gci, 'ax'):
        divider = make_axes_locatable(gci.ax)
    if hasattr(gci, 'axes'):
        divider = make_axes_locatable(gci.axes)
    elif hasattr(gci, 'get_axes'):
        divider = make_axes_locatable(gci.get_axes())

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

    updateColorBar(cbar, **kwargs)

    return cbar


def createColorBarOnly(cMin=1, cMax=100, logScale=False, cMap=None, nLevs=5,
                       label=None, orientation='horizontal', savefig=None,
                       **kwargs):
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
    updateColorBar(cbar, cMin=cMin, cMax=cMax, nLevs=nLevs, label=label
                   **kwargs)

    if savefig is not None:
        saveFigure(fig, savefig)

    return fig


def prettyFloat(v):
    """Return a pretty string for a given value suitable for graphical output."""
    if abs(round(v)-v) < 1e-4 and abs(v) < 1e3:
        return str(int(v))
    elif abs(v) == 0.0:
        return "0"
    elif abs(v) > 1e3 or abs(v) <= 1e-3:
        return str("%.1e" % v)
    elif abs(v) < 1e-2:
        return str("%.4f" % v)
    elif abs(v) < 1e-1:
        return str("%.3f" % v)
    elif abs(v) < 1e0:
        return str("%.2f" % v)
    elif abs(v) < 1e1:
        return str("%.1f" % v)
    elif abs(v) < 1e2:
        return str("%.1f" % v)
    else:
        return str("%.0f" % v)


def setCbarLevels(cbar, cMin=None, cMax=None, nLevs=5):
    """TODO Documentme."""
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

    #print(cbarLevels)
    # FIXME: [10.1, 10.2, 10.3] mapped to [10 10 10]

    cbarLevelsString = []
    for i in cbarLevels:
        cbarLevelsString.append(prettyFloat(i))

    if hasattr(cbar, 'mappable'):
        cbar.mappable.set_clim(vmin=cMin, vmax=cMax)
        cbar.set_clim(cMin, cMax)

    #print('setCbarLevels.ticks ********************', cbarLevels)
    #print('setCbarLevels.ticklabels ********************', cbarLevelsString)
    cbar.set_ticks(cbarLevels)
    cbar.set_ticklabels(cbarLevelsString)

    cbar.draw_all()


def setMappableValues(mappable, dataIn):
    """Change the data values for a given mapable."""
    data = dataIn
    if not isinstance(data, np.ma.core.MaskedArray):
        data = np.array(dataIn)

    # set bad value color to white
    if mappable.get_cmap() is not None:
        mappable.get_cmap().set_bad([1.0, 1.0, 1.0, 0.0])

    mappable.set_array(data)


def setMappableData(mappable, dataIn, cMin=None, cMax=None, logScale=False):
    """Change the data values for a given mappable.

    DEPRECATED
    """
    #DEPRECATED
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

#   print("set mappable data, log: ", logScale, "cmin: ", cMin, "cmax: ", cMax)
    mappable.set_array(data)
    # mappable.set_level(10)
    mappable.set_clim(cMin, cMax)


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

    if hasattr(patches, 'ax'):
        updateAxes(patches.ax)
    elif hasattr(patches, 'get_axes'):
        updateAxes(patches.get_axes())

parameters = {'xp':
              [22.674387857633945, 11.221508276482126, -14.356589454756971,
               -47.18817758739222, -34.59001004812521, -6.0516291196352654],
              'yp': [-20.102530541012214, -33.08246073298429,
                     -42.24476439790574, -5.595549738219887, 42.5065445026178,
                     40.13395157135497],
              'min_JK': 18.8671875,
              'max_JK': 92.5}

viridis_data = [
    [0.26700401, 0.00487433, 0.32941519], [0.26851048, 0.00960483, 0.33542652],
    [0.26994384, 0.01462494, 0.34137895], [0.27130489, 0.01994186, 0.34726862],
    [0.27259384, 0.02556309, 0.35309303], [0.27380934, 0.03149748, 0.35885256],
    [0.27495242, 0.03775181, 0.36454323], [0.27602238, 0.04416723, 0.37016418],
    [0.27701840, 0.05034437, 0.37571452], [0.27794143, 0.05632444, 0.38119074],
    [0.27879067, 0.06214536, 0.38659204], [0.27956550, 0.06783587, 0.39191723],
    [0.28026658, 0.07341724, 0.39716349], [0.28089358, 0.07890703, 0.40232944],
    [0.28144581, 0.08431970, 0.40741404], [0.28192358, 0.08966622, 0.41241521],
    [0.28232739, 0.09495545, 0.41733086], [0.28265633, 0.10019576, 0.42216032],
    [0.28291049, 0.10539345, 0.42690202], [0.28309095, 0.11055307, 0.43155375],
    [0.28319704, 0.11567966, 0.43611482], [0.28322882, 0.12077701, 0.44058404],
    [0.28318684, 0.12584799, 0.44496000], [0.28307200, 0.13089477, 0.44924127],
    [0.28288389, 0.13592005, 0.45342734], [0.28262297, 0.14092556, 0.45751726],
    [0.28229037, 0.14591233, 0.46150995], [0.28188676, 0.15088147, 0.46540474],
    [0.28141228, 0.15583425, 0.46920128], [0.28086773, 0.16077132, 0.47289909],
    [0.28025468, 0.16569272, 0.47649762], [0.27957399, 0.17059884, 0.47999675],
    [0.27882618, 0.17549020, 0.48339654], [0.27801236, 0.18036684, 0.48669702],
    [0.27713437, 0.18522836, 0.48989831], [0.27619376, 0.19007447, 0.49300074],
    [0.27519116, 0.19490540, 0.49600488], [0.27412802, 0.19972086, 0.49891131],
    [0.27300596, 0.20452049, 0.50172076], [0.27182812, 0.20930306, 0.50443413],
    [0.27059473, 0.21406899, 0.50705243], [0.26930756, 0.21881782, 0.50957678],
    [0.26796846, 0.22354911, 0.51200840], [0.26657984, 0.22826210, 0.51434870],
    [0.26514450, 0.23295593, 0.51659930], [0.26366320, 0.23763078, 0.51876163],
    [0.26213801, 0.24228619, 0.52083736], [0.26057103, 0.24692170, 0.52282822],
    [0.25896451, 0.25153685, 0.52473609], [0.25732244, 0.25613040, 0.52656332],
    [0.25564519, 0.26070284, 0.52831152], [0.25393498, 0.26525384, 0.52998273],
    [0.25219404, 0.26978306, 0.53157905], [0.25042462, 0.27429024, 0.53310261],
    [0.24862899, 0.27877509, 0.53455561], [0.24681140, 0.28323662, 0.53594093],
    [0.24497208, 0.28767547, 0.53726018], [0.24311324, 0.29209154, 0.53851561],
    [0.24123708, 0.29648471, 0.53970946], [0.23934575, 0.30085494, 0.54084398],
    [0.23744138, 0.30520222, 0.54192140], [0.23552606, 0.30952657, 0.54294396],
    [0.23360277, 0.31382773, 0.54391424], [0.23167350, 0.31810580, 0.54483444],
    [0.22973926, 0.32236127, 0.54570633], [0.22780192, 0.32659432, 0.54653200],
    [0.22586330, 0.33080515, 0.54731353], [0.22392515, 0.33499400, 0.54805291],
    [0.22198915, 0.33916114, 0.54875211], [0.22005691, 0.34330688, 0.54941304],
    [0.21812995, 0.34743154, 0.55003755], [0.21620971, 0.35153548, 0.55062743],
    [0.21429757, 0.35561907, 0.55118440], [0.21239477, 0.35968273, 0.55171011],
    [0.21050310, 0.36372671, 0.55220646], [0.20862342, 0.36775151, 0.55267486],
    [0.20675628, 0.37175775, 0.55311653], [0.20490257, 0.37574589, 0.55353282],
    [0.20306309, 0.37971644, 0.55392505], [0.20123854, 0.38366989, 0.55429441],
    [0.19942950, 0.38760678, 0.55464205], [0.19763650, 0.39152762, 0.55496905],
    [0.19585993, 0.39543297, 0.55527637], [0.19410009, 0.39932336, 0.55556494],
    [0.19235719, 0.40319934, 0.55583559], [0.19063135, 0.40706148, 0.55608907],
    [0.18892259, 0.41091033, 0.55632606], [0.18723083, 0.41474645, 0.55654717],
    [0.18555593, 0.41857040, 0.55675292], [0.18389763, 0.42238275, 0.55694377],
    [0.18225561, 0.42618405, 0.55712010], [0.18062949, 0.42997486, 0.55728221],
    [0.17901879, 0.43375572, 0.55743035], [0.17742298, 0.43752720, 0.55756466],
    [0.17584148, 0.44128981, 0.55768526], [0.17427363, 0.44504410, 0.55779216],
    [0.17271876, 0.44879060, 0.55788532], [0.17117615, 0.45252980, 0.55796464],
    [0.16964573, 0.45626209, 0.55803034], [0.16812641, 0.45998802, 0.55808199],
    [0.16661710, 0.46370813, 0.55811913], [0.16511703, 0.46742290, 0.55814141],
    [0.16362543, 0.47113278, 0.55814842], [0.16214155, 0.47483821, 0.55813967],
    [0.16066467, 0.47853961, 0.55811466], [0.15919413, 0.48223740, 0.55807280],
    [0.15772933, 0.48593197, 0.55801347], [0.15626973, 0.48962370, 0.55793600],
    [0.15481488, 0.49331293, 0.55783967], [0.15336445, 0.49700003, 0.55772371],
    [0.15191820, 0.50068529, 0.55758733], [0.15047605, 0.50436904, 0.55742968],
    [0.14903918, 0.50805136, 0.55725050], [0.14760731, 0.51173263, 0.55704861],
    [0.14618026, 0.51541316, 0.55682271], [0.14475863, 0.51909319, 0.55657181],
    [0.14334327, 0.52277292, 0.55629491], [0.14193527, 0.52645254, 0.55599097],
    [0.14053599, 0.53013219, 0.55565893], [0.13914708, 0.53381201, 0.55529773],
    [0.13777048, 0.53749213, 0.55490625], [0.13640850, 0.54117264, 0.55448339],
    [0.13506561, 0.54485335, 0.55402906], [0.13374299, 0.54853458, 0.55354108],
    [0.13244401, 0.55221637, 0.55301828], [0.13117249, 0.55589872, 0.55245948],
    [0.12993270, 0.55958162, 0.55186354], [0.12872938, 0.56326503, 0.55122927],
    [0.12756771, 0.56694891, 0.55055551], [0.12645338, 0.57063316, 0.54984110],
    [0.12539383, 0.57431754, 0.54908564], [0.12439474, 0.57800205, 0.54828740],
    [0.12346281, 0.58168661, 0.54744498], [0.12260562, 0.58537105, 0.54655722],
    [0.12183122, 0.58905521, 0.54562298], [0.12114807, 0.59273889, 0.54464114],
    [0.12056501, 0.59642187, 0.54361058], [0.12009154, 0.60010387, 0.54253043],
    [0.11973756, 0.60378459, 0.54139999], [0.11951163, 0.60746388, 0.54021751],
    [0.11942341, 0.61114146, 0.53898192], [0.11948255, 0.61481702, 0.53769219],
    [0.11969858, 0.61849025, 0.53634733], [0.12008079, 0.62216081, 0.53494633],
    [0.12063824, 0.62582833, 0.53348834], [0.12137972, 0.62949242, 0.53197275],
    [0.12231244, 0.63315277, 0.53039808], [0.12344358, 0.63680899, 0.52876343],
    [0.12477953, 0.64046069, 0.52706792], [0.12632581, 0.64410744, 0.52531069],
    [0.12808703, 0.64774881, 0.52349092], [0.13006688, 0.65138436, 0.52160791],
    [0.13226797, 0.65501363, 0.51966086], [0.13469183, 0.65863619, 0.51764880],
    [0.13733921, 0.66225157, 0.51557101], [0.14020991, 0.66585927, 0.51342680],
    [0.14330291, 0.66945881, 0.51121549], [0.14661640, 0.67304968, 0.50893644],
    [0.15014782, 0.67663139, 0.50658890], [0.15389405, 0.68020343, 0.50417217],
    [0.15785146, 0.68376525, 0.50168574], [0.16201598, 0.68731632, 0.49912906],
    [0.16638320, 0.69085611, 0.49650163], [0.17094840, 0.69438405, 0.49380294],
    [0.17570671, 0.69789960, 0.49103252], [0.18065314, 0.70140222, 0.48818938],
    [0.18578266, 0.70489133, 0.48527326], [0.19109018, 0.70836635, 0.48228395],
    [0.19657063, 0.71182668, 0.47922108], [0.20221902, 0.71527175, 0.47608431],
    [0.20803045, 0.71870095, 0.47287330], [0.21400015, 0.72211371, 0.46958774],
    [0.22012381, 0.72550945, 0.46622638], [0.22639690, 0.72888753, 0.46278934],
    [0.23281498, 0.73224735, 0.45927675], [0.23937390, 0.73558828, 0.45568838],
    [0.24606968, 0.73890972, 0.45202405], [0.25289851, 0.74221104, 0.44828355],
    [0.25985676, 0.74549162, 0.44446673], [0.26694127, 0.74875084, 0.44057284],
    [0.27414922, 0.75198807, 0.43660090], [0.28147681, 0.75520266, 0.43255207],
    [0.28892102, 0.75839399, 0.42842626], [0.29647899, 0.76156142, 0.42422341],
    [0.30414796, 0.76470433, 0.41994346], [0.31192534, 0.76782207, 0.41558638],
    [0.31980860, 0.77091403, 0.41115215], [0.32779580, 0.77397953, 0.40664011],
    [0.33588539, 0.77701790, 0.40204917], [0.34407411, 0.78002855, 0.39738103],
    [0.35235985, 0.78301086, 0.39263579], [0.36074053, 0.78596419, 0.38781353],
    [0.36921420, 0.78888793, 0.38291438], [0.37777892, 0.79178146, 0.37793850],
    [0.38643282, 0.79464415, 0.37288606], [0.39517408, 0.79747541, 0.36775726],
    [0.40400101, 0.80027461, 0.36255223], [0.41291350, 0.80304099, 0.35726893],
    [0.42190813, 0.80577412, 0.35191009], [0.43098317, 0.80847343, 0.34647607],
    [0.44013691, 0.81113836, 0.34096730], [0.44936763, 0.81376835, 0.33538426],
    [0.45867362, 0.81636288, 0.32972749], [0.46805314, 0.81892143, 0.32399761],
    [0.47750446, 0.82144351, 0.31819529], [0.48702580, 0.82392862, 0.31232133],
    [0.49661536, 0.82637633, 0.30637661], [0.50627130, 0.82878621, 0.30036211],
    [0.51599182, 0.83115784, 0.29427888], [0.52577622, 0.83349064, 0.28812650],
    [0.53562110, 0.83578452, 0.28190832], [0.54552440, 0.83803918, 0.27562602],
    [0.55548397, 0.84025437, 0.26928147], [0.56549760, 0.84242990, 0.26287683],
    [0.57556297, 0.84456561, 0.25641457], [0.58567772, 0.84666139, 0.24989748],
    [0.59583934, 0.84871722, 0.24332878], [0.60604528, 0.85073310, 0.23671214],
    [0.61629283, 0.85270912, 0.23005179], [0.62657923, 0.85464543, 0.22335258],
    [0.63690157, 0.85654226, 0.21662012], [0.64725685, 0.85839991, 0.20986086],
    [0.65764197, 0.86021878, 0.20308229], [0.66805369, 0.86199932, 0.19629307],
    [0.67848868, 0.86374211, 0.18950326], [0.68894351, 0.86544779, 0.18272455],
    [0.69941463, 0.86711711, 0.17597055], [0.70989842, 0.86875092, 0.16925712],
    [0.72039115, 0.87035015, 0.16260273], [0.73088902, 0.87191584, 0.15602894],
    [0.74138803, 0.87344918, 0.14956101], [0.75188414, 0.87495143, 0.14322828],
    [0.76237342, 0.87642392, 0.13706449], [0.77285183, 0.87786808, 0.13110864],
    [0.78331535, 0.87928545, 0.12540538], [0.79375994, 0.88067763, 0.12000532],
    [0.80418159, 0.88204632, 0.11496505], [0.81457634, 0.88339329, 0.11034678],
    [0.82494028, 0.88472036, 0.10621724], [0.83526959, 0.88602943, 0.10264590],
    [0.84556056, 0.88732243, 0.09970219], [0.85580960, 0.88860134, 0.09745186],
    [0.86601325, 0.88986815, 0.09595277], [0.87616824, 0.89112487, 0.09525046],
    [0.88627146, 0.89237353, 0.09537439], [0.89632002, 0.89361614, 0.09633538],
    [0.90631121, 0.89485467, 0.09812496], [0.91624212, 0.89609127, 0.10071680],
    [0.92610579, 0.89732977, 0.10407067], [0.93590444, 0.89857040, 0.10813094],
    [0.94563626, 0.89981500, 0.11283773], [0.95529972, 0.90106534, 0.11812832],
    [0.96489353, 0.90232311, 0.12394051], [0.97441665, 0.90358991, 0.13021494],
    [0.98386829, 0.90486726, 0.13689671], [0.99324789, 0.90615657, 0.14393620]
]

viridis = LinearSegmentedColormap.from_list('viridis', viridis_data)
mpl.cm.register_cmap('viridis', viridis)
viridis_data.reverse()
viridis_r = LinearSegmentedColormap.from_list('viridis', viridis_data)
mpl.cm.register_cmap('viridis', viridis)
