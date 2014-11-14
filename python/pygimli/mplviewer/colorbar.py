# -*- coding: utf-8 -*-

import pygimli as pg

import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
import matplotlib.ticker as ticker
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from mpl_toolkits.axes_grid1 import make_axes_locatable

import math

cdict = {'red': ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)),
         'green': ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)),
         'blue': ((0.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0))}

blueRedCMap = mpl.colors.LinearSegmentedColormap('my_colormap', cdict, 256)

def autolevel(z, N, logscale=None):
    """
        Create N levels for the data array z based on matplotlib ticker. 
    """
    locator = None
    if logscale:
        locator = ticker.LogLocator()
    else:
        locator = ticker.MaxNLocator(N + 1)
    
    zmin = min(z)
    zmax = max(z)
       
    return locator.tick_values(zmin, zmax)
        
    # For line contours, drop levels outside the data range.
    return lev[(lev > zmin) & (lev < zmax)]
# def autolevel()

def cmapFromName(cmapname, ncols=256, bad=[1.0, 1.0, 1.0, 0.0]):
    """
    """
    cmap = mpl.cm.get_cmap('jet', ncols)
    
    if cmapname is not None:
        if cmapname == 'b2r':
            cmap = mpl.colors.LinearSegmentedColormap(
                'my_colormap', cdict, ncols)
        else:
            try:
                cmap = mpl.cm.get_cmap(cmapname, ncols)
            except:
                print(("could not retrieve colormap ", cmapname))

    cmap.set_bad(bad)
    return cmap


def findAndMaskBestClim(dataIn, cMin=None, cMax=None,
                        dropColLimitsPerc=5, logScale=True):
    """What is this."""
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

    (Nhist, xHist) = np.histogram(data, bins=100)

    if not cMin:
        cMin = xHist[dropColLimitsPerc]
        if (logScale):
            cMin = pow(10.0, cMin)

    if not cMax:
        cMax = xHist[100 - dropColLimitsPerc]
        if (logScale):
            cMax = pow(10.0, cMax)

    if (logScale):
        data = pow(10.0, data)

    data[np.where(data < cMin)] = cMin
    data[np.where(data > cMax)] = cMax

    return data, cMin, cMax


def createLogLevs(vMin, vMax, nLevs):
    vMinLog = np.log10(vMin)
    vMaxLog = np.log10(vMax)

    lev_exp = list(range(0, nLevs))
    lev_exp[0] = vMinLog
    dxLog = (vMaxLog - vMinLog) / (nLevs - 1)

    for i in range(nLevs - 1):
        lev_exp[i + 1] = lev_exp[i] + dxLog

    levs = np.power(10, lev_exp)

    return levs


def createLinLevs(vMin, vMax, nLevs):
    levs = list(range(0, nLevs))
    levs[0] = vMin
    dx = (float(vMax) - float(vMin)) / (nLevs - 1)

    for i in range(nLevs - 1):
        levs[i + 1] = levs[i] + dx

    return levs


def createColorbar2(patches, cMin=None, cMax=None,
                    nLevs=5, label=None, orientation='horizontal'):
    
    DEPRECATED
    cbarTarget = plt
    
    if hasattr(patches, 'ax'):
        cbarTarget = patches.ax

    cax = mpl.colorbar.make_axes(cbarTarget,
                                 orientation=orientation,
                                 aspect=50)

    # print cax
    cbar = mpl.colorbar.Colorbar(cax[0], patches,
                                 orientation=orientation)

#    if cMin is None:
#        cMin= patches.zmin
#    if cMax is None:
#        cMax= patches.zmax

    # setCbarLevels( cbar, cMin, cMax, nLevs )

    if label is not None:
        cbar.set_label(label)

    return cbar


def createColorbar(patches, cMin=None, cMax=None, nLevs=5,
                   label=None, orientation='horizontal', *args, **kwargs):
    cbarTarget = plt
    cax = None
    divider = None
    #if hasattr(patches, 'figure'):
        #cbarTarget = patches.figure

    #print( patches)

    if hasattr(patches, 'ax'):
        divider = make_axes_locatable(patches.ax)
    elif hasattr(patches, 'get_axes'):
        divider = make_axes_locatable(patches.get_axes())
    
    if divider:
        if orientation == 'horizontal':
            cax = divider.append_axes("bottom", size=0.25, pad=0.65)
        else:
            cax = divider.append_axes("right", size=0.25, pad=0.05)
            #cax = divider.append_axes("right", size="5%", pad=0.05)
            #cbar3 = plt.colorbar(im3, cax=cax3)
        
    #print(patches,  cax)
    cbar = cbarTarget.colorbar(patches, cax=cax,
                               orientation=orientation)

    setCbarLevels(cbar, cMin, cMax, nLevs)

    if label is not None:
        cbar.set_label(label)

    return cbar


def setCbarLevels(cbar, cMin=None, cMax=None, nLevs=5):

    #print "setCbarLevels", cMin, cMax

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
        cbarLevels = createLogLevs(cMin, cMax, nLevs)
    else:
        cbarLevels = createLinLevs(cMin, cMax, nLevs)

    #print cbarLevels
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
        cbar.mappable.set_clim(cMin + 1e-6, cMax)

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
    
    data = dataIn
    
    if type(data) != np.ma.core.MaskedArray:
        data = np.array(dataIn)
    
    if logScale and data.min() <= 0:
        data = np.ma.masked_array(data, data <= 0.0)

    # set bad value color to white
    if mappable.get_cmap() is not None:
        mappable.get_cmap().set_bad([1.0, 1.0, 1.0, 0.0 ])

    if not cMin:
        cMin = data.min()
    if not cMax:
        cMax = data.max()

    #print("set mappable data, log: ", logScale, "cmin: ", cMin, "cmax: ", cMax)

    if cMin > 0.0 and logScale:
        mappable.set_norm(mpl.colors.LogNorm())
    else:
        mappable.set_norm(mpl.colors.Normalize())

    mappable.set_array(data)
    #mappable.set_level(10)
    mappable.set_clim(cMin, cMax)

def addCoverageAlpha(patches, coverage, dropTolerance=0.4):
    """
        Add alpha values to the colors of a polygon collection.
    
    Parameters
    ----------

    patches : 2D mpl mappable 
        
    coverage : array
        coverage values. Maximum coverage mean no opaqueness.

    dropTolerance : float
        minimum coverage
    
    Usage
    -----
    
        
    """
    
    patches.set_antialiaseds(True)
    patches.set_linewidth(0.000)

    # generate individual color values here
    patches.update_scalarmappable()

    cols = patches.get_facecolor()

    C = np.asarray(coverage)
    #print(np.min(C), np.max(C))
            
    if (np.min(C) < 0.) | (np.max(C) > 1.) | (np.max(C) < 0.5):
        
        nn, hh = np.histogram(C, 50)
        nnn = nn.cumsum(axis = 0) / float(len(C))
        
        #print("min-max nnn ", min(nnn), max(nnn))
        mi = hh[min(np.where(nnn > 0.02)[0])]

        if min(nnn) > 0.4:
            ma = max(C)
        else:
            ma = hh[max(np.where(nnn < 0.4)[0])]

            #mi = hh[min(np.where(nnn > 0.2)[0])]
            #ma = hh[max(np.where(nnn < 0.7)[0])]
            C = (C - mi) / (ma - mi)
            C[np.where(C < 0.)] = 0.0
            C[np.where(C > 1.)] = 1.0

            # add alpha value to the color values
            cols[:, 3] = C

            patches._facecolors = cols
            #patches._edgecolor = 'None'

            # delete patch data to avoid automatically rewrite of _facecolors
            patches._A = None

# def addCoverageAlpha(...)
