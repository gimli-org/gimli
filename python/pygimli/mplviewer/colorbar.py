# -*- coding: utf-8 -*-

import pygimli as g

import matplotlib as mpl
import matplotlib.ticker as ticker
import matplotlib.colors as colors
import matplotlib.cbook as cbook
import pylab
import numpy as np
from numpy import array, ma
import math

cdict = { 'red': ((0.0, 0.0, 0.0),
                  (0.5, 1.0, 1.0),
                  (1.0, 1.0, 1.0)),
        'green': ((0.0, 0.0, 0.0),
                  (0.5, 1.0, 1.0),
                  (1.0, 0.0, 0.0)),
         'blue': ((0.0, 1.0, 1.0),
                  (0.5, 1.0, 1.0),
                  (1.0, 0.0, 0.0))}

blueRedCMap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)


def findAndMaskBestClim( dataIn, cMin = None, cMax = None, dropColLimitsPerc = 5, logScale = True ):
    if type( dataIn ) == g.RVector:
        data = asarray( dataIn )
    elif type( dataIn ) == list:
        data = array( dataIn )
    else:
        data = array( dataIn )

    if ( min( data ) < 0 ):
        logScale = False;
    if ( logScale ):
        data = np.log10( data )

    (Nhist, xHist) = np.histogram( data, bins = 100 );

    if not cMin:
        cMin = xHist[ dropColLimitsPerc ];
        if ( logScale ):
            cMin = pow( 10.0, cMin )

    if not cMax:
        cMax = xHist[ 100 - dropColLimitsPerc ];
        if ( logScale ):
            cMax = pow( 10.0, cMax )

    if ( logScale ):
        data = pow( 10.0, data )

    data[ np.where( data < cMin ) ] = cMin
    data[ np.where( data > cMax ) ] = cMax

    return data, cMin, cMax

def createLogLevs( vMin, vMax, nLevs ):
    vMinLog = np.log10( vMin )
    vMaxLog = np.log10( vMax )

    lev_exp = range( 0, nLevs )
    lev_exp[ 0 ] = vMinLog
    dxLog = ( vMaxLog - vMinLog ) / ( nLevs - 1 );

    for i in range( nLevs -1):
        lev_exp[ i + 1 ] = lev_exp[ i ] + dxLog

    levs = np.power( 10, lev_exp )

    return levs

def createLinLevs( vMin, vMax, nLevs ):
    levs = range( 0, nLevs )
    levs[ 0 ] = vMin
    dx = ( float(vMax) - float(vMin) ) / ( nLevs - 1 );

    for i in range( nLevs -1):
        levs[ i + 1 ] = levs[ i ] + dx
    
    return levs

def createColorbar2( patches, cMin = None, cMax = None, nLevs = 5, label = None, orientation = 'horizontal' ):
    cbarTarget = pylab
    if hasattr( patches, 'ax' ):
        cbarTarget = patches.ax

    cax = mpl.colorbar.make_axes( cbarTarget
                                    , orientation = orientation
                                    , aspect      = 50 )

    #print cax
    cbar = mpl.colorbar.Colorbar( cax[0],  patches
                        , orientation = orientation
                        #
                        )

#    if cMin is None:
#        cMin= patches.zmin
#    if cMax is None:
#        cMax= patches.zmax

    #setCbarLevels( cbar, cMin, cMax, nLevs )
    
    if label is not None:
        cbar.set_label( label )

    return cbar

def createColorbar( patches, cMin = None, cMax = None, nLevs = 5, label = None, orientation = 'horizontal' ):
    cbarTarget = pylab
    
    if hasattr( patches, 'figure' ):
        cbarTarget = patches.figure

    cbar = cbarTarget.colorbar( patches
                                , orientation = orientation
                                , aspect      = 50
                              )
                    
    setCbarLevels( cbar, cMin, cMax, nLevs )
    
    if label is not None:
        cbar.set_label( label )

    return cbar

def setCbarLevels( cbar, cMin = None, cMax = None, nLevs = 5 ):

    print "setCbarLevels", cMin, cMax
    
    if cMin is None: cMin = cbar.get_clim()[ 0 ]
    if cMax is None: cMax = cbar.get_clim()[ 1 ]

    #cbar.add_lines( patches )
#    if cMin is None and hasattr( patches, 'zmin' ):
#        cMin = patches.zmin
#    if cMax is None and hasattr( patches, 'zmin' ):
#        cMax = patches.zmax
        
#    if cMin is None and hasattr( patches, 'get_clim' ):
#        cMin = patches.get_clim()[0]
#    if cMax is None and hasattr( patches, 'get_clim' ):
#        cMax = patches.get_clim()[1]

    #if isinstance( cbar.mappable.norm, mpl.colors.LogNorm ) and cMin <= 0:
        #cMin = cbar.get_clim()[ 0 ]
        #cMax = cbar.get_clim()[ 1 ]
    #else:
        ##print "lin: cMin", cMin, "cMax", cMax
        #pass
    if cMin == cMax:
        cMin*=0.999
        cMax*=1.001


    norm = None
    if hasattr( cbar, 'mappable' ):
        norm = cbar.mappable.norm
    elif hasattr( cbar, 'norm' ):
        norm = cbar.norm
        cMin = norm.vmin
        cMax = norm.vmax
        
    if isinstance( norm, mpl.colors.LogNorm ):
        cbarLevels = createLogLevs( cMin, cMax, nLevs )
    else:
        cbarLevels = createLinLevs( cMin, cMax, nLevs )


    print cbarLevels
    cbarLevelsString=[]
    for i in cbarLevels:
        if abs(i) == 0.0:
            cbarLevelsString.append( "0" )
        elif abs(i) > 1e4 or abs(i) <= 1e-4:
            cbarLevelsString.append( "%.1e" %i )
        elif abs(i) < 1e-3:
            cbarLevelsString.append( "%.5f" %i )
        elif abs(i) < 1e-2:
            cbarLevelsString.append( "%.4f" %i )
        elif abs(i) < 1e-1:
            cbarLevelsString.append( "%.3f" %i )
        elif abs(i) < 1e0:
            cbarLevelsString.append( "%.2f" %i )
        elif abs(i) < 1e1:
            cbarLevelsString.append( "%.1f" %i )
        else :
            cbarLevelsString.append( "%.0f" %i )

    if hasattr( cbar, 'mappable' ):
        cbar.mappable.set_clim( cMin + 1e-6, cMax )

    #print ticks, cbarLevels, cbarLevelsString
    #if cbar.orientation == 'vertical':
    cbar.set_ticks( cbarLevels )
    cbar.set_ticklabels( cbarLevelsString )

    #print cbar._ticker()

#    else:
#        cbar.ax.set_xticks( ticks )
#        cbar.ax.set_xticklabels( cbarLevelsString )


    cbar.draw_all()

def setMappableData( mappable, dataIn, cMin = None, cMax = None, logScale = False ):
    data = np.asarray( dataIn )
    
    if logScale and data.min() <= 0:
        data = ma.masked_where( data <= 0.0, data )

    # set bad value color to white
    if mappable.get_cmap() is not None:
        mappable.get_cmap().set_bad( [1.0, 1.0, 1.0, 0.0 ] )

    if not cMin: cMin = data.min()
    if not cMax: cMax = data.max()

    #print "set mappable data, log: ", logScale, "cmin: ", cMin, "cmax: ", cMax

    if cMin > 0.0 and logScale:
        mappable.set_norm( mpl.colors.LogNorm() )
    else:
        mappable.set_norm( mpl.colors.Normalize() )

    mappable.set_array( data )
    mappable.set_clim( cMin, cMax)
