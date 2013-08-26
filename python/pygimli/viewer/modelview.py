# -*- coding: utf-8 -*-
''' pygimli model viewer functions '''

import pylab as P
import numpy as N
import pygimli as g
from matplotlib.patches import Rectangle
#from math import sqrt, floor, ceil

def showmymatrix(A, x, y, dx=2, dy=1, xlab=None, ylab=None, cbar=None):
    '''
        Pls. insert short Docu here
    '''
    P.imshow(A, interpolation='nearest')
    P.xticks(N.arange(0, len(x), dx), ["%g" % rndig(xi, 2) for xi in x]) #,b
    P.yticks(N.arange(0, len(y), dy), ["%g" % rndig(yi, 2) for yi in y]) #,a
    P.ylim((len(y) - 0.5, - 0.5))
    
    if xlab is not None: P.xlabel(xlab)
    if ylab is not None: P.ylabel(ylab)
    
    P.axis('auto') 
    
    if cbar is not None: P.colorbar(orientation=cbar)
    
    return
    
def drawModel1D( ax, thickness, values,  plotfunction = 'plot', xlabel = 'Resistivity $[\Omega$ m$]$'
                , *args, **kwargs ):
    '''
        Draw 1d block model into axis ax defined by values and thickness vectors using plotfunction
    '''
    
    nLayers = len( thickness ) + 1
    px = P.zeros( nLayers * 2 )
    pz = P.zeros( nLayers * 2 )
    z1 = P.cumsum( thickness )
    
    for i in range( nLayers ):
        px[ 2 * i ]     = values[ i ]
        px[ 2 * i + 1 ] = values[ i ]
        
        if i == nLayers - 1:
            pz[ 2 * i + 1 ]  = z1[ i - 1 ] * 1.2
            pass
        else:
            pz[ 2 * i + 1 ] = z1[ i ]
            pz[ 2 * i + 2 ] = z1[ i ]
        
    if plotfunction == 'loglog' or plotfunction == 'semilogy':
        pz[ 0 ] = thickness[ 0 ] * 0.8
        
    try:
        plot = getattr( ax, plotfunction )
        plot( px, pz, *args, **kwargs )
    except Exception as e:
        print e
    
    ax.set_ylabel( 'Depth [m]' )
    ax.set_xlabel( xlabel )
    ax.set_ylim( pz[-1], pz[0] )
    ax.grid()
    return ax
    
#def draw1dmodel(... )
    
def draw1dmodel(x, thk=None, xlab=None, zlab="z in m", islog=True, z0=0):
    print "STYLE_WARNING!!!!!!!!!!! don't use this call. Use show1dmodel instead."
    show1dmodel(x, thk, xlab, zlab, islog, z0 )
    
def show1dmodel(x, thk=None, xlab=None, zlab="z in m", islog=True, z0=0):
    """draw 1d block model defined by value and thickness vectors."""
    if xlab is None: 
        xlab = "$\\rho$ in $\\Omega$m"

    if thk is None: #gimli blockmodel (thk+x together) given
        nl = int( N.floor( ( len(x) - 1 ) / 2. ) ) + 1
        thk = N.asarray(x)[:nl-1]
        x = N.asarray(x)[nl-1:nl*2-1]
    
    z1 = N.concatenate( ( [0], N.cumsum( thk ) ) ) + z0
    z = N.concatenate( ( z1, [z1[-1] * 1.2] ) )
    nl = len(x) #x.size()
    px = N.zeros( ( nl * 2, 1 ) )
    pz = N.zeros( ( nl * 2 , 1 ) )
    for i in range( nl ):
        px[2*i] = x[i]
        px[2*i+1] = x[i]
        pz[2*i+1] = z[i+1]
        if i < nl - 1:
            pz[2*i+2] = z[i+1]

#    P.cla()
    if islog:
        P.semilogx( px, pz )
    else:        
        P.plot( px, pz )
    P.grid(which='both')
    P.xlim( ( N.min(x) * 0.9, N.max(x) * 1.1 ) )
    P.ylim( ( max(z1) * 1.15 , 0. ) )
    P.xlabel(xlab)
    P.ylabel(zlab)
    P.show()
    return

def showStitchedModels(models, x=None, cmin=None, cmax=None,
                       islog=True, title=None):
    """ show several 1d block models as (stitched) section """
    if x is None:
        x = P.arange( len(models) )
        
    nlay = int( P.floor( ( len(models[0]) - 1 ) / 2. ) ) + 1
    if cmin is None or cmax is None:
        cmin = 1e9
        cmax = 1e-9
        for model in models:
            res = P.asarray(model)[nlay-1:nlay*2-1]
            cmin = min( cmin, min(res) )
            cmax = max( cmax, max(res) )
            
        print "cmin=", cmin, " cmax=", cmax
        
    dx = P.diff(x)
    dx = P.hstack( (dx, dx[-1]) )
    x1 = x - dx/2
    ax = P.gcf().add_subplot(111)
    ax.cla()
    mapsize = 64
    cmap = jetmap( mapsize )
    P.plot( x, x * 0., 'k.' )
    maxz = 0.
    for i, mod in enumerate( models ):
        mod1 = P.asarray(mod)
        res = mod1[nlay-1:]
        if islog:
            res = P.log( res )
            cmi = P.log( cmin )
            cma = P.log( cmax )
        else:
            cmi = cmin
            cma = cmax
            
            
        thk = mod1[:nlay-1]
        thk = P.hstack( (thk, thk[-1]) )
        z = P.hstack( (0., P.cumsum(thk)) )
        maxz = max( maxz, z[-1] )
        nres = ( res - cmi ) / ( cma - cmi )
        cind = N.around( nres * mapsize )
        cind[ cind >= mapsize ] = mapsize - 1
        cind[ cind < 0 ] = 0
        for j in range( len(thk) ):
            fc = cmap[ cind[j], : ]
            rect = Rectangle( ( x1[i], z[j] ), dx[i], thk[j], fc=fc )
            P.gca().add_patch(rect)
    
    ax.set_ylim( ( maxz, 0. ) )
    ax.set_xlim( ( x1[0], x1[-1] + dx[-1] ) )
    if title is not None:
        P.title( title )

    P.draw()
    return

def showfdemsounding(freq, inphase, quadrat, response=None, npl=2):
    """ show FDEM sounding as real(inphase) and imaginary (quadrature)
        fields normalized by the (purely real) free air solution """
    nf = len(freq)
    fig = P.figure(1)
    fig.clf()
    ax1 = fig.add_subplot(1, npl, npl-1)
    P.semilogy( inphase, freq, 'x-' )
    if response is not None:
        P.semilogy( P.asarray(response)[:nf], freq, 'x-' )
        
    P.grid( which='both' )
    ax2 = fig.add_subplot(1, npl, npl)
    P.semilogy( quadrat, freq, 'x-' )
    if response is not None:
        P.semilogy( P.asarray( response )[nf:], freq, 'x-')
    P.grid( which='both' )
    fig.show()
    ax = [ ax1, ax2 ]
    if npl > 2: 
        ax3 = fig.add_subplot(1, npl, 1)
        ax.append( ax3 )

    return ax
