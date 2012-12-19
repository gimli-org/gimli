# -*- coding: utf-8 -*-
''' pygimli base functions '''

import pylab as P
import numpy as N
import pygimli as g
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.cm import jet
from pygimli.mplviewer import setMappableData
#from math import sqrt, floor, ceil

def gmat2numpy(mat):
    """convert pygimli matrix into numpy.array."""
    nmat = N.zeros( ( len( mat ), len( mat[ 0 ] ) ) )
    for i, row in enumerate( mat ): 
        nmat[i] = row
    return nmat

def numpy2gmat(nmat):
    """convert numpy.array into pygimli RMatrix."""
    gmat = g.RMatrix()
    for arr in nmat:
        gmat.push_back( g.asvector( arr ) )
    return gmat

def rndig(a, ndig=3):
    """round float using a number of counting digits."""
    if N.abs( a ) < 1e-4:
        return a
    else:
        return N.around( a, ndig - int( N.ceil( N.log10( N.abs( a ) + 1e-4 ) ) ) )

def num2str(a):
    s=[]
    for ai in a:
        s.append( '%g' % rndig(ai) )
    return s

def inthist( a, vals, bins=None, islog=False ):
    if bins is None:
        bins = N.min( ( N.round( len(a) / 20 ) , 10 ) )

    if islog:
        hists, edges = N.histogram( N.log( a ), bins=bins )
    else:        
        hists, edges = N.histogram( a, bins=bins )

    cums = N.cumsum( N.hstack( (0., hists) ) ) / N.sum( hists ) * 100.
    out = N.interp( vals, cums, edges )
    if islog:
        return N.exp( out )
    else:    
        return out

def interperc( a, trimval=3.0, islog=False, bins=None ):
    return inthist( a, N.array( [ trimval, 100. - trimval ] ), bins=bins, islog=islog )

def jetmap(m=64):
    """ jet color map """
    n = int( P.ceil(m/4) )
    u = P.hstack( ( ( P.arange( 1. * n ) + 1 ) / n, P.ones( n - 1 ),
                   P.arange( 1. * n, 0, -1 ) / n ) )
    g1 = P.arange( len(u), dtype=int ) + n / 2
    r1 = g1 + n
    b1 = g1 - n
    gg = g1[ g1 < n * 4 ]
    rr = r1[ r1 < n * 4 ]
    bb = b1[ b1 >= 0 ]
    J = P.zeros( ( n * 4, 3 ) )
    J[rr, 0] = u[:len(rr)]
    J[gg, 1] = u[:len(gg)]
    J[bb, 2] = u[-len(bb):]
    return J

def showmymatrix(A,x,y,dx=2,dy=1,xlab=None,ylab=None,cbar=None):
    P.imshow(A,interpolation='nearest')
    P.xticks(N.arange(0,len(x),dx),["%g" % rndig(xi,2) for xi in x]) #,b
    P.yticks(N.arange(0,len(y),dy),["%g" % rndig(yi,2) for yi in y]) #,a
    P.ylim((len(y)-0.5,-0.5))
    if xlab is not None: P.xlabel(xlab)
    if ylab is not None: P.ylabel(ylab)
    P.axis('auto') 
    if cbar is not None: P.colorbar(orientation=cbar)
    return
    
def draw1dmodel(x, thk=None, xlab=None, zlab="z in m", islog=True, fs=14, z0=0, **kwargs):
    """draw 1d block model defined by value and thickness vectors."""
#    if xlab is None: 
#        xlab = "$\\rho$ in $\\Omega$m"

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
    li = []
    if islog:
        li = P.semilogx( px, pz, **kwargs )
    else:        
        li = P.plot( px, pz, **kwargs )
        
    P.gca().xaxis.set_label_position('top')

    locs = P.xticks()[0]
    if len( locs ) < 2:
        locs = N.hstack( ( min(x), locs, max( x ) ) )
    elif len( locs ) < 5:
        locs[0] = max( locs[0], min(x) )
        locs[-1] = min( locs[-1], max(x) )
    
    a = []
    for l in locs:
        a.append( '%g' % rndig(l) )

    P.xticks( locs, a, fontsize=fs )
    P.yticks(fontsize = fs)

    P.xlim( ( N.min(x) * 0.9, N.max(x) * 1.1 ) )
    P.ylim( ( max(z1) * 1.15 , 0. ) )
    if xlab is not None: P.xlabel(xlab,fontsize=fs)
    if zlab is not None: P.ylabel(zlab,fontsize=fs)
    P.grid(which='both')
    #P.show()
    return li

def draw1dmodelLU( x, xL, xU, thk=None, **kwargs ):
    """ draw 1d model with lower and upper bounds """
    draw1dmodel(x,thk,color='red',**kwargs)
    for i in range( len(x) ):
        x1 = x.copy()
        x1[i] = xL[i]
        draw1dmodel(x1,thk,color='blue')
        x1[i] = xU[i]
        draw1dmodel(x1,thk,color='blue')

    li = draw1dmodel(x,thk,color='red',**kwargs)
    P.xlim( ( min( xL ) * 0.9, max( xU ) * 1.1 ) )
    return li

def showStitchedModels(models, x=None, cmin=None, cmax=None,
                       islog=True, title=None):
    """ show several 1d block models as (stitched) section """
    if x is None:
        x = P.arange( len(models) )
        
    nlay = int( P.floor( ( len(models[0]) + 1 ) / 2. ) )
        
    ax = P.gcf().add_subplot(111)
    ax.cla()

    dxmed2 = P.median( P.diff(x) ) / 2.
    vals = P.zeros( (len(models),nlay) )
    patches = []
    maxz = 0.
    for i, imod in enumerate( models ):
        if isinstance( imod, g.RVector ):
            vals[i,:] = imod(nlay-1,2*nlay-1)
            thk = P.asarray( imod( 0, nlay-1 ) )
        else:
            vals[i,:] = imod[nlay-1:2*nlay-1]
            thk = imod[:nlay-1]
        
        thk = P.hstack( (thk, thk[-1]) )
        z = P.hstack( (0., P.cumsum(thk)) )
        maxz = max( maxz, z[-1] )

        for j in range( nlay ):    
            rect = Rectangle( ( x[i]-dxmed2, z[j] ), dxmed2*2, thk[j] )
            patches.append( rect )

    p = PatchCollection(patches, cmap=jet, linewidths=0)
    
    if cmin is not None:
        p.set_clim( cmin, cmax )
    
    #p.set_array( P.log10( vals.ravel() ) )
    setMappableData( p, vals.ravel(), logScale = True )
    ax.add_collection(p)
    
    ax.set_ylim( ( maxz, 0. ) )
    ax.set_xlim( ( min(x)-dxmed2, max(x)+dxmed2 ) )
    if title is not None:
        P.title( title )

    P.colorbar(p, orientation='horizontal',aspect=50,pad=0.1)

    P.draw()
    return

def showStitchedModelsOld(models, x=None, cmin=None, cmax=None,
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
    P.plot( x, P.zeros( len( x ) ), 'k.' )
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
