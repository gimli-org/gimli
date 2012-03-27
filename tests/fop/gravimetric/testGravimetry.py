#!/usr/bin/env python

import pygimli as g
from pygimli.mplviewer import drawMesh, drawModel,drawSelectedMeshBoundaries

import pylab as P
import numpy as np

def lineIntegralZ( p1, p2 ):
    '''
        WonBevis1987
    '''
   
    x1 = p1[ 0 ]; z1 = p1[ 1 ]
    x2 = p2[ 0 ]; z2 = p2[ 1 ]
    
    if x1 == 0. and z1 == 0.:
        return 0.0
    if x2 == 0. and z2 == 0.:
        return 0.0
        
    theta1 = np.arctan2( z1, x1 )
    theta2 = np.arctan2( z2, x2 )

    r1 = np.sqrt( x1*x1 + z1*z1 )
    r2 = np.sqrt( x2*x2 + z2*z2 )
    
    Z = 0.0
    
    if np.sign( z1 ) != np.sign( z2 ):
        if ( x1*z2 < x2*z1 ) and z2 >=0.0: 
            theta1 = theta1 + 2.* np.pi
        
        if ( x1*z2 > x2*z1 ) and z1 >=0.0: 
            theta2 = theta2 + 2.* np.pi
        
        if x1*z2 == x2*z1:
            Z = 0.
            
    if x1 == x2: 
        #case 3
        Z = x1 * np.log( r2 / r1 )
                
    else: #default
        B = ( z2 - z1 ) / ( x2 - x1 )
        A = ( ( x2 - x1 ) * ( x1 * z2 - x2 * z1 ) ) / ( ( x2 - x1 )**2.0 + ( z2 - z1 )**2.0  )

        Z = A * ( ( theta1 - theta2 ) + B * np.log( r2 / r1 ) )
    
    return Z

def calcPolydgdz( pnts, poly, density ):
    '''
       Calculate gravimetric response at given points for a polygon with relative density change. 
       pnts must be numbered clockwise. Else change the sign of the result.
       Return values are in mGal.
    '''
    gz = g.RVector( len( pnts ), 0.0 )
    
    for i, p in enumerate( pnts ):
        for j in range( len( poly ) ):
            a = poly[ j ]
            b = poly[ (j+1)%len( poly ) ]
            gz[ i ] += g.lineIntegraldGdz( a - p , b - p )
            #gz[ i ] += lineIntegralZ( a - p , b - p )
    return gz * density * 2.0 * 6.67384e-11 * 1e5
# def calcPolydgdz()
    
    
def calcGBounds( pos, mesh, rho ):
    '''
    '''
    G = g.RMatrix( len( pos ), mesh.cellCount() )
        
    for i, p in enumerate( pos ):
        for cId, b in enumerate( mesh.boundaries() ):
            A = b.node( 0 )
            B = b.node( 1 )
            #Z = lineIntegralZ( A.pos() - p , B.pos() - p )
            Z = g.lineIntegraldGdz( A.pos() - p , B.pos() - p )
            
            if b.leftCell(): 
                G[ i ][ b.leftCell().id() ] = G[ i ][ b.leftCell().id() ] - Z
            if b.rightCell(): 
                G[ i ][ b.rightCell().id() ] = G[ i ][ b.rightCell().id() ] + Z
   
    return G * rho * 2.0 * 6.67384e-11 * 1e5, G

def calcGCells( pos, mesh, rho ):
    '''
    '''
    G = g.RMatrix( len( pos ), mesh.cellCount() )
    
    for i, p in enumerate( pos ):
        for cId, c in enumerate( mesh.cells() ):
            Z = 0
            for j in range( c.nodeCount() ):
                
                A = c.node( j )
                B = c.node( (j+1)%c.nodeCount() )
                Z += g.lineIntegraldGdz( A.pos() - p , B.pos() - p )
                #Z += lineIntegralZ( A.pos() - p , B.pos() - p )
                
            # negative Z because all cells are numbered counterclockwise
            G[ i ][ c.id() ] = -Z
    
    return G * rho * 2.0 * 6.67384e-11 * 1e5, G

mesh = g.Mesh( 'mesh/world2d.bms' )
print mesh

xMin = mesh.boundingBox( ).min()[0]
yMax = mesh.boundingBox( ).max()[0]
x = P.arange( xMin, yMax, 1. );

mesh.createNeighbourInfos()
rho = g.RVector( len( mesh.cellAttributes() ), 1. ) * 2000.0 
rho.setVal( 0.0, g.find( mesh.cellAttributes() == 1.0 ) )

swatch = g.Stopwatch( True )
pnts = []
spnts = g.stdVectorRVector3()
for i in x:
    pnts.append( g.RVector3( i, 0.0001 ) )
    spnts.append( g.RVector3( i, 0.0001 ) )
    
#gzC, GC = calcGCells( pnts, mesh, rho )
gzC = g.calcGCells( spnts , mesh, rho )
print "calcGCells",  swatch.duration( True )
#gzB, GB = calcGBounds( pnts, mesh, rho )
gzB = g.calcGBounds( spnts , mesh, rho )
print "calcGBounds", swatch.duration( True )

gZ_Mesh = gzB

fig = P.figure()
left, width = 0.1, 0.6

rect1 = [left, 0.75, width, 0.15]
rect2 = [left, 0.05, width, 0.69]

ax1 = fig.add_axes(rect1)
ax2 = fig.add_axes(rect2,sharex=ax1) 


# sphere analytical solution
radius = 2.0
dDensity = 2000.0
depth = 5.0 
G = 6.67384e-11
gAna = G * 2. * P.pi * radius **2. * dDensity * depth / ( x**2.  + depth**2. )
gAna2 = G * 2. * P.pi * radius **2. * dDensity * depth / ( (x-5.)**2.  + depth**2. )
gAna = gAna+ gAna2
gAna *= 1e5 # mGal

ax1.plot( x, gAna, '-x', label= 'Analytisch' )
ax1.plot( x, gZ_Mesh, label= 'WonBevis1987-mesh' )

print gAna / gZ_Mesh

#rho=GB[0]/mesh.cellSizes()

gci = drawModel( ax2, mesh, rho )

drawSelectedMeshBoundaries( ax2, mesh.findBoundaryByMarker( 0 )
                                , color = ( 1.0, 1.0, 1.0, 1.0 )
                                , linewidth = 0.3 )
drawSelectedMeshBoundaries( ax2, mesh.findBoundaryByMarker( 1 )
                                , color = ( 1.0, 1.0, 1.0, 1.0 )
                                , linewidth = 0.3 )
       
# sphere polygone
poly1 = g.stdVectorRVector3()
poly2 = g.stdVectorRVector3()
nSegment=124
for i in range( nSegment ):
    xp = np.sin( (i+1) * ( 2. * np.pi ) / nSegment )
    yp = np.cos( (i+1) * ( 2. * np.pi ) / nSegment )
    poly1.append( g.RVector3( xp * radius, yp * radius - depth ) )
    poly2.append( g.RVector3( xp * radius + 5., yp * radius - depth ) )

gZ_Poly  = calcPolydgdz( spnts, poly1, 2000 )
gZ_Poly += calcPolydgdz( spnts, poly2, 2000 )

ax1.plot( x, gZ_Poly, label= 'WonBevis1987-Poly' )
ax2.plot( g.x( poly1 ), g.y( poly1 ), color = 'red' )
ax2.plot( g.x( poly2 ), g.y( poly2 ), color = 'red' )

ax2.plot( g.x( spnts ), g.y( spnts ), marker = 'x', color = 'black' )

# test some special case
for i, p in enumerate( poly1 ):
    poly1[i] = g.RVector3( poly1[i] - g.RVector3( 5.0, -6. ) )
ax2.plot( g.x( poly1 ), g.y( poly1 ), color = 'green' )

gz  = calcPolydgdz( spnts, poly1, 2000 )
ax1.plot( x, gz, label= 'Special Case', color = 'green' )
ax1.set_ylabel( 'dg/dz [mGal]' )
ax2.set_ylabel( 'Tiefe [m]' )

ax1.legend()
ax2.set_xlim( [ x[0], x[-1] ] )
ax2.grid()



P.show()

