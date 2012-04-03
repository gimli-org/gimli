#!/usr/bin/env python

import pygimli as g
from pygimli.mplviewer import drawMesh, drawModel,drawSelectedMeshBoundaries

import pylab as P
import numpy as np

def getaxes():
    fig = P.figure()
    left, width = 0.1, 0.6

    rect1 = [left, 0.75, width, 0.15]
    rect2 = [left, 0.05, width, 0.69]

    ax1 = fig.add_axes(rect1)
    ax2 = fig.add_axes(rect2,sharex=ax1) 
    return ax1,ax2
    
def analyticalCircle2D( spnts, radius, pos, dDensity ):
    '''
        calculate horizontal component of gravimetric field in mGal
        sphere radius = radius in [Meter]
        sphere center at pos at (x, -depth) 
        dDensity in [kg/m^3]
    '''
    depth = -pos[ 1 ] 
    G = 6.67384e-11
    
    gz = g.RVector( len(spnts), 0.0 )
    for i, p in enumerate( spnts ):
        gz[ i ] = G * 2. * P.pi * radius ** 2. * dDensity * depth / ( (p[0]-pos[0])**2.  + depth**2. ) * 1e5
        
    return gz
#def analyticalSphere(  ):

def analyticalSphere( spnts, radius , pos, dDensity ):
    '''
        calculate horizontal component of gravimetric field in mGal
        sphere radius = radius in [Meter]
        sphere center location at (0.0, 0.0, -depth) 
        dDensity in [kg/m^3]
    '''
    depth = -pos[ 2 ] 
    G = 6.67384e-11
    V = 4./3. * np.pi * radius ** 3.
    M = V * dDensity
    
    k = 1e5 * G * M * depth
    
    gz = g.RVector( len(spnts), 0.0 )
    for i, p in enumerate( spnts ):
        #r = (p-pos).abs()
        #gdr = G * M / ( r **2.) * 1e5
        #gdz = gdr * depth / r
        #gz[ i ] = G * V * dDensity * depth / ( r**3.) * 1e5
        gz[ i ] = k / ( (p[0]-pos[0])**2. + depth**2. )**(3./2.)
       
    return gz
#def analyticalSphere(  ):
    
    
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

def calcGCells( pos, mesh, rho, nInt = 0 ):
    '''
    '''
    G = g.RMatrix( len( pos ), mesh.cellCount() )
    rules = g.IntegrationRules()
    
    for i, p in enumerate( pos ):
        for cId, c in enumerate( mesh.cells() ):
            Z = 0.

            if nInt == 0:
                for j in range( c.nodeCount() ):
                
                    A = c.node( j )
                    B = c.node( (j+1)%c.nodeCount() )
                    # negative Z because all cells are numbered counterclockwise
                    Z -= 2.0 * g.lineIntegraldGdz( A.pos() - p , B.pos() - p )
                    #Z += lineIntegralZ( A.pos() - p , B.pos() - p )
            else:

                for j, t in enumerate( rules.abscissa( c, nInt ) ):
                    w = rules.weights( c, nInt )[ j ]
                    Z += 2.0 * w * functor( c.shape().xyz( t ), p )
                #for j, t in enumerate( rules.quaAbscissa( nInt ) ):
                    #w = rules.quaWeights( nInt )[ j ]
                    #Z += 2.0 * w * functor( c.shape().xyz( t ), p )

                Z *= c.jacobianDeterminant()

            # negative Z because all cells are numbered counterclockwise
            G[ i ][ c.id() ] = Z
    
    return G * rho * 6.67384e-11 * 1e5, G

def test2d():    
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
    
    gzC, GC = calcGCells( pnts, mesh, rho, 1 )
    #gzC = g.calcGCells( spnts , mesh, rho, 1 )
    print "calcGCells",  swatch.duration( True )
    #gzB, GB = calcGBounds( pnts, mesh, rho )
    gzB = g.calcGBounds( spnts , mesh, rho )
    print "calcGBounds", swatch.duration( True )

    gZ_Mesh = gzC

    
    ax1, ax2 = getaxes()

    # sphere analytical solution
    gAna  = analyticalCircle2D( spnts, radius = 2.0, pos = g.RVector3( 0.0, -5.0 ), dDensity = 2000 )
    gAna2 = analyticalCircle2D( spnts, radius = 2.0, pos = g.RVector3( 5.0, -5.0 ), dDensity = 2000 )
    
    gAna = gAna + gAna2

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
    radius = 2.
    depth = 5.
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

#def test2d


def lineIntegralZ3D( p1, p2, n ):
    '''
        Holstein1999
    '''
    
    t = ( p2 - p1 ) / ( p2 - p1 ).abs()
    v = p1.dot( n )
    
    h = p1.dot( t.cross( n ) ) 
    
    r0 = np.sqrt( v**2 + h**2 )
    l1 = p1.dot( t )
    l2 = p2.dot( t )
    
    r1  = np.sqrt( v**2 + h**2 + l1**2 )
    r2  = np.sqrt( v**2 + h**2 + l2**2 )
    
    L = l2 - l1 # edgelength
    
    # vertex method
    ch = h * np.sign( l1 + l2 ) * np.log( ( r2 + abs( l2 ) ) / ( r1 + abs( l1 ) ) )  #(51)
    atanLambda = np.arctan2( 2.0 * h * L, ( r1 + r2 + L ) * ( r1 + r2 - L ) + 2.* abs(v) * ( r2 + r1 ) ) #(39)
    
    # line method
    #A = abs( L ) / ( r2 + r1 ) # (46)
    #S = 0.5 * ( r1 + r2 - abs( L ) * A ) # (47)
    
    #atanLambda = np.arctan2(  ( h * A ), ( abs(v) * S ) ) #(48)
    #ch = 2. * h * arctanh( L ) # (55)
        
    return ch - 2. * abs( v ) * atanLambda # (50)

    
def createHolstein1999Model( ):
    mesh = g.Mesh( 3 )
    
    mesh.createNode( g.RVector3(  0,  0, 0 ) ) # dummy
    mesh.createNode( g.RVector3(  10,  10, -12 ) )
    mesh.createNode( g.RVector3(  10, -10, -12 ) )
    mesh.createNode( g.RVector3( -10, -10, -12 ) )
    mesh.createNode( g.RVector3( -10,  10, -12 ) )
    mesh.createNode( g.RVector3( -20,  30, -12 ) )
    mesh.createNode( g.RVector3(  30,  30, -12 ) )
    mesh.createNode( g.RVector3(  20,  20, -22 ) )
    mesh.createNode( g.RVector3(  20, -30, -22 ) )
    mesh.createNode( g.RVector3( -20, -30, -22 ) )
    mesh.createNode( g.RVector3( -20,  20, -22 ) )

    facets = [ [1, 6, 5, 4, 3, 2],
              [1, 2, 8, 7],
              [2, 3, 9, 8],
              [3, 4, 10, 9],
              [4, 5, 10],
              [5, 6, 7, 10],
              [1, 7, 6],
              [7, 8, 9, 10] ]
              
    return mesh, facets
# def createHolstein1999Model( ... )
    
def test3d():
    '''
    '''
    
    # test1. 
    mesh = g.Mesh( 3 )
    mesh.importSTL( 'sphere.stl' )
    mesh.scale( g.RVector3( 4.0, 4.0, 4.0 ) )
    mesh.translate( g.RVector3( 0.0, 0.0, -5.0 ) )
    print mesh
    mesh.exportVTK('sphere.vtk' )
    
    alpha = 24
    v = 8.2
    epsilon = 1e-12
    gamma = ( 100 * epsilon )** 1./v
    d = alpha / ( gamma * np.sqrt( 2 ) )
    d = 0
    p = g.RVector3( d, d, 0.0 )
    
    x = P.arange( -15, 15, 1. );
    spnts = g.stdVectorRVector3()
    
    for i in x:
        spnts.append( g.RVector3( i, 0.000 ) )
        
    gz = g.RVector( len(x), 0.0 )
    for i, p in enumerate( spnts ):   
        
        for face in mesh.boundaries():
            Z = 0
            norm = face.norm()
            # cos( n, z )
                
            for j in range( face.nodeCount( ) ):
                A = face.node( j )
                B = face.node( (j+1)%face.nodeCount( ) )
                Z += lineIntegralZ3D( A.pos() -p , B.pos() -p, norm )
          
            gz[ i ] += Z
    
        
    gz = gz * 2000 * 6.67384e-11 * 1e5
    
    gAna = analyticalSphere( spnts, radius = 2.0, pos = g.RVector3( 0.0, 0.0, -5.0 ), dDensity = 2000 )
    
    ax1, ax2 = getaxes()
    ax1.plot( x, gAna, '-x', label = 'Analytical Sphere' )
    
    ax1.plot( x, gz, label = 'polyhedral' )
    
    ax1.legend()
    ax2.set_xlim( [ x[0], x[-1] ] )
    ax2.grid()
    
    
# def test3d()


def functor( x, p ):
    '''
        
    '''
    
    r = x.dist( p )
    
    return (p[1]-x[1]) / (r**2.)
# def rz( pos, P )

mesh = g.Mesh(2)
#mesh.load('mesh/world2d.bms')
A = mesh.createNode( g.RVector3( -1., -1. ) )
B = mesh.createNode( g.RVector3( -2., -1. ) )
C = mesh.createNode( g.RVector3( -2., -2. ) )
D = mesh.createNode( g.RVector3( -1., -2. ) )

mesh.createTriangle( A, B, C)
#mesh.createQuadrangle( A, B, C, D)
mesh.scale( g.RVector3( 3., 3. ) )
#mesh.createTriangle( A, B, D)

print mesh.cellSizes()

x = np.arange( -10, 10, 1. );
rho = g.RVector( len( mesh.cellAttributes() ), 1. ) * 2000.0
print rho
pnts = g.stdVectorRVector3()

for i in x:
    pnts.append( g.RVector3( i, 0.0001 ) )

gzNum = []
gzNum.append( calcGCells( pnts, mesh, rho, 0 )[0] )

#P.plot(x, gzNum[0], label = str(0) )

for i in range( 1, 10 ):
    gzNum.append( calcGCells( pnts, mesh, rho, i )[0] )

    err = g.abs(gzNum[i]/gzNum[0]-1.)*100.
    P.semilogy(x, err, label = str(i) )
    
    P.plot(x, gzNum[i], label = str(i) )

P.legend()
#P.plot(x, p1 )


test2d()
#test3d()

P.show()
