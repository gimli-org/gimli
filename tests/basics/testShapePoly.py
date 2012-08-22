#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as g

import pylab as P
from mpl_toolkits.mplot3d import Axes3D

def findShapeFunction( uN, dim, nCoeff, pnts, pascale, serendipity, powCombination ): 

    fop = g.PolynomialModelling( dim, nCoeff, pnts )
    fop.setPascalsStyle( pascale )
    fop.setSerendipityStyle( serendipity )
    fop.setPowCombinationTmp( powCombination )

    fop2 = g.PolynomialModelling( 2, 3, pnts )
    fop2.setPascalsStyle( True )
    xy = g.RPolynomialFunction( g.RVector( fop2.polynomialFunction().size() ) )
    xy.fill( fop2.startModel() )
    print "xy:", xy
    
    start = g.RVector( 27 )
    start.setVal( (xy).coeff()[0:9], 0, 9 )
    start.setVal( (xy).coeff()[0:9], 9, 18 )
    start.setVal( (xy).coeff()[0:9], 18, 27 )
    start[ 18 + 2 ] = 0
    start[ 18 + 4 ] = 0
    start[ 18 + 6 ] = 0
    print start
    fop.setStartModel( start ) 
    #print fop.setStartModel( )
    tmp = g.RPolynomialFunction( g.RVector( fop.polynomialFunction().size() ) )
    tmp.fill( fop.startModel() )
    print "base:", tmp
    

    inv = g.RInversion( uN, fop, False, False)
    inv.setRelativeError( 0.0 )
    inv.setLambda( 0 )
    inv.stopAtChi1( False )
    inv.setCGLSTolerance( 1e-40 )
    inv.setMaxIter( 20 )
    
    coeff = inv.run()
    
    
    N = fop.polynomialFunction()
    
    return N


def findMatrix( uN, dim, nCoeff, ptns, pascale, serendipity ):
    fop = g.PolynomialModelling( dim, nCoeff, pnts )
    fop.setPascalsStyle( pascale )
    fop.setSerendipityStyle( serendipity )
    
    tmp = g.RPolynomialFunction( g.RVector( fop.polynomialFunction().size() ) )
    tmp.fill( fop.startModel() )
    print "base:", tmp
    
    G = P.zeros( ( len( pnts), len( tmp.elements() ) ) )
    
    for i, p in enumerate( pnts ):
        for j, e in enumerate( tmp.elements() ):
            G[i,j] = e( p )
    
    GI = P.inv( G )
    coef = P.dot( GI, uN )
    
    print coef
    print P.dot( G, coef )
    tmp.fill( g.asvector( coef ) )
    print coef, tmp
    return tmp
     
def showSF( N, label = '' ):
    fig = P.figure()
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    ax2 = fig.add_subplot(1, 2, 2)
    
    Nx = 21
    Ny = 21
    nLevels = 12

    tix = P.linspace( 0.0, 1.0, Nx )
    tiy = P.linspace( 0.0, 1.0, Ny )
    (X,Y) = P.meshgrid( tix, tiy )
    z = g.RVector( len( X.flat ) )
    
    for i, x in enumerate( X.flat ):
        p = g.RVector3( X.flat[ i ], Y.flat[ i ] )
        z[ i ] = N(p)
    
    Z = P.ma.masked_where( z == -99., z )
    Z = Z.reshape( Ny, Nx )
    
    ax2.contourf( X, Y, Z )
    ax2.set_aspect( 'equal' )
    surf = ax1.plot_surface( X, Y, Z, rstride = 1, cstride = 1, cmap=P.cm.jet,linewidth=0 )
    
    ax2.set_title( label + N.__str__() )
    fig.colorbar( surf )
    
e1 = [ 1, -1 ] # 1-r 
e2 = [ 0,  1 ] #   r

E2_1R = g.RPolynomialFunction( g.asvector( e1 ) )
E2_2R = g.RPolynomialFunction( g.asvector( e2 ) )
E2_1S = g.RPolynomialFunction( g.RVector(), g.asvector( e1 ) )
E2_2S = g.RPolynomialFunction( g.RVector(), g.asvector( e2 ) )
E2_1T = g.RPolynomialFunction( g.RVector(), g.RVector(), g.asvector( e1 ) )
E2_2T = g.RPolynomialFunction( g.RVector(), g.RVector(), g.asvector( e2 ) )

E3_1R = E2_1R * ( E2_1R - E2_2R )
E3_2R = E2_2R * ( E2_2R - E2_1R )
E3_3R = 4. * ( E2_1R * E2_2R )

E3_1S = E2_1S * ( E2_1S - E2_2S )
E3_2S = E2_2S * ( E2_2S - E2_1S )
E3_3S = 4. * ( E2_1S * E2_2S )

E3_1T = E2_1T * ( E2_1T - E2_2T )
E3_2T = E2_2T * ( E2_2T - E2_1T )
E3_3T = 4. * ( E2_1T * E2_2T )

T3_1 = g.RPolynomialFunction( g.asvector( e1 ), g.asvector( e1 ) )
T3_2 = g.RPolynomialFunction( g.asvector( e2 ), g.RVector() )
T3_3 = g.RPolynomialFunction( g.RVector(), g.asvector( e2 ) )

T6_1 = T3_1 * ( 2.0 * T3_1 + -1.0 )
T6_2 = T3_2 * ( 2.0 * T3_2 + -1.0 )
T6_3 = T3_3 * ( 2.0 * T3_3 + -1.0 )
T6_4 = T3_1 * T3_2 * 4.0
T6_5 = T3_2 * T3_3 * 4.0
T6_6 = T3_3 * T3_1 * 4.0

#showSF( T )

Q4_1 = E2_1R * E2_1S
Q4_2 = E2_2R * E2_1S
Q4_3 = E2_2R * E2_2S
Q4_4 = E2_1R * E2_2S

Q8_5 = E3_3R * E2_1S
Q8_6 = E2_2R * E3_3S
Q8_7 = E3_3R * E2_2S
Q8_8 = E2_1R * E3_3S

Q8_1 = Q4_1 - Q8_8*0.5 - Q8_5*0.5
Q8_2 = Q4_2 - Q8_5*0.5 - Q8_6*0.5
Q8_3 = Q4_3 - Q8_6*0.5 - Q8_7*0.5
Q8_4 = Q4_4 - Q8_7*0.5 - Q8_8*0.5

T4_2 = g.RPolynomialFunction( g.asvector( e2 ), g.RVector() )
T4_3 = g.RPolynomialFunction( g.RVector(), g.asvector( e2 ) )
T4_4 = g.RPolynomialFunction( g.RVector(), g.RVector(), g.asvector( e2 ) )
T4_1 = -(-1.0 + T4_2 + T4_3 + T4_4)

T10_1 = T4_1 * ( 2.0 * T4_1 + -1.0 )
T10_2 = T4_2 * ( 2.0 * T4_2 + -1.0 )
T10_3 = T4_3 * ( 2.0 * T4_3 + -1.0 )
T10_4 = T4_4 * ( 2.0 * T4_4 + -1.0 )

T10_5  = 4.0 * T4_1 * T4_2;
T10_6  = 4.0 * T4_2 * T4_3;
T10_7  = 4.0 * T4_3 * T4_1;
T10_8  = 4.0 * T4_1 * T4_4;
T10_9  = 4.0 * T4_2 * T4_4;
T10_10 = 4.0 * T4_3 * T4_4;

H8_1 = E2_1R * E2_1S * E2_1T
H8_2 = E2_2R * E2_1S * E2_1T
H8_3 = E2_2R * E2_2S * E2_1T
H8_4 = E2_1R * E2_2S * E2_1T
H8_5 = E2_1R * E2_1S * E2_2T
H8_6 = E2_2R * E2_1S * E2_2T
H8_7 = E2_2R * E2_2S * E2_2T
H8_8 = E2_1R * E2_2S * E2_2T

P15_6 = T3_3 * E2_2T * ( 2. * ( T3_3 + E2_2T ) + -3. )
P15_7 = T6_4 * E2_1T
P15_8 = T6_5 * E2_1T
P15_9 = T6_6 * E2_1T
P15_10 = T6_4 * E2_2T
P15_11 = T6_5 * E2_2T
P15_12 = T6_6 * E2_2T
P15_13 = T3_1 * E3_3T
P15_14 = T3_2 * E3_3T
P15_15 = T3_3 * E3_3T

pnts = g.stdVectorRVector3()
    
# quad
pnts.append( g.RVector3( 0.0, 0.0 ) )
pnts.append( g.RVector3( 1.0, 0.0 ) )
pnts.append( g.RVector3( 1.0, 1.0 ) )
pnts.append( g.RVector3( 0.0, 1.0 ) )
pnts.append( g.RVector3( 0.5, 0.5 ) )
#pnts.append( g.RVector3( 0.5, 0.0 ) )
#pnts.append( g.RVector3( 1.0, 0.5 ) )
#pnts.append( g.RVector3( 0.5, 1.0 ) )
#pnts.append( g.RVector3( 0.0, 0.5 ) )

## tri
#pnts.append( g.RVector3( 0.0, 0.0 ) )
#pnts.append( g.RVector3( 1.0, 0.0 ) )
#pnts.append( g.RVector3( 0.0, 1.0 ) )
#pnts.append( g.RVector3( 0.5, 0.0 ) )
#pnts.append( g.RVector3( 0.5, 0.5 ) )
#pnts.append( g.RVector3( 0.0, 0.5 ) )

# tet
#pnts.append( g.RVector3( 0.0, 0.0, 0.0 ) )
#pnts.append( g.RVector3( 1.0, 0.0, 0.0 ) )
#pnts.append( g.RVector3( 0.0, 1.0, 0.0 ) )
#pnts.append( g.RVector3( 0.0, 0.0, 1.0 ) )
#pnts.append( ( pnts[ 0 ] + pnts[ 1 ] ) / 2.0 )
#pnts.append( ( pnts[ 1 ] + pnts[ 2 ] ) / 2.0 )
#pnts.append( ( pnts[ 2 ] + pnts[ 0 ] ) / 2.0 )
#pnts.append( ( pnts[ 0 ] + pnts[ 3 ] ) / 2.0 )
#pnts.append( ( pnts[ 1 ] + pnts[ 3 ] ) / 2.0 )
#pnts.append( ( pnts[ 2 ] + pnts[ 3 ] ) / 2.0 )

## hex
#pnts.append( g.RVector3( 0.0, 0.0, 0.0 ) )
#pnts.append( g.RVector3( 1.0, 0.0, 0.0 ) )
#pnts.append( g.RVector3( 1.0, 1.0, 0.0 ) )
#pnts.append( g.RVector3( 0.0, 1.0, 0.0 ) )
#pnts.append( g.RVector3( 0.0, 0.0, 1.0 ) )
#pnts.append( g.RVector3( 1.0, 0.0, 1.0 ) )
#pnts.append( g.RVector3( 1.0, 1.0, 1.0 ) )
#pnts.append( g.RVector3( 0.0, 1.0, 1.0 ) )
#pnts.append( (pnts[ 0 ] + pnts[ 1 ]) / 2.0 )
#pnts.append( (pnts[ 1 ] + pnts[ 2 ]) / 2.0 )
#pnts.append( (pnts[ 2 ] + pnts[ 3 ]) / 2.0 )
#pnts.append( (pnts[ 3 ] + pnts[ 0 ]) / 2.0 )
#pnts.append( (pnts[ 4 ] + pnts[ 5 ]) / 2.0 )
#pnts.append( (pnts[ 5 ] + pnts[ 6 ]) / 2.0 )
#pnts.append( (pnts[ 6 ] + pnts[ 7 ]) / 2.0 )
#pnts.append( (pnts[ 7 ] + pnts[ 4 ]) / 2.0 )
#pnts.append( (pnts[ 0 ] + pnts[ 4 ]) / 2.0 )
#pnts.append( (pnts[ 1 ] + pnts[ 5 ]) / 2.0 )
#pnts.append( (pnts[ 2 ] + pnts[ 6 ]) / 2.0 )
#pnts.append( (pnts[ 3 ] + pnts[ 7 ]) / 2.0 )

# pri
#pnts.append( g.RVector3( 0.0, 0.0, 0.0 ) )
#pnts.append( g.RVector3( 1.0, 0.0, 0.0 ) )
#pnts.append( g.RVector3( 0.0, 1.0, 0.0 ) )
#pnts.append( g.RVector3( 0.0, 0.0, 1.0 ) )
#pnts.append( g.RVector3( 1.0, 0.0, 1.0 ) )
#pnts.append( g.RVector3( 0.0, 1.0, 1.0 ) )

#pnts.append( (pnts[ 0 ] + pnts[ 1 ]) / 2.0 )
#pnts.append( (pnts[ 1 ] + pnts[ 2 ]) / 2.0 )
#pnts.append( (pnts[ 2 ] + pnts[ 0 ]) / 2.0 )
#pnts.append( (pnts[ 3 ] + pnts[ 4 ]) / 2.0 )
#pnts.append( (pnts[ 4 ] + pnts[ 5 ]) / 2.0 )
#pnts.append( (pnts[ 5 ] + pnts[ 3 ]) / 2.0 )
#pnts.append( (pnts[ 0 ] + pnts[ 3 ]) / 2.0 )
#pnts.append( (pnts[ 1 ] + pnts[ 4 ]) / 2.0 )
#pnts.append( (pnts[ 2 ] + pnts[ 5 ]) / 2.0 )

uN = g.RVector( len(pnts) )
uN[ 0 ] = 1
dim = 2
nCoeff = 3
pascale = False
serendipity = False
powCombination = 0

N = findShapeFunction( uN, dim, nCoeff, pnts, pascale, serendipity, powCombination )
#N = findMatrix( uN, dim, nCoeff, pnts, pascale, serendipity )

test = P15_6

#test = H8_1
print test

for p in pnts:
    print p, N(p), test(p)

print 'fopp', N
print 'test', test

print 'H8_1', H8_1
print 'H8_1', Q4_1 * E2_1T

print 'H8_2', H8_2
print 'H8_2', Q4_2 * E2_1T

print 'H8_5', H8_5
print 'H8_5', Q4_1 * E2_2T

print 'T10_8', T10_8
print 'T10_5', T10_5
print 'T10_1', T10_1
print 'T4_1', T4_1

print 'T6_2', T6_2
print 'T6_1', T6_1

print 'T3_1', T3_1

print 'Q4_1', Q4_1

print 'Q8_1', Q8_1


showSF( N, 'fop.f' )
showSF( test, 'test' )

#showSF( Q4_1, 'Q4_1' )
#showSF( Q8_6, 'Q8_6' )
#showSF( Q8_1, 'Q8_1' )
#showSF( Q8_3, 'Q8_3' )
#showSF( Q8_4, 'Q8_4' )



P.show()

