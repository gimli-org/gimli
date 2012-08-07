#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from os import system
import pygimli as g

from pygimli.mplviewer import drawMeshBoundaries
import numpy as np
import pylab as P

from mpl_toolkits.mplot3d import axes3d

def drawShapes( ax, mesh, u ):
    '''
    '''
    
    ax.set_aspect( 'equal' )

    Nx = 21
    Ny = 21
    nLevels = 12

    tix = np.linspace( -1.0, 1.0, Nx )
    tiy = np.linspace( -1.0, 1.0, Ny )
    (X,Y) = np.meshgrid( tix, tiy )

    uc = g.RVector( len( X.flat ) )
    pnts = []
    
    c = mesh.cell( 0 )
    
    imax = g.find( u == max( u ) )[0]
    
    for i in range( c.nodeCount() ):
        print c.rst( i )
    print imax
    print c.createShapeFunctions()[ imax ]
    print "dx", c.createShapeFunctions()[ imax ].derive( 0 )
    print "dy", c.createShapeFunctions()[ imax ].derive( 1 )
    
    # draw nodes
    for i in range( c.nodeCount() ):
        col = 'black'
        if i == imax:
            col = 'red'
            
        ax.plot( c.node( i ).pos()[0], c.node( i ).pos()[1], '.', markersize = 12, linewidth=0, color=col )
        
    # draw boundary
    drawMeshBoundaries( ax, mesh )
    ptns = []
    grads = []

    swatch = g.Stopwatch( True )
    for i, x in enumerate( X.flat ):
        p = c.shape().xyz( g.RVector3( X.flat[ i ], Y.flat[ i ] ) )
        
        X.flat[ i ] = p[ 0 ]
        Y.flat[ i ] = p[ 1 ]
        
        #ax.plot( p[0], p[1], '.', zorder=10, color='black', markersize = 1 )
        
        if not c.shape().isInside( p ): 
            uc[ i ] = -99.0
            continue
            
        uc[ i ] = c.pot( p, u ) 
        
        gr = c.grad( p, u )
        ptns.append( p )
        grads.append( gr )

    print swatch.duration( True )
    
    for i, p in enumerate( ptns ):
        #print p, grads[i]
        ax.quiver( p[0], p[1], grads[i][0], grads[i][1], zorder=10)

    Z = np.ma.masked_where( uc == -99., uc )
    Z = Z.reshape( Ny, Nx )
    cs = ax.contourf( X, Y, Z, nLevels )
    
# def drawShapes( ... )

def show( mesh ):   
    fig = P.figure()
    
    for i in range( mesh.nodeCount() ):
        u = g.RVector( mesh.nodeCount() )    
        u[ i ] = 1.0;
        ax = fig.add_subplot( 3, 3, i + 1)
        drawShapes( ax, mesh, u )
      
def triangle( p2 = False ):
    mesh = g.Mesh()
    n1 = mesh.createNode( 0.0, 0.0, 0.0 )
    n2 = mesh.createNode( 1.0, 0.0, 0.0 )
    #n3 = mesh.createNode( 1.0, 1.0, 0.0 )
    n4 = mesh.createNode( 0.0, 1.0, 0.0 )
    mesh.createTriangle( n1, n2, n4 )
    mesh.rotate( g.RVector3( 0.0, 0.0, 45 ) )
    #mesh.createTriangle( n2, n3, n4 )
    
    if p2:
        mesh = mesh.createP2Mesh( )
        
    show( mesh )    
    
def quadrangle( p2 = False ):
    mesh = g.Mesh()
    n1 = mesh.createNode( 0.0, 0.0, 0.0 )
    n2 = mesh.createNode( 1.0, 0.0, 0.0 )
    n3 = mesh.createNode( 1.0, 1.0, 0.0 )
    n4 = mesh.createNode( 0.0, 1.0, 0.0 )
    mesh.createQuadrangle( n1, n2, n3, n4 )
    mesh.rotate( g.RVector3( 0.0, 0.0, 45 ) )
    mesh.scale( g.RVector3( 1.0, 0.5, 1.0 ) )
    
    if p2:
        mesh = mesh.createP2Mesh( )
        
    show( mesh )    
   
triangle( p2 = False )
triangle( p2 = True )
quadrangle( p2 = False )
quadrangle( p2 = True )

P.show()
