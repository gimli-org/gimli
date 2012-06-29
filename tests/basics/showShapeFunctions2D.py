#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from os import system
import pygimli as g

from pygimli.mplviewer import drawMeshBoundaries
import numpy as np
import pylab as P

from mpl_toolkits.mplot3d import axes3d

def drawShapes( axes, mesh, u ):
    axes.set_aspect( 'equal' )

    Nx = 20
    Ny = 20
    nLevels = 12

    tix = np.linspace( mesh.xmin(), mesh.xmax(), Nx )
    tiy = np.linspace( mesh.ymin(), mesh.ymax(), Ny )
    (X,Y) = np.meshgrid( tix, tiy )

    swatch = g.Stopwatch( True )
    print "interpolate prep t = ", swatch.duration( True )
    print u
    z = g.interpolate( mesh, u,   g.asvector( X.flat[:].tolist() )
                                , g.asvector( Y.flat[:].tolist() )
                                , g.RVector( len( Y.flat[:] ), 0.0 )
                                 )
    z = np.asarray ( z )
    print "interpolate t = ", swatch.duration( True )

    Z=z
    #Z = np.ma.masked_where( z <= 0.0, z )
    Z = Z.reshape( Ny, Nx )

    cs = axes.contourf( X, Y, Z, nLevels )
    drawMeshBoundaries( axes, mesh )

    tixG = np.linspace( 0, 1.0, 10 )
    tiyG = np.linspace( 0, 1.0, 10 )
    (X,Y) = np.meshgrid( tixG, tiyG )

    for i, x in enumerate( X.flat ):
        p = g.RVector3( X.flat[ i ], Y.flat[ i ])
        cell = mesh.findCell( p )
        if cell is None:
            continue
        gr = cell.grad( p,u )
        #print gr
        ux=gr[0]; uy=gr[1]
        #ux,uy = grad( mesh, p, u )
        axes.quiver( p[0], p[1], ux, uy )

    Zg = np.gradient( Z, 0.1, 0.1 );
    #axes.quiver( X, Y, Zg[1],Zg[0])

    axes.set_xlim( mesh.xmin(), mesh.xmax() );
    axes.set_ylim( mesh.ymin(), mesh.ymax() );


def show( mesh ):   
    fig = P.figure()
    for i in range( mesh.nodeCount() ):
        u = g.RVector( mesh.nodeCount() )    
        u[ i ] = 1.0;
        at = fig.add_subplot( 3, 3, i + 1)
        drawShapes( at, mesh, u )
      
def triangle( p2 = False ):
    mesh = g.Mesh()
    n1 = mesh.createNode( 0.0, 0.0, 0.0 )
    n2 = mesh.createNode( 1.0, 0.0, 0.0 )
    n3 = mesh.createNode( 1.0, 1.0, 0.0 )
    n4 = mesh.createNode( 0.0, 1.0, 0.0 )
    mesh.createTriangle( n1, n2, n4 )
    mesh.createTriangle( n2, n3, n4 )
    
    if p2:
        tmp = g.Mesh( mesh )
        mesh = g.Mesh(); mesh.createP2Mesh( tmp )
        
    show( mesh );    
           
    
def quadrangle( p2 = False ):
    mesh = g.Mesh()
    n1 = mesh.createNode( 0.0, 0.0, 0.0 )
    n2 = mesh.createNode( 1.0, 0.0, 0.0 )
    n3 = mesh.createNode( 1.0, 1.0, 0.0 )
    n4 = mesh.createNode( 0.0, 1.0, 0.0 )
    mesh.createQuadrangle( n1, n2, n3, n4 )

    if p2:
        tmp = g.Mesh( mesh )
        mesh = g.Mesh(); mesh.createP2Mesh( tmp )
        
    show( mesh );    
   

triangle( p2 = False )
triangle( p2 = True )
quadrangle( p2 = False )
quadrangle( p2 = True )

P.show()
