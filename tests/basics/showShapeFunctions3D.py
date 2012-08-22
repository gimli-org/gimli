#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import pygimli as g
from pygimli.viewer import showMesh

import pylab as P

import numpy as np

from mpl_toolkits.mplot3d import Axes3D

def drawShapes( ax, mesh, u ):
    #ax.set_aspect( 'equal' )
    N = 11
    
    mesh3 = g.createMesh3D( g.asvector( np.linspace( 0, 1, N ) ), 
                            g.asvector( np.linspace( 0, 1, N ) ), 
                            g.asvector( np.linspace( 0, 1, N ) ) )
                            
    uc = g.RVector( mesh3.nodeCount(  ) )
    
    grads = g.stdVectorRVector3( )
    pnts = g.stdVectorRVector3( )
    
    c = mesh.cell( 0 )
    imax=g.find( u == max( u ) )[0]
    
    N = c.createShapeFunctions()[ imax ] 
    print imax, N
#    print imax, c.shape().createShapeFunctions()[imax]
    for i in range( c.nodeCount() ):
        print c.rst( i ), N( c.rst( i ) )
    
    # draw nodes
    for i in range( c.nodeCount() ):
        col = 'black'
        if i == imax:
            col = 'red'
            
        #ax.plot( [c.node( i ).pos()[0], c.node( i ).pos()[0] ], 
                 #[c.node( i ).pos()[1], c.node( i ).pos()[1] ],
                 #[c.node( i ).pos()[2], c.node( i ).pos()[2] ],
                #'.', markersize = 15, linewidth=0, color=col )
    
    newNode = []
    
    for i in range( mesh3.nodeCount(  ) ):
        p = c.shape().xyz( mesh3.node( i ).pos() )
        newNode.append( p )
        
        #ax.plot( p[0], p[1], '.', zorder=10, color='black', markersize = 1 )
        
        if not c.shape().isInside( p ): 
            uc[ i ] = -99.0
            grads.append( g.RVector3( 0.0, 0.0 ) )
            continue
          
        uc[ i ] = c.pot( p, u ) 
        gr = c.grad( p, u ).normalise()
        grads.append( gr )
        
        #ax.plot( [ p[ 0 ], p[ 0 ] + gr[ 0 ]*0.1 ], 
                 #[ p[ 1 ], p[ 1 ] + gr[ 1 ]*0.1 ], 
                 #[ p[ 2 ], p[ 2 ] + gr[ 2 ]*0.1 ], '-', color='black' )
    
        #pnts.append( p )
            
    #print len(pnts)
    #Z = np.ma.masked_where( uc == -99., uc )
    ##ax.plot( g.x(pnts), g.y(pnts), g.z(pnts), '.' )
    
    for i, n in enumerate( mesh3.nodes() ):
        n.setPos( newNode[ i ] )
        
    mesh3.addExportData( 'u', uc.setVal( 0.0, g.find( uc == -99 ) ) )
    name = 'cell' + str( c.nodeCount() ) + '-' + str( imax )
    print "write ", name
    mesh3.exportVTK( name, grads )
    
def show( mesh ):   
    fig = P.figure()
    
    for i in range( mesh.nodeCount() ):
        u = g.RVector( mesh.nodeCount() )    
        u[ i ] = 1.0;
        
        ax = fig.gca()
#        ax = fig.add_subplot( 4, 3, i + 1, projection='3d')
        
        drawShapes( ax, mesh, u )    
    
def tetrahedron( p2 ):
    
    mesh = g.Mesh( 3 )
    n0 = mesh.createNode( 0.0, 0.0, 0.0 )
    n1 = mesh.createNode( 1.0, 0.0, 0.0 )
    n2 = mesh.createNode( 0.0, 1.0, 0.0 )
    n3 = mesh.createNode( 0.0, 0.0, 1.0 )

    mesh.createTetrahedron( n0, n1, n2, n3 )
    mesh.rotate( g.RVector3( 10.0, 34.0, 45 ) )
    mesh.scale( g.RVector3( 0.9, 0.8, 0.7 ) )

    if p2: mesh = mesh.createP2()
    
    show( mesh )

def hexahedron( p2 ):
    """
    """
    mesh = g.createMesh3D( g.asvector( np.linspace( 0, 1, 2 ) ), 
                           g.asvector( np.linspace( 0, 1, 2 ) ), 
                           g.asvector( np.linspace( 0, 1, 2 ) ) )

    mesh.node(0).setPos( g.RVector3( 0.0, 0.0, 0.5 ) )
    #mesh.rotate( g.RVector3( 10.0, 34.0, 45 ) )
    #mesh.scale( g.RVector3( 0.9, 0.8, 0.7 ) )

    #mesh.rotate( g.RVector3( 10.0, 34.0, 45 ) )
    #mesh.scale( g.RVector3( 1.0, 1.0, 5.0 ) )

    if p2: mesh = mesh.createP2()

    show( mesh )

def prism( p2 ):
    
    mesh = g.Mesh( 3 )
    n0 = mesh.createNode( 0.0, 0.0, 0.0 )
    n1 = mesh.createNode( 1.0, 0.0, 0.0 )
    n2 = mesh.createNode( 0.0, 1.0, 0.0 )
    n4 = mesh.createNode( 0.0, 0.0, 1.0 )
    n5 = mesh.createNode( 1.0, 0.0, 1.0 )
    n6 = mesh.createNode( 0.0, 1.0, 1.4 )
        
    ns = g.stdVectorNodes()
    for n in mesh.nodes():
        ns.append( n )
    
    mesh.createCell( ns, 0 )
    mesh.rotate( g.RVector3( 10.0, 34.0, 45 ) )
    mesh.scale( g.RVector3( 0.9, 0.8, 0.7 ) )
    
    if p2: 
        mesh = mesh.createP2()
        mesh.exportVTK( "prism15" )
    else:
        mesh.exportVTK( "prism6" )
   
    show( mesh )

    
#tetrahedron( p2 = False )
#tetrahedron( p2 = True )
#hexahedron( p2 = False )
hexahedron( p2 = True )
#prism( p2 = False )
#prism( p2 = True )


#mesh = g.createMesh3D( g.asvector( np.linspace( 0, 1, 5 ) ), 
                        #g.asvector( np.linspace( 0, 1, 5 ) ), 
                           #g.asvector( np.linspace( 0, 1, 5 ) ) )
                           

##mesh = g.Mesh(2)
##n0 = mesh.createNode( 0.0, 0.0, 0.0 )
##n1 = mesh.createNode( 1.0, 0.0, 0.0 )
##n2 = mesh.createNode( 0.0, 1.0, 0.0 )


##mesh.createTriangle( n0, n1, n2)
#mesh.createNeighbourInfos()

#print mesh
##mesh = g.createMesh2D( g.asvector( np.linspace( 0, 1, 2 ) ), 
                       ##g.asvector( np.linspace( 0, 1, 2 ) ) )

#c = mesh.findCell( g.RVector3( 1.0, 1.0, 0.8 ), False )

#print c
    