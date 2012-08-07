#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as g
from pygimli.viewer import showMesh

import pylab as P

import numpy as np

def testCreateTriPrismMesh():
    poly = g.Mesh( 2 )
    n0 = poly.createNode( 0.0, 0.0, 0. )
    n1 = poly.createNode( 1.0, 0.0, 0. )
    n2 = poly.createNode( 0.0, 1.0, 0. )
    n3 = poly.createNode( 1.0, 1.0, 0. )
    poly.createEdge( n0, n1 )
    poly.createEdge( n1, n3 )
    poly.createEdge( n3, n2 )
    poly.createEdge( n2, n0 )
    
    mesh2 = g.Mesh( 2 )
    g.TriangleWrapper( poly, mesh2, "-pzeAfa0.01q34" );

    mesh3 = g.createMesh3D( mesh2, g.asvector( np.arange(0, -1, -0.1 ) ) )
    mesh3.setCellAttributes( g.asvector( range(0, mesh3.cellCount() ) ) )

    mesh3.save( "prism" )
    mesh3.exportVTK( "prism" )
# def testCreateTriPrismMesh( ... )
    
def testShapefunctions( ):
    poly = g.Mesh( 3 )
    n0 = poly.createNode( 0.0, 0.0, 0.0 )
    n1 = poly.createNode( 1.0, 0.0, 0.0 )
    n2 = poly.createNode( 0.0, 1.0, 0.0 )
    n3 = poly.createNode( 0.0, 0.0, 0.5 )

    poly.createTetrahedron( n0, n1, n2, n3 )
    prism = poly.createP2Mesh()
            
    
    u = g.RVector( prism.nodeCount(), 0.0 ) 
    u[ 5 ] = 1.0
    prism.addExportData( "u", u )
    
    #n3 = poly.createNode( 1.0, 1.0, 0.0 )
    #poly.createTriangle( n0, n1, n3 )
    #prism = g.createMesh3D( poly, g.asvector( np.linspace( 0, -1, 2 ) ) )
    #prism.exportVTK( "prism" )

    
    
    #mesh2 = g.createMesh2D( g.asvector( np.linspace( 0, 1, 10 ) ), 
                            #g.asvector( np.linspace( 0, 1, 10 ) ) )
                                
    #mesh2 = mesh2.createH2Mesh()
                            
    #g.interpolate( prism, mesh2 )
        
    #ax = g.viewer.showMesh( mesh2, mesh2.exportData('u'), filled = True, showLater = True )
        
    mesh3 = g.createMesh3D( g.asvector( np.linspace( 0, 1, 11 ) ), 
                            g.asvector( np.linspace( 0, 1, 11 ) ), 
                            g.asvector( np.linspace( 0, 1, 11 ) ) )

    grads = g.stdVectorRVector3( )
    c = prism.cell( 0 )
    uc = g.RVector( mesh3.nodeCount() )
    
    for n in mesh3.nodes():
        p = c.shape().xyz( n.pos() )
        
        if not c.shape().isInside( p ):
            grads.append( g.RVector3( 0.0, 0.0, 0.0 ) )
            uc[ n.id() ] = 0.0
            continue
        
        uc[ n.id() ] = c.pot( p, u )
        print uc[ n.id() ]
        gr = c.grad( p, u )
        grads.append( gr )
        
    g.interpolate( prism, mesh3 )
    mesh3.addExportData( 'ua', uc )

    mesh3.exportVTK( "prismHex", grads )
    
    
    P.show()
    
    
testShapefunctions( )