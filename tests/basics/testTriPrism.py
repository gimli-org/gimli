#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as g
import numpy as np

poly = g.Mesh( 2 )
n0 = poly.createNode( 0.0, 0.0, 0. )
n1 = poly.createNode( 1.0, 0.0, 0. )
n2 = poly.createNode( 0.0, 1.0, 0. )
n3 = poly.createNode( 1.0, 1.0, 0. )
poly.createEdge( n0, n1 )
poly.createEdge( n1, n3 )
poly.createEdge( n3, n2 )
poly.createEdge( n2, n0 )
    
mesh2= g.Mesh( 2 )
g.TriangleWrapper( poly, mesh2, "-pzeAfa0.01q34" );

mesh3 = g.createMesh3D( mesh2, g.asvector( np.arange(0, -1, -0.1 ) ) )
mesh3.setCellAttributes( g.asvector( range(0, mesh3.cellCount() ) ) )

print mesh3

mesh3.save( "prism" )
mesh3.exportVTK( "prism" )