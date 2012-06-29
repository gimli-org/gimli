# -*- coding: iso-8859-1 -*-
from os import system
import pygimli as g

from pygimli.mplviewer import drawMeshBoundaries
import numpy as np
import pylab as P

def drawShapes( axes, edge, u ):
    axes.set_aspect( 'equal' )

    Nx = 20

    tix = np.linspace( edge.node(0).pos()[0], edge.node(1).pos()[0], Nx )

    z = []
    for x in tix:
        z.append( edge.interpolate( g.RVector3( x,0.0 ), u ) )

    axes.plot( tix, z )

    tixG = np.linspace( 0, 1.0, 10 )

    for i, x in enumerate( tixG ):
        p = g.RVector3( tixG[ i ], 0.0 )
        gr = edge.grad( p,u )
        ##print gr
        ux=gr[0]; uy=gr[1]
        ##ux,uy = grad( mesh, p, u )
        axes.quiver( p[0], 0.5, ux, uy )

    axes.set_aspect('equal')
    #axes.set_xlim( mesh.xmin(), mesh.xmax() );

fig = P.figure()

mesh1 = g.Mesh( )
n1 = mesh1.createNode( 0.0, 0.0, 0.0 )
n2 = mesh1.createNode( 1.0, 0.0, 0.0 )
edge1 = g.EdgeCell( n1, n2 )
u1 = g.RVector( mesh1.nodeCount() )

mesh2 = g.Mesh( )
n1 = mesh2.createNode( 0.0, 0.0, 0.0 )
n2 = mesh2.createNode( 1.0, 0.0, 0.0 )
n3 = mesh2.createNode( 0.5, 0.0, 0.0 )
nodes = g.stdVectorNodes()
nodes.append(n1); nodes.append(n2); nodes.append(n3)
edge2 = g.Edge3Cell( nodes )

for i in range( 2 ):
    u = g.RVector( mesh1.nodeCount() )
    u[ i ] = 1.0;
    at = fig.add_subplot(2,3,1+i)
    drawShapes( at, edge1, u )

for i in range( 3 ):
    u = g.RVector( mesh2.nodeCount() )
    u[ i ] = 1.0;
    at = fig.add_subplot(2,3,4+i)
    drawShapes( at, edge2, u )

P.show()