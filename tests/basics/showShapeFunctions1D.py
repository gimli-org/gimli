# -*- coding: iso-8859-1 -*-
from os import system
import pygimli as g

from pygimli.mplviewer import drawMeshBoundaries
import numpy as np
import pylab as P

def drawShapes( ax, edge, u ):
    Nx = 41
    X = np.linspace( -2, 2, Nx )
    
    uc = g.RVector( len( X.flat ) )
    
    for i, x in enumerate( X ):
        p = edge.shape().xyz( g.RVector3( X.flat[ i ], 0.0) )
        
        X.flat[ i ] = p[ 0 ]
        
        if not edge.shape().isInside( p ): 
            uc[ i ] = -99.0
            continue

        gr = edge.grad( p, u )
        uc[ i ] = edge.pot( p, u )
        ax.quiver( p[0], uc[i], gr[0], 0 )
        
        
    Z = np.ma.masked_where( uc == -99., uc )
    ax.plot( X, Z )
    ax.set_aspect('equal')
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