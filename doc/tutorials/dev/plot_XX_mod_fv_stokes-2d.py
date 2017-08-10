#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
import sys

import pygimli as pg
import pygimli.solver as solver

from pygimli.meshtools import polyCreateDefaultEdges_
from pygimli.meshtools import createMesh

import matplotlib.pyplot as plt
import numpy as np

#from solverFVM import solveStokes_NEEDNAME

def modelCavity0(maxArea=0.0025):
    mesh = pg.createGrid(x=np.linspace(.0, 1.0, 17),
                         y=np.linspace(.0, 1.0, 17))
    mesh = pg.meshtools.refineQuad2Tri(mesh, style=2)
    
    velBoundary=[[1,[0.0, 0.0]],
                 [2,[0.0, 0.0]],
                 [3,[1.0, 0.0]],
                 [4,[0.0, 0.0]],
                ]

    c = mesh.findCell((0.0, 0.0))
    for b in range(c.boundaryCount()):
        if c.boundary(b).marker()==4:
            c.boundary(b).setMarker(7)
    
    preBoundary=[[7, 0.0],]
    
    a = pg.RVector(mesh.cellCount(), 1.0)
    return mesh, velBoundary, preBoundary, a, 100

def modelCavity1(maxArea=0.0025):
    boundary = []
    boundary.append([-1.0, -1.0])
    boundary.append([ -0.5, -1.0])
    boundary.append([ -0.5, -0.7])
    boundary.append([ 0.5, -0.7])
    boundary.append([ 0.5, -1.0])
    boundary.append([ 1.0, -1.0])
    boundary.append([ 1.0, 1.0])
    boundary.append([-1.0, 1.0])
    
    poly = pg.Mesh(2)
    nodes = [poly.createNode(b) for b in boundary]

    polyCreateDefaultEdges_(poly, boundaryMarker=[4,4,4,4,4,2,3,1])
    mesh = createMesh(poly, quality=33.4, area=maxArea, smooth=[0,10])

    # Diffusions coefficient, viscosity
 
    b7 = mesh.findBoundaryByMarker(1)[0]
    for b in mesh.findBoundaryByMarker(1):
        if b.center()[1] < b.center()[1]:
            b7 = b
    b7.setMarker(7)

    velBoundary=[ [1,[0.0, 0.0]],
                [2,[0.0, 0.0]],
                [3,[1.0, 0.0]],
                [4,[0.0, 0.0]],
                [7,[0.0, 0.0]]]

    preBoundary=[[7,0.0]]
    
    a = pg.RVector(mesh.cellCount(), 10000.0)
    return mesh, velBoundary, preBoundary, a, 100

def modelCavity2(area, refine=True):
    boundary = []
    boundary.append([-1.0,  0.0]) #0
    boundary.append([-1.0, -1.0]) #1
    
    boundary.append([-0.2, -1.0]) #2
    boundary.append([-0.2, -0.8]) #3
    boundary.append([ 0.2, -0.8]) #4
    boundary.append([ 0.2, -1.0]) #5
    
    boundary.append([ 1.0, -1.0]) #6
    boundary.append([ 1.0,  0.0]) #7
    
    boundary.append([ 0.2, 0.0]) #8
    boundary.append([ 0.2, -0.2]) #9
    boundary.append([-0.2, -0.2]) #10
    boundary.append([-0.2, 0.0]) #11
    
    poly = pg.Mesh(2)
    nodes = [poly.createNode(b) for b in boundary]

    eMarker = [1, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2]
    for i in range(len(nodes)): 
        poly.createEdge(nodes[i], nodes[(i+1)%len(nodes)], eMarker[i])
    
    if refine:
        poly.createNode([-0.21, -0.99])
    mesh = createMesh(poly, quality=34.0, area=area, smooth=[0,10])

    # Diffusions coefficient, viscosity
    
    b7 = mesh.findBoundaryByMarker(2)[0]
    for b in mesh.findBoundaryByMarker(1):
        if b.center()[1] < b.center()[1]:
            b7 = b
    b7.setMarker(4)

    velBoundary=[[1,[1.0, 0.0]],
                [2,[0.0, 0.0]],
                [3,[1.0, 0.0]],
                [4,[0.0, 0.0]],
                ]
    preBoundary=[[4,0.0]]
    
    a = pg.RVector(mesh.cellCount(), 1.0)
    return mesh, velBoundary, preBoundary, a, 50000

def modelPipe():
    boundary = []
    #left inflow
    boundary.append([0.1, 0.0]) # 0 
    boundary.append([0.0, 0.0]) # 1
    
    boundary.append([0.0, -1.0])
    boundary.append([0.0, -1.1])
    
    boundary.append([1.0, -1.1])
    boundary.append([1.0, -1.0])
    
    #right outflow
    boundary.append([1.0, 0.0]) # 1
    boundary.append([0.9, 0.0]) # 0
    
    #closing
    boundary.append([0.9, -1.0]) # 0
    boundary.append([0.1, -1.0]) # 0 
    
    poly = pg.Mesh(2)
    nodes = [poly.createNode(b) for b in boundary]
    
    for i in range(len(nodes)): poly.createEdge(nodes[i], nodes[(i+1)%len(nodes)], 1) 
    poly.boundaries()[0].setMarker(2)
        
    for b in poly.boundaries():
        if b.norm()[1]==1.0 and b.center()[0]==0.95:
            b.setMarker(3)
        if b.norm()[1]==-1.0 and b.center()[0]==0.5:
            b.setMarker(4)
    
    mesh = createMesh(poly, quality=34.0, area=0.0005, smooth=[0,10])
    
    velBoundary=[[1, [0.0,  0.0]],
                 [2, [0.0, -1.0]],
                 [3, [0.0,  1.0]],
                 [4, [0.0,  0.0]]
                 ]

    #preBoundary=None
    preBoundary=[[4, 0.0]]
    
    a=1
    return mesh, velBoundary, preBoundary, a, 400

def modelPlume(maxArea=0.1):
    boundary = []
    
    boundary.append([-1.,   0.0])#0
    boundary.append([-1., -1.0])#1
    
    boundary.append([-0.1,  -1.0])#2
    boundary.append([ 0.1,  -1.0])#3
    
    boundary.append([ 1., -1.0])#4
    boundary.append([ 1.,    0.0])#5
    boundary.append([ 0.1,    0.0])#6
    boundary.append([ -0.1,   0.0])#7
    
    poly = pg.Mesh(2)
    nodes = [poly.createNode(b) for b in boundary]
    
    polyCreateDefaultEdges_(poly, boundaryMarker=[1,2,3,2,4,5,6,7])
    
    #poly.createEdge(nodes[0], nodes[1], 1) # left
    #poly.createEdge(nodes[1], nodes[2], 2) # bottom1
    #poly.createEdge(nodes[2], nodes[3], 3) # bottom2
    #poly.createEdge(nodes[3], nodes[4], 2) # bottom3
    #poly.createEdge(nodes[4], nodes[5], 4) # right
    #poly.createEdge(nodes[5], nodes[6], 5) # top1
    #poly.createEdge(nodes[6], nodes[7], 6) # topcenter
    #poly.createEdge(nodes[7], nodes[0], 7) # top2
    
    mesh = createMesh(poly, quality=33.4, area=maxArea,
                      smooth=[0,10], verbose=False)
    
    velBoundary=[#[1, [0.0,  0.0]],
 #                [2, [0.0,  0.0]],
                 #[3, [0.0,  0.0]],
                 #[4, [0.0,  0.0]],
                 [5, [1.0,  0.0]],
                 [6, [0.0,  0.0]],
                 [7, [-1.0,  0.0]]
                 ]

    b = mesh.findBoundaryByMarker(2)[0]
    b.setMarker(8)
    #preBoundary=None
    preBoundary=[[8, 0.0]]
    
    a=1
    return mesh, velBoundary, preBoundary, a, 100



modelBuilder = modelCavity0
#modelBuilder = modelCavity1
#modelBuilder = modelCavity2
#modelBuilder = modelPipe
#modelBuilder = modelPlume


swatchG = pg.Stopwatch(True)
swatch = pg.Stopwatch(True)

nSteps = 1
multigridArea = (10.**(np.linspace(np.log10(0.001), np.log10(0.1), nSteps)))[::-1]
print(multigridArea)

pre = None
vel = None
    
for i in range(0, len(multigridArea)):
    if i == 0:
        mesh, velBoundary, preBoundary, a, maxIter = modelBuilder(multigridArea[0])
        pre = np.zeros(mesh.cellCount())
        vel = np.zeros((mesh.cellCount(), 3))
    else:
        #mesh1 = mesh.createH2()
        mesh1, velBoundary, preBoundary, a, maxIter = modelBuilder(multigridArea[i])
        a = pg.RVector(mesh1.cellCount(), 1.0)
         
        pre = pg.interpolate(mesh, pre, mesh1.cellCenter())
        vx0 = pg.interpolate(mesh, vel[:,0], mesh1.cellCenter())
        vy0 = pg.interpolate(mesh, vel[:,1], mesh1.cellCenter())
        vel = np.vstack([vx0, vy0]).T
        mesh = mesh1
    
    print("Cells: ", mesh.cellCount(), multigridArea[i])
    vel, pre, pCNorm, divVNorm = pg.solver.solveStokes(mesh, 
                                                       viscosity=a,
                                                       velBoundary=velBoundary, 
                                                       preBoundary=preBoundary,
                                                       pre0=pre,
                                                       vel0=vel,
                                                       maxIter=maxIter,
                                                       tol=1e-2,
                                                       verbose=1)
    print(" Time: ", swatch.duration(True))
    
   

print("OverallTime:", swatchG.duration(True))

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)
ax, cbar = pg.show(mesh, 
                   data=pg.cellDataToPointData(mesh, pre),
                   logScale=False, colorBar=True, axes=ax1)
cbar.ax.set_xlabel('Pressure in ??')
meshC, velBoundary, preBoundary, a, maxIter = modelBuilder(0.01)
pg.show(mesh, data=vel, coarseMesh=meshC, axes=ax1)
pg.show(mesh, axes=ax1)

ax2 = fig.add_subplot(1, 2, 2)
ax, cbar = pg.show(mesh, 
                   data=pg.cellDataToPointData(mesh,
                                    np.sqrt(vel[:,0]*vel[:,0] +vel[:,1]*vel[:,1])),
                   logScale=False, colorBar=True, axes=ax2)
cbar.ax.set_xlabel('Geschwindigkeit in m$/$s')
pg.show(mesh, data=vel, coarseMesh=meshC, axes=ax2)
pg.show(mesh, axes=ax2)

        
plt.figure()
plt.semilogy(pCNorm, label='norm')
plt.semilogy(divVNorm, label='norm div v')
plt.legend()

pg.wait()
#drawMesh(ax, mesh)