#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
import sys

import pygimli as pg
import pygimli.solver as solver
from pygimli.viewer import show

from pygimli.meshtools import createMesh

import matplotlib.pyplot as plt
import numpy as np

from solverFVM import solveStokes_NEEDNAME

from solverFVM import solveFiniteVolume, createFVPostProzessMesh, diffusionConvectionKernel, cellToFace
from solverFVM import boundaryToCellDistances
from solverFVM import cellDataToCellGrad, cellDataToCellGrad2, divergence
from solverFVM import cellDataToBoundaryDataMatrix, cellDataToBoundaryGrad


def modelCavity():
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

    poly.createEdge(nodes[0], nodes[1], 4) # bottom
    poly.createEdge(nodes[1], nodes[2], 4) # bottom
    poly.createEdge(nodes[2], nodes[3], 4) # bottom
    poly.createEdge(nodes[3], nodes[4], 4) # bottom
    poly.createEdge(nodes[4], nodes[5], 4) # bottom

    poly.createEdge(nodes[5], nodes[6], 2) # right
    poly.createEdge(nodes[6], nodes[7], 3) # top
    poly.createEdge(nodes[7], nodes[0], 1) # left

    mesh = createMesh(poly, quality=33.4, area=0.0025, smooth=[0,10])

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

def modelPlume():
    boundary = []
    
    boundary.append([-100.,    0.0])#0
    boundary.append([-100., -100.0])#1
    
    boundary.append([-10.,  -100.0])#2
    boundary.append([ 10.,  -100.0])#3
    
    boundary.append([ 100., -100.0])#4
    boundary.append([ 100.,    0.0])#5
    boundary.append([ 10. ,    0.0])#6
    boundary.append([ -10. ,   0.0])#7
    
    
    poly = pg.Mesh(2)
    nodes = [poly.createNode(b) for b in boundary]
    
    poly.createEdge(nodes[0], nodes[1], 1) # left
    poly.createEdge(nodes[1], nodes[2], 2) # bottom1
    poly.createEdge(nodes[2], nodes[3], 3) # bottom2
    poly.createEdge(nodes[3], nodes[4], 2) # bottom3
    poly.createEdge(nodes[4], nodes[5], 4) # right
    poly.createEdge(nodes[5], nodes[6], 5) # top1
    poly.createEdge(nodes[6], nodes[7], 6) # topcenter
    poly.createEdge(nodes[7], nodes[0], 7) # top2
    
    mesh = createMesh(poly, quality=33.4, area=20., smooth=[0,10], verbose=False)
    
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

#mesh, velBoundary, preBoundary, a, maxIter = modelCavity()
mesh, velBoundary, preBoundary, a, maxIter = modelCavity2(0.10001101937239093957)
#mesh, velBoundary, preBoundary, a, maxIter= modelPipe()
#mesh, velBoundary, preBoundary, a, maxIter= modelPlume()


modelBuilder = modelCavity2

swatchG = pg.Stopwatch(True)
swatch = pg.Stopwatch(True)

nSteps = 3
aRange = (10.**(np.linspace(np.log10(0.0005), np.log10(0.1), nSteps)))[::-1]
print( aRange)

pre = np.zeros(mesh.cellCount())
vel = np.zeros((mesh.cellCount(), 3))
    
for i in range(0, len(aRange)):
    if i == 0:
        mesh, velBoundary, preBoundary, a, maxIter = modelCavity2(aRange[0])
    else:
        #mesh1 = mesh.createH2()
        mesh1, velBoundary, preBoundary, a, maxIter = modelBuilder(aRange[i])
        a = pg.RVector(mesh1.cellCount(), 1.0)
         
        pre = pg.interpolate(mesh, pre, mesh1.cellCenter())
        vx0 = pg.interpolate(mesh, vel[:,0], mesh1.cellCenter())
        vy0 = pg.interpolate(mesh, vel[:,1], mesh1.cellCenter())
        vel = np.vstack([vx0, vy0]).T
        mesh = mesh1
    
    print("Cells: ", mesh.cellCount(), aRange[i])
    vel, pre, pCNorm, divVNorm = solveStokes_NEEDNAME(mesh, velBoundary, preBoundary,
                                            viscosity=a,
                                            pre0=pre,
                                            vel0=vel,
                                            maxIter=maxIter,
                                            tol=1e-2,
                                            verbose=1)
    print(" Time: ", swatch.duration(True))
    
   

print("OverallTime:", swatchG.duration(True))
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
plt.ion()
ax, cbar = show(mesh, data=pg.cellDataToPointData(mesh, np.sqrt(vel[:,0]*vel[:,0] +vel[:,1]*vel[:,1])),
     logScale=False, colorBar=True, axes=ax1, cbar='b2r')
cbar.ax.set_xlabel('Geschwindigkeit in m$/$s')
     
meshC, velBoundary, preBoundary, a, maxIter = modelCavity2(0.001, False)

show(mesh, data=vel, coarseMesh=meshC, axes=ax1)
#show(mesh, data=vel, axes=ax1)

show(meshC, axes=ax1)
#show(mesh, axes=ax1)
        
plt.figure()
plt.semilogy(pCNorm, label='norm')
plt.semilogy(divVNorm, label='norm')
plt.legend()

plt.ioff()
plt.show()
#drawMesh(ax, mesh)