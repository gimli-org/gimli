#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import pygimli as pg
import pygimli.solver as solver
from pygimli.viewer import show, showMesh
from pygimli.mplviewer import drawMesh, drawModel, drawField, drawStreams
from pygimli.meshtools import createMesh
from solverFVM import solveFiniteVolume, createFVPostProzessMesh, diffusionConvectionKernel

import matplotlib.pyplot as plt
import numpy as np

swatch = pg.Stopwatch(True)

x = np.linspace(-1.0, 1.0, 41)
y = np.linspace( 0.0, 1.0, 21)
dx = x[1] - x[0]
dy = y[1] - y[0]
print(dx,dy)

grid = pg.createGrid(x=x, y=y)

# force vector per cell
f = pg.RVector(grid.cellCount(), 0.0)
#f[grid.findCell([-0.85, 0.55]).id()]=10.0

# velocity per cell [x-direction, y-direction]
vC = np.array(list(map(lambda p_: [ (2.*p_[1] *(1.0 -p_[0]**2)),
                                   (-2.*p_[0] *(1.0 -p_[1]**2))],
                       grid.cellCenters())))

vB = np.array(list(map(lambda p_: [ (2.*p_[1] *(1.0 -p_[0]**2)),
                                   (-2.*p_[0] *(1.0 -p_[1]**2))],
                       grid.boundaryCenters())))

ud0=10
udN=0

bounds = grid.findBoundaryByMarker(4)
for b in bounds:
    if b.center()[0] > 0: b.setMarker(5)

#for scheme in ['CDS', 'UDS', 'ES', 'HS', 'PS']:
print(vC.shape)

neumannBC = [[5, 0]]
Peclet = 5000 # vel/diffus
scheme='UDS'
dirichletBC = [[1, lambda b_: 1. - np.tanh(10.)],
               [2, lambda b_: 1. - np.tanh(10.)],
               [3, lambda b_: 1. - np.tanh(10.)],
               [4, lambda b_: 1. + np.tanh(10 * (2 * b_.center()[0] + 1.))],
               ]

print ("start", swatch.duration())
uGrid = solveFiniteVolume(grid, a=1./Peclet, f=f, vel=vB,
                          uBoundary=dirichletBC,
                          duBoundary=neumannBC,
                          scheme=scheme)
            

pg.showLater(1)
ax1,cb = show(grid)


drawStreams(ax1, grid, vC, coarseMesh=pg.createGrid(x=np.linspace(-1.0, 1.0, 41),
                                                    y=np.linspace(0.0, 1.0, 21)),
            )
print(swatch.duration())

#grid2, u2 = createFVPostProzessMesh(grid, u, dirichletBC)
show(grid, data=uGrid, 
     logScale=False, interpolate=1, tri=1,
     colorBar=True, axes=ax1)

#unstructured
boundary = []
boundary.append([-1.0, 0.0])
boundary.append([ 0.0, 0.0])
boundary.append([ 1.0, 0.0])
boundary.append([ 1.0, 1.0])
boundary.append([-1.0, 1.0])
    
poly = pg.Mesh(2)
nodes = [poly.createNode(b) for b in boundary]

poly.createEdge(nodes[0], nodes[1], 4) # dirichlet (inflow)
poly.createEdge(nodes[1], nodes[2], 5) # hom neumann (outflow)
poly.createEdge(nodes[2], nodes[3], 3) # hom dirichlet (isolation)
poly.createEdge(nodes[3], nodes[4], 2) # hom dirichlet (isolation)
poly.createEdge(nodes[4], nodes[0], 1) # hom dirichlet (isolation)
    
mesh = createMesh(poly, quality=34, area=0.002, smooth=[0,10])
print(mesh)

ax2,cb = show(mesh)

vC = np.array(list(map(lambda p_: [ (2.*p_[1] *(1.0 -p_[0]**2)),
                                   (-2.*p_[0] *(1.0 -p_[1]**2))],
                       mesh.cellCenters())))

vB = np.array(list(map(lambda p_: [ (2.*p_[1] *(1.0 -p_[0]**2)),
                                   (-2.*p_[0] *(1.0 -p_[1]**2))],
                       mesh.boundaryCenters())))
drawStreams(ax2, mesh, vC)

f = pg.RVector(mesh.cellCount(), 0.0)
##f[mesh.findCell([-0.75, 0.75]).id()]=1000.0


uMesh = solveFiniteVolume(mesh, a=1./Peclet, f=f, vel=vB, 
                          uBoundary=dirichletBC,
                          duBoundary=neumannBC,
                          scheme=scheme)
show(mesh, data=uMesh, 
     logScale=False, interpolate=1, tri=1,
     colorBar=True, axes=ax2)

xi = np.linspace(0.2, 0.8, 81)
uG = pg.interpolate(grid, uGrid, x=x)
uM = pg.interpolate(mesh, uMesh, x=x)
plt.figure()
plt.plot(x, uG, label='outlet ' + scheme + '(grid)')
plt.plot(x, uM, label='outlet ' + scheme + '(mesh)')
#plt.plot(x, 1. + np.tanh(10 * (2 * x + 1.)) )
plt.plot(-x, 1. + np.tanh(10 * (2 * x + 1.)), label='inlet-exact')
plt.legend()
plt.xlim([0,1])


pg.showNow()