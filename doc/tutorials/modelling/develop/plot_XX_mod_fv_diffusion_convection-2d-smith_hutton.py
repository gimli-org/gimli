#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import pygimli as pg
import pygimli.solver as solver
from pygimli.viewer import show, showMesh
from pygimli.mplviewer import drawMesh, drawModel, drawField, drawStreamLines2
from pygimli.meshtools import createMesh
from solverFVM import solveFiniteVolume, createFVPostProzessMesh, diffusionConvectionKernel

import matplotlib.pyplot as plt
import numpy as np

   
x = np.linspace(-1.0, 1.0, 41)
y = np.linspace( 0.0, 1.0, 21)
dx = x[1] - x[0]
dy = y[1] - y[0]
print(dx,dy)

grid = pg.createGrid(x=x, y=y)
N = grid.cellCount()

# force vector per cell
f = pg.RVector(grid.cellCount(), 0.0)
f[40*15+3]=10.0
# diffusions coefficient
a = pg.RVector(grid.cellCount(), 0.1)

# velocity per cell [x-direction, y-direction]

xN = pg.x(grid.positions())
yN = pg.y(grid.positions())
vN = np.vstack(((2.*yN *(1.0 -xN*xN)),
               (-2.*xN *(1.0 -yN*yN)))).T

xC = pg.x(grid.cellCenters())
yC = pg.y(grid.cellCenters())
vC = np.vstack((( 2.*yC *(1.0 -xC*xC)),
                (-2.*xC *(1.0 -yC*yC)))).T

ud0=10
udN=0

ax1,cb = show(grid, showLater=True)

swatch = pg.Stopwatch(True)
drawStreamLines2(ax1, grid, vC)

print(swatch.duration())

bounds = grid.findBoundaryByMarker(4)
for b in bounds:
    if b.center()[0]>0: b.setMarker(0)
    
bounds[7].setMarker(5)

for b in bounds:
    ax1.plot(b.center()[0], b.center()[1], 'o')


#for scheme in ['CDS', 'UDS', 'ES', 'HS', 'PS']:
print(vC.shape)
u = solveFiniteVolume(grid, a=a*0.1, f=f, v=vC, 
                      uBoundary=[[1,0],[2,0],[3,0],[4,0],[5,1]],
                      scheme='PS')
                           
show(grid, data=u, 
     cMin=0, cMax=0.3,
     logScale=False,
     colorBar=True, showLater=True, axes=ax1)


#unstructured
boundary = []
boundary.append([-1.0, 0.0])
boundary.append([ 0.0, 0.0])
boundary.append([ 1.0, 0.0])
boundary.append([ 1.0, 1.0])
boundary.append([-1.0, 1.0])
    
poly = pg.Mesh(2)
nodes = [poly.createNode(b) for b in boundary]

poly.createEdge(nodes[0], nodes[1], 2) # dirichlet (inflow)
poly.createEdge(nodes[1], nodes[2], 0) # hom neumann (outflow)
poly.createEdge(nodes[2], nodes[3], 1) # hom dirichlet (isolation)
poly.createEdge(nodes[3], nodes[4], 1) # hom dirichlet (isolation)
poly.createEdge(nodes[4], nodes[0], 1) # hom dirichlet (isolation)

for b in poly.boundaries():
    print(b.center(), b.marker())
    
mesh = createMesh(poly, quality=34, area=0.002, smooth=[0,10])

ax2,cb = show(mesh, showLater=True)

xC = pg.x(mesh.cellCenters())
yC = pg.y(mesh.cellCenters())
vC = np.vstack((( 2.*yC *(1.0 -xC*xC)),
                (-2.*xC *(1.0 -yC*yC)))).T
drawStreamLines2(ax2, mesh, vC)

bounds = mesh.findBoundaryByMarker(1,99)

for i, b in enumerate(bounds):
    ax2.plot(b.center()[0], b.center()[1], 'o')
    if b.center()[1] == 0.0 and abs(b.center()[0]+0.7) < 0.02:
        b.setMarker(5)
    
u = solveFiniteVolume(mesh,
                      a=pg.RVector(mesh.cellCount(), 0.01),
                      f=pg.RVector(mesh.cellCount(), 0.0),
                      v=vC, 
                      uBoundary=[[1,0],[2,0],[3,0],[4,0],[5,1]],
                      scheme='PS')

show(mesh, data=u, 
     cMin=0, cMax=0.3,
     logScale=False,
     colorBar=True, showLater=True, axes=ax2)


plt.show()
