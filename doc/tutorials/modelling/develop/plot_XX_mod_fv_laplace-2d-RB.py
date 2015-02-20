#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import pygimli as pg
import pygimli.solver as solver
from pygimli.viewer import show
from pygimli.mplviewer import drawMesh
from pygimli.meshtools import createMesh

import matplotlib.pyplot as plt
import numpy as np

from solverFVM import solveFiniteVolume, createFVPostProzessMesh

def mixedBC(boundary, userData):
    sourcePos = userData['sourcePos']
    k = userData['k']
    r1 = boundary.center() - sourcePos
    # Mirror on surface at depth=0
    r2 = boundary.center() - pg.RVector3(1.0, -1.0, 1.0) * sourcePos
    r1A = r1.abs()
    r2A = r2.abs()

    n = boundary.norm()
    # need rho here !!!!!!!!!!!!!!!!!!!!!!!!!!!1
 
    if k == 0:
        return ((r2A * r2A) * abs(r1.dot(n)) / r1A + \
                (r1A * r1A) * abs(r2.dot(n)) / r2A) / \
                (r2A * r1A * (r1A + r2A))
                   
    elif r1A > 1e-12 and r2A > 1e-12:
        return k * ((r1.dot(n)) / r1A * pg.besselK1(r1A * k) +
                    (r2.dot(n)) / r2A * pg.besselK1(r2A * k)) / \
        (pg.besselK0(r1A * k) + pg.besselK0(r2A * k))
    else:
        return 0.

def pointSource(cell, f, userData):
    sourcePos = userData['sourcePos']

    if cell.shape().isInside(sourcePos):
        f.setVal(cell.N(cell.shape().rst(sourcePos)), cell.ids())
        
grid = pg.createGrid(x=np.linspace(-10.0, 10.0, 50.),
                     y=np.linspace(-15.0,   .0, 50.))

#a = pg.RVector(grid.cellCount(),1)
a = np.ones(grid.cellCount())
for c in grid.cells():
    if c.center()[1] < -5:
        a[c.id()] = 0.1        

#grid = grid.createH2()
#grid = grid.createP2()

sourcePosA = [-5.0, -4.0]
sourcePosB = [ 5.0, -4.0]

neumannBC = [[1, mixedBC], #left boundary
             [2, mixedBC], #right boundary
             ]

dirichletBC = [[4, 0]] #bottom boundary
k = 0.1

uFE = solver.solveFEM(grid, a=a, b=k*k, f=pointSource,
                      duBoundary=neumannBC,
                      uBoundary=dirichletBC,
                      userData={'sourcePos': sourcePosA, 'k': k},
                      verbose=True)

pg.showLater(1)
ax1, cb = show(grid, uFE, 
               cMin=0, cMax=1, nLevs=10, colorBar=True)
drawMesh(ax1, grid)


f = pg.RVector(grid.cellCount(), 0)
f[grid.findCell(sourcePosA).id()]=1.0/grid.findCell(sourcePosA).size()

uFV = solveFiniteVolume(grid, a=a, f=f, 
                        duBoundary=neumannBC,
                        uBoundary=dirichletBC,
                        userData={'sourcePos': sourcePosA, 'k': k},
                        verbose=True)

#print('FVM:', swatch.duration(True))

ax2, cb = show(grid, uFV, interpolate=1, tri=1,
               cMin=0, cMax=1, nLevs=10, colorBar=True)
drawMesh(ax2, grid)

pg.showNow()
