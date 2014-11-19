#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import pygimli as pg
import pygimli.solver as solver
from pygimli.viewer import show
from pygimli.mplviewer import drawMesh, drawStreamLines2
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
        
grid = pg.createGrid(x=np.linspace(-10.0, 10.0, 50), y=np.linspace(-15.0, .0, 50))

#grid = grid.createH2()
#grid = grid.createP2()

sourcePosA = [-5.0, -4.0]
sourcePosB = [ 5.0, -4.0]

neumannBC = [[1, mixedBC], #left boundary
             [2, mixedBC], #right boundary
             ]

dirichletBC = [[4, 0]] #bottom boundary
k = 0
u = solver.solvePoisson(grid, a=1, b=k*k, f=pointSource,
     #                   duBoundary=neumannBC,
                        uBoundary=dirichletBC,
                        userData={'sourcePos': sourcePosA, 'k': k},
                        verbose=True)


ax1, cb = show(grid, u, colorBar=True, filled=True, showLater=True, nLevs=10)
drawMesh(ax1, grid)


f = pg.RVector(grid.cellCount(), 0)

#f[grid.findCell(sourcePosA).id()]=10.0
f[grid.findCell(sourcePosA).id()]=1.0/grid.findCell(sourcePosA).size()

u = solveFiniteVolume(grid, a=1, f=f, 
                      uBoundary=dirichletBC)

#print('FVM:', swatch.duration(True))

ax2, cb = show(grid, u, 
               colorBar=True, logScale=False, filled=True, showLater=True, nLevs=10)
drawMesh(ax2, grid)


mesh2, boundSortIdx = createFVPostProzessMesh(grid)


ax3,cbar = show(mesh2, data=pg.cat(u, np.zeros(len(boundSortIdx))),
                filled=True, colorBar=True, showLater=True, nLevs=10)

#print('---:', swatch.duration(True))
show(mesh2, axes=ax3, showLater=True)









#mesh2, boundSortIdx = createFVPostProzessMesh(mesh, bounds)

#print('---:', swatch.duration(True))
#ax3,cbar = showMesh(mesh2, data=pg.cat(u, uDirichlet[boundSortIdx]),
                    #filled=True, 
                    #cMin=0, cMax=10, colorBar=True, showLater=True)

#drawMesh(ax3, mesh2)



##ax3,cbar = showMesh(mesh,
                    ##data=np.array(list(map(lambda p_: np.sin(np.arctan2(p_[1],p_[0]))/p_.abs(), mesh.cellCenter()))),
                    ##cMin=0, cMax=10, logScale=False,
                    ##showLater=True)
##drawMesh(ax3, mesh)

plt.show()
#drawMesh(ax, grid)