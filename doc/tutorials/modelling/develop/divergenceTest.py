#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import sys

import pygimli as pg
from pygimli.viewer import show

import matplotlib.pyplot as plt
import numpy as np

from solverFVM import cellDataToCellGrad, cellDataToCellGrad2
from solverFVM import cellDataToBoundaryGrad

from fipy.meshes import Grid2D
from fipy.variables.cellVariable import CellVariable


def divergenceCell(c, F):
    ret = 0
    for bi in range(c.boundaryCount()):
        b = pg.findBoundary(c.boundaryNodes(bi))
        # print(b.norm(c).dot(F[b.id()]))
        ret += b.norm(c).dot(F[b.id()]) * b.size()

    return ret


def divergence(mesh, F):
    div = pg.RVector(mesh.cellCount())

    for c in mesh.cells():
        div[c.id()] = divergenceCell(c, F)

    return div

grid = pg.createGrid(x=np.arange(3. + 1), y=np.arange(2. + 1))
pot = np.arange(3. * 2)
print(grid, pot)

plt.ion()
show(grid, pot)

pN = pg.cellDataToPointData(grid, pot)

ax, cbar = show(grid, pN)
ax, cbar = show(grid, axes=ax)

vF = np.zeros((grid.boundaryCount(), 2))
cellGrad = cellDataToCellGrad(grid, pot)
print(cellGrad)
cellGrad = cellDataToCellGrad2(grid, pot)
print(cellGrad)

boundGrad = cellDataToBoundaryGrad(grid, pot)
boundGrad2 = cellDataToBoundaryGrad(grid, pot, 1)

for c in grid.cells():
    gr = cellGrad[c.id()]
#    gr = c.grad(c.center(), pN)
    ax.arrow(c.center()[0], c.center()[1], gr[0], gr[1])

for b in grid.boundaries():
    gr = boundGrad[b.id()]
#    gr = c.grad(c.center(), pN)
    ax.arrow(b.center()[0], b.center()[1], gr[0], gr[1], color='red')

for b in grid.boundaries():
    gr = boundGrad2[b.id()]
#    gr = c.grad(c.center(), pN)
    ax.arrow(b.center()[0], b.center()[1], gr[0], gr[1], color='green')

ax.set_xlim([-0.5, 3.5])
ax.set_ylim([-0.5, 3.1])
# F = grid.grad(pot)
# print(F)
print("div:", divergence(grid, boundGrad))

sys.path.append('/home/carsten/src/fipy')
mesh = Grid2D(nx=3, ny=2)
val = np.arange(3 * 2.)
var = CellVariable(mesh=mesh, value=val)
print(var.faceGrad._calcValueNoInline())
print('-____________')

# print(var.faceGrad)
print(var.faceGrad.divergence)
# print(var.__pos__)
# [ 4.  3.  2. -2. -3. -4.]

ax, cbar = show(grid, np.array(var))
show(grid, axes=ax)
X, Y = mesh.getCellCenters()
gr = var.grad.numericValue
m = 1
for i in range(len(X)):
    ax.arrow(X[i], Y[i], gr[0][i], gr[1][i])


gr = var.faceGrad
# gr = var.faceGrad._calcValueInline()
X, Y = mesh.getFaceCenters()
m = 1
for i in range(len(X)):
    ax.arrow(X[i], Y[i], gr[0][i], gr[1][i], color='red')


# from fipy.variables.faceGradVariable import _FaceGradVariable
# gr = _FaceGradVariable(var)
# for i in range(len(X)):
#    ax.arrow(X[i], Y[i], gr[0][i], gr[1][i], color='green')


ax.set_xlim([-0.5, 3.5])
ax.set_ylim([-0.5, 3.1])


plt.ioff()
plt.show()
