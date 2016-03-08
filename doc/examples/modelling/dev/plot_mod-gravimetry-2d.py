#!/usr/bin/env python

"""
Simple gravimetry field solution

Calculate dgdz at a profile for a cylindrical heterogeneity with different approaches

"""

import numpy as np

import pygimli as pg
from pygimli.viewer import *
from pygimli.solver import *
from pygimli.meshtools import *

from pygimli.physics.gravimetry import gradUCylinderHoriz, solveGravimetry

radius = 2 # [m]
depth = 5 # [m]
x = np.arange(-20, 20.1, 1)
pnts = np.array([x, np.zeros(len(x))]).T
pos = [0, -depth]
dRho = 100

ax1 = pg.plt.subplot(2,1,1)
ax2 = pg.plt.subplot(2,1,2)

# analytical first 
ax1.plot(x, gradUCylinderHoriz(pnts, radius, dRho, pos)[:,1], label='Analytical')

# polygon 
circ = createCircle([0, -depth], radius=radius, marker=2, segments=16)
ax1.plot(x, solveGravimetry(circ, dRho, pnts, complete=False), label='Poly')

# complete cell based mesh
world = createWorld(start=[-20, 0], end=[20, -10], marker=1)
mesh = createMesh([world, circ])
dRhoC = pg.solver.parseMapToCellArray([[1, 0.0], [2, dRho]], mesh)
ax1.plot(x, solveGravimetry(mesh, dRhoC, pnts), label='Mesh')

pg.show([world,  circ], axes=ax2)

# Finite Element way
world = createWorld(start=[-200, 200], end=[200, -200], marker=1)
# add some nodes to the measurement points to increase accuracy a bit
[world.createNode(x_,0.0,  1) for x_ in x]

circ = createCircle([0, -depth], radius=radius, area=0.1, marker=2, segments=16)
mesh = createMesh([world, circ], quality=34)
mesh = mesh.createP2()
print(mesh)

#pg.show(mesh,  axes=ax2)

density=solver.parseMapToCellArray([[1, 0.0], [2, dRho]], mesh)

# Calculate gravimetry potential u
u = solve(mesh, a=1, f=density, uB=[[-2,0], [-1,0]])

# Calculate gradient of gravimetry potential du/d(x,z)
duz = np.zeros(len(pnts))

for i, p in enumerate(pnts):
    c = mesh.findCell(p)
    g = c.grad(p, u)
    
    duz[i] = -g[1] * 4 * np.pi * pg.physics.constants.GmGal# wo kommen die 4 pi her?

ax1.plot(x, duz, label='duz-FEM')
ax1.legend()

pg.wait()
