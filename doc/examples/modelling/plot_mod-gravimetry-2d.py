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

radius = 2. # [m]
depth = 5. # [m]
pos = [0., -depth]
dRho = 100

x = np.arange(-20, 20.1, 1)
pnts = np.array([x, np.zeros(len(x))]).T

ax1 = pg.plt.subplot(2, 1, 1)
ax2 = pg.plt.subplot(2, 1, 2)

# analytical first 
ax1.plot(x, gradUCylinderHoriz(pnts, radius, dRho, pos)[:,1], '-vb', label='Analytical')

# polygon 
circ = createCircle([0, -depth], radius=radius, marker=2, area=0.1, segments=16)
ax1.plot(x, solveGravimetry(circ, dRho, pnts, complete=False), label='Integration: Polygon ')

# complete cell based mesh
world = createWorld(start=[-20, 0], end=[20, -10], marker=1)
mesh = createMesh([world, circ])
dRhoC = pg.solver.parseMapToCellArray([[1, 0.0], [2, dRho]], mesh)
ax1.plot(x, solveGravimetry(mesh, dRhoC, pnts), label='Integration: Mesh')
pg.show([world,  circ], axes=ax2)

# Finite Element way
world = createWorld(start=[-200, 200], end=[200, -200], marker=1)

# add some nodes to the measurement points to increase accuracy a bit
[world.createNode(x_, 0.0,  1) for x_ in x]
mesh = createMesh([world, circ], quality=34)
mesh = mesh.createP2()

#pg.show(mesh,  axes=ax2)
density=pg.solver.parseMapToCellArray([[1, 0.0], [2, dRho]], mesh)

# Calculate gravimetry potential u
u = pg.solver.solve(mesh, a=1, f=density, uB=[[-2,0], [-1,0]])

# Calculate gradient of gravimetry potential du/d(x,z)
dudz = np.zeros(len(pnts))

for i, p in enumerate(pnts):
    c = mesh.findCell(p)
    g = c.grad(p, u)
    dudz[i] = -g[1] * 4. * np.pi * pg.physics.constants.GmGal# wo kommen die 4 pi her?

ax1.plot(x, dudz, label=r'FEM: $\frac{\partial u}{\partial z}$')
ax1.set_ylabel(r'$\frac{\partial u}{\partial z}$ [mGal]')
ax1.set_xlabel('$x$-coordinate [m]')
ax1.grid()
ax1.legend()

ax2.set_aspect(1)
ax2.set_ylabel('$z$-coordinate [m]')
ax2.set_xlabel('$x$-coordinate [m]')
ax2.set_ylim((-9,1))
ax2.plot(x, x*0,  'bv')

pg.wait()
