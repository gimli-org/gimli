#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
Simple gravimetric field solution.

Calculate for the gravimetric potential :math:`u`

.. math::
    \frac{\partial u}{\partial z}

at a measurement profile for a cylindrical heterogeneity with different approaches.
"""

import numpy as np
import pygimli as pg
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

###############################################################################
# Analytical solution first 
ax1.plot(x, gradUCylinderHoriz(pnts, radius, dRho, pos)[:,1], 
         '-vb', label='Analytical')

###############################################################################
# Integration for a 2D polygon after :cite:`WonBevis1987`
circ = createCircle([0, -depth], radius=radius, marker=2, area=0.1, segments=16)
ax1.plot(x, solveGravimetry(circ, dRho, pnts, complete=False),
         label='Integration: Polygon ')

###############################################################################
# Integration for complete 2D mesh after :cite:`WonBevis1987`
world = createWorld(start=[-20, 0], end=[20, -10], marker=1)
mesh = createMesh([world, circ])
dRhoC = pg.solver.parseMapToCellArray([[1, 0.0], [2, dRho]], mesh)
ax1.plot(x, solveGravimetry(mesh, dRhoC, pnts), label='Integration: Mesh')
pg.show([world,  circ], axes=ax2)

###############################################################################
# Finite Element solution for :math:`u`
world = createWorld(start=[-200, 200], end=[200, -200], marker=1)

# Add some nodes to the measurement points to increase the accuracy a bit
[world.createNode(x_, 0.0,  1) for x_ in x]
mesh = createMesh([world, circ], quality=34)
mesh = mesh.createP2()

density = pg.solver.parseMapToCellArray([[1, 0.0], [2, dRho]], mesh)
u = pg.solver.solve(mesh, a=1, f=density, uB=[[-2,0], [-1,0]])

# Calculate gradient of gravimetric potential 
# :math:`\frac{\partial u}{\partial (x,z)}`
dudz = np.zeros(len(pnts))

for i, p in enumerate(pnts):
    c = mesh.findCell(p)
    g = c.grad(p, u)
    dudz[i] = -g[1] * 4. * np.pi * pg.physics.constants.GmGal # why 4 pi here?
ax1.plot(x, dudz, label=r'FEM: $\frac{\partial u}{\partial z}$')

###############################################################################
# Finishing the plots

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
