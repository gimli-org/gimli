#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
2D gravity synthetic modelling and inversion
============================================

Based on the GIMLi issue #532.

In the following, we will build the model, create synthetic data, and do
inversion using a depth-weighting function as outlined in the paper.

"""
# %%%
import numpy as np
import pygimli as pg
import pygimli.meshtools as mt
# from pygimli.viewer import pv
from pygimli.physics.gravimetry import GravityModelling2D

# %%
# Synthetic model and data generation
# -----------------------------------
# We create a rectangular modelling domain (50x20m) with a flat anomaly  
# in a depth of about 5m.
#

x = np.arange(-25, 25.1, .5)
pnts = np.array([x, np.zeros(len(x))]).T

world = mt.createWorld(start=[-25, 0], end=[25, -20],
                       marker=1)
rect = mt.createRectangle(start=[-6, -3.5], end=[6, -6.0],
                          marker=2, area=0.1)

geom = world + rect
pg.show(geom, markers=True)
mesh = mt.createMesh(geom, quality=33, area=0.2)

# %%%
# We initialize the forward response by passing mesh and measuring points.
# Additionally, we map a density to the cell markers to build a model vector.
#

fop = GravityModelling2D(mesh=mesh, points=pnts)
dRho = pg.solver.parseMapToCellArray([[1, 0.0], [2, 300]], mesh)
g = fop.response(dRho)

# %%%
# The model response is then plotted along with the model
#

fig, ax = pg.plt.subplots(ncols=1, nrows=2, sharex=True)
ax[0].plot(x, g)
ax[0].set_ylabel(r'$\frac{\partial u}{\partial z}$ [mGal]')
ax[0].grid()
ax[0].legend()

pg.show(mesh, dRho, ax=ax[1])
ax[1].plot(x, x*0, 'bv')
ax[1].set_xlabel('$x$-coordinate [m]')
ax[1].set_ylabel('$z$-coordinate [m]')
ax[1].set_ylim((-9, 1))
ax[1].set_xlim((-25, 25))


# %%%
# We define an absolute error and add some Gaussian noise.
#

error = 0.001
g += np.random.randn(len(g)) * error


# %%%
# For inversion, we create a new mesh from the rectangular domain and setup a
# new instance of the modelling operator.
#

mesh = mt.createMesh(world, quality=33, area=1)
fop = GravityModelling2D(mesh=mesh, points=pnts)

# %%%
# Depth weighting
# ---------------
#
# In the paper of Li & Oldenburg (1996), they propose a depth weighting of the
# constraints with the formula
#
# .. math::
#
#
#    w_z = \frac{1}{(z+z_0)^{3/2}}
#

# depth weighting
bz = np.array([b.center().z() for b in mesh.boundaries() if not b.outside()])
z0 = 5
wz = 10 / (bz+z0)**1.5

# %%%
# Inversion
# ---------
#
# The inversion is rather straightforward using the standard inversion
# framework :py:class:`pygimli.Inversion`.
#

inv = pg.Inversion(fop=fop)
inv.setRegularization(limits=[-2000, 2000], trans="Cot",
                      correlationLengths=[10, 2])
inv.setConstraintWeight(wz)
rho = inv.run(g, absoluteError=error, lam=1e4, verbose=True)

# %%%
# Visualization
# -------------
#
# For showing the model, we again use the threshold filter with a value of
# 0.02. For comparison with the synthetic model, we plot the latter as a
# wireframe.
#

fig, ax = pg.plt.subplots(ncols=1, nrows=2, sharex=True)
ax[0].plot(x, g)
ax[0].set_ylabel(r'$\frac{\partial u}{\partial z}$ [mGal]')
ax[0].grid()
ax[0].legend()

pg.show(mesh, rho, ax=ax[1], logScale=False)
ax[1].plot(x, x*0, 'bv')
ax[1].set_xlabel('$x$-coordinate [m]')
ax[1].set_ylabel('$z$-coordinate [m]')
ax[1].set_ylim((-9, 1))
ax[1].set_xlim((-25, 25))


# %%%
# References
# ----------
#
# -  Li, Y. & Oldenburg, D. (1996): 3-D inversion of magnetic data.
#    Geophysics 61(2), 394-408.
# -  Holstein, H., Sherratt, E.M., Reid, A.B.  (2007): Gravimagnetic field
#    tensor gradiometry formulas for uniform polyhedra, SEG Ext. Abstr.
#
