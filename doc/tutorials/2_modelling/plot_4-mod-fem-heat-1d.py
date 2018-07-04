#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""

Heat equation in 1D
-------------------

Assume isotropic and homogeneous heat equation in one dimension:

.. math::

    \Delta u(t,x) - check(-) \frac{\partial u(t,x)}{\partial t} & = f(t,x)\\
    u(0,x) & = \sin(\pi x)\in x=\Omega \\
* Handling of time discretization

As showcase we assume the homogeneous heat equation on isotropic and 
homogeneous media in one dimension:

We will solve this for :math:`(t,x) \in [0,1]
\text{s} \times \Omega=[0,1]\text{m}`
temporal :math:`k=0.04\text{s}` & spatial discretization :math:`h=0.1\text{m}`

See: :py:mod:`pygimli.viewer`

"""

import pygimli as pg
import pygimli.solver as solver
import matplotlib.pyplot as plt
import numpy as np

grid = pg.createGrid(x=np.linspace(0.0, 1.0, 100))
times = np.arange(0, 1.0, 0.04)

dirichletBC = [[1, 0],  # top
               [2, 0]]  # bottom

probeID = int(grid.nodeCount() / 2)

###############################################################################
# For this case we have an analytical solution:
#
# .. math::
#
#     u(t,x) = \e^{-\pi^2 t} \sin(\pi x)
#
#

def uAna(t, x):
    return np.exp(-np.pi**2. * t) * np.sin(np.pi * x)

plt.plot(times, uAna(times, grid.node(probeID).pos()[0]), label='Analytical')

dof = grid.nodeCount()
u = np.zeros((len(times), dof))
u[0, :] = list(map(lambda r: np.sin(np.pi * r[0]), grid.positions()))

dt = times[1] - times[0]
A = solver.createStiffnessMatrix(grid, np.ones(grid.cellCount()))
M = solver.createMassMatrix(grid, np.ones(grid.cellCount()))

ut = pg.RVector(dof, 0.0)
rhs = pg.RVector(dof, 0.0)
b = pg.RVector(dof, 0.0)
theta = 0

boundUdir = solver.parseArgToBoundaries(dirichletBC, grid)

for n in range(1, len(times)):
    b = (M - A * dt) * u[n - 1] + rhs * dt
    S = M

    solver.assembleDirichletBC(S, boundUdir, rhs=b)

    solve = pg.LinSolver(S)
    solve.solve(b, ut)

    u[n, :] = ut

plt.plot(times, u[:, probeID], label='Explicit Euler')

theta = 1

for n in range(1, len(times)):
    b = (M + A * (dt*(theta - 1.0))) * u[n-1] + \
        rhs * (dt*(1.0 - theta)) + \
        rhs * dt * theta

    b = M * u[n-1] + rhs * dt

    S = M + A * dt

    solver.assembleDirichletBC(S, boundUdir, rhs=b)

    solve = pg.LinSolver(S)
    solve.solve(b, ut)

    u[n, :] = ut

plt.plot(times, u[:, probeID], label='Implicit Euler')

u = solver.solve(grid, times=times, theta=0.5,
                 u0=lambda node: np.sin(np.pi * node[0]),
                 bc={'Dirichlet':dirichletBC})

plt.plot(times, u[:, probeID], label='Crank-Nicolson')

plt.xlabel("t[s] at x = " + str(round(grid.node(probeID).pos()[0], 2)))
plt.ylabel("u")
plt.ylim(0.0, 1.0)
plt.legend()
plt.grid()

###############################################################################
# Explicit Euler scheme is unstable at progressing time.

plt.show()
