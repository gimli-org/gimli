#!/ussr/bin/env python
# -*- coding: utf-8 -*-
"""
    
Heat equation in 1D
-------------------

Isotropic and homogeneous heat equation in one dimension with test case:

.. math::

    \Delta u(t,x) + \frac{\partial u(t,x)}{\partial t} & = f(t,x)\\
    u(0,x) & = \sin(\pi x)\in x=\Omega \\
    u(t,x) & = 0 \in x=\partial\Omega

We will solve this for :math:`(t,x) \in [0,1]\text{s} \times \Omega=[0,1]\text{m}`
time step :math:`k=0.04\text{s}` and spatial discretization :math:`h=0.1\text{m}`

"""

import pygimli as g
from pygimli.solver import solvePoisson
import matplotlib.pyplot as plt
import numpy as np

grid = g.createGrid(x=np.linspace(0.0, 1.0, 10))
times = np.arange(0, 1.0, 0.04)
            
dirichletBC = [[1, 0], # top
               [2, 0]] #bottom

probeID = grid.nodeCount() / 2

"""
Fortunately we have an analytical solution:

.. math::

    u(t,x) = \e^{-\pi^2 t} \sin(\pi x)
    

"""

def uAna(t, x):
    return np.exp(-np.pi**2. * t) * np.sin(np.pi * x)

plt.plot(times, uAna(times, grid.node(probeID).pos()[0]), label='Analytical')

"""
    

"""

u = solvePoisson(grid, times=times, theta=0.0,
                 u0=lambda r: np.sin(np.pi * r[0]),
                 uBoundary=dirichletBC)
plt.plot(times, u[:, probeID], label='explicit Euler')

"""
    

"""
u = solvePoisson(grid, times=times, theta=0.5,
                 u0=lambda r: np.sin(np.pi * r[0]),
                 uBoundary=dirichletBC)
plt.plot(times, u[:, probeID], label='Crank-Nicolson')

"""
    

"""
u = solvePoisson(grid, times=times, theta=1.0,
                 u0=lambda r: np.sin(np.pi * r[0]),
                 uBoundary=dirichletBC)
plt.plot(times, u[:, probeID], label='implicit Euler')

plt.xlabel("t[s] at x = " + str(round(grid.node(probeID).pos()[0],2)))
plt.ylabel("u")
plt.ylim(0.0, 1.0)
plt.xlim(0.0, 0.5)
plt.legend()
plt.grid()

plt.show()

"""
.. image:: PLOT2RST.current_figure
    :scale: 75

Explicit Euler scheme is unstable at hight times.   
"""

