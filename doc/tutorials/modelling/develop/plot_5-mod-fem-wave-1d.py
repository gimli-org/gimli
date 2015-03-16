#!/ussr/bin/env python
# -*- coding: utf-8 -*-
"""
Wave equation in 1D
-------------------

Simple isotropic and homogeneous full time dependent wave equation in one dimension
is defined by:

.. math::
    
    \frac{\partial^2 u}{\partial t^2}- \Delta u & = f(t,x)\quad\in\Omega \\
    u(x,0) & = u_0 \quad\in\Omega \\
    \frac{\partial u(x,0)}{\partial t} & = u_1 \quad\in\Omega \\
    u(x,t) & = g \quad\in\partial\Omega
        
We want to find the approximate solution :math:`\arr{u} \approx u(\vec{r},t)`
with :math:`\vec{r} \in \Omega` and for :math:`t\ge 0`
We will solve this for :math:`(t,x) \in [0,1]\\u{s} \times \Omega=[0,1]\\u{m}`
time step :math:`k=0.04\\u{s}` and spatial discretization :math:`h=0.1\\u{m}`

To make our life easier we transform the problem by using the substitution
:math:`\dfrac{\partial u}{\partial t} = v`.
By doing this, we get a system of equations with vanished 2nd 
derivation of time:

.. math::

    \frac{\partial v}{\partial t} - \Delta u & = f(x,t) \\
    v - \frac{\partial u}{\partial t} & = 0 \\
    u(x,0) & = u_0 \in \Omega \\
    u(x,t) & = g \in \partial\Omega \\
    v(x,0) & = u_1 \in \Omega \\
    v(x,t) & = \frac{\partial g}{\partial t} \in \partial\Omega     

The discretization in time we perform with the generally theta scheme
(see example 3):
The superscripts :math:`{}^{(n)}` and :math:`{}^{(n-1)}` denotes 
the :math:`n`-th and :math:`(n-1)`-th time step for :math:`u,v` and :math:`f`,
e.g., :math:`u^{(n)} = u(x,nk)` for the discretization in time with
:math:`k = \Delta t > 0`

.. math::
    
    \frac{\partial v}{\partial t} & = \frac{v^{(n)}-v^{(n-1)}}{k} = 
    (1-\theta)\left(\Delta u^{(n-1)} + f^{(n-1)}\right) + 
    \theta\left(\Delta u^{(n)} + f^{(n)}\right) \\
    \frac{\partial u}{\partial t} & = \frac{u^{(n)}-u^{(n-1)}}{k} = 
    (1-\theta)v^{(n-1)} + \theta v^{(n)}

Now booth equations can be rearranged to solve find a solution for
:math:`\arr{u}` that only depends on results from prior time steps.
:math:`\arr{v}^{(n)}` for the :math:`n`-th time step in a first solution of 
the linear equation system followed by a second
calculation for :math:`\arr{u}^{(n)}`:

.. math::
    
    u^{(n)} - k^2\theta^2 \Delta u^{(n)} & = u^{(n-1)} + k v^{(n-1)} + 
    k^2\theta\left[(1-\theta)\Delta u^{(n-1)} + 
    (1-\theta)f^{(n-1)} + \theta f^{(n)}\right] \\
    v & = v^{(n-1)} + k\theta\Delta u^{(n)} +
    k\left[(1-\theta)\Delta u^{(n-1)} + 
    (1-\theta)f^{(n-1)} + \theta f^{(n)}\right]

Which can be solved in the usual way by the method of weighted residuals
(see tutorial 2 and 3).
After integration by parts and spatial discretization we get 
a system of 2 equation for :math:`\arr{u}` and :math:`\arr{v}`, respectively:

.. math::
    
    (\arr{M} + k^2\theta^2\arr{A})\arr{u}^{(n)} & = 
    \arr{M}\arr{u}^{(n-1)} + k\arr{M}\arr{v}^{(n-1)} -
    k^2\theta\left[(1-\theta)(\arr{A}\arr{u}^{(n-1)} -
    \arr{F}^{(n-1)}) - \theta \arr{F}^{(n)}\right] \\
    \arr{M}\arr{v}^{(n)} & = 
    \arr{M}\arr{v}^{(n-1)} -
    k\theta\arr{A}\arr{u}^{(n)} - 
    k\left[(1-\theta)(\arr{A} \arr{u}^{(n-1)} -
    \arr{F}^{(n-1)}) -
    \theta \arr{F}^{(n)}\right]
    
Both systems of equations can by assembled and solved for each time step.

"""


import pygimli as pg
import pygimli.solver as solver

import matplotlib.pyplot as plt
import numpy as np


dx = 0.005
x = np.arange(0, 1., dx )

"""
    Start modelling code
"""

h = 1./20.
grid = pg.createGrid(x=np.arange(-1.0, 1.0+h, h))
#grid = grid.createP2()
times = np.arange(0, 2., (h)*0.5)
print(grid)

dirichletBC = [[1, 0], # top
               [2, 0]] #bottom

u = np.zeros((len(times), grid.nodeCount()))
v = np.zeros((len(times), grid.nodeCount()))
k = times[1]-times[0]

u[0,0] = 0.0
v[0,0] = 0.0
e = np.zeros(len(times))
A = solver.createStiffnessMatrix(grid)
M = solver.createMassMatrix(grid)
F = solver.assembleForceVector(grid, 0)
theta = 0.5

for n in range(1,len(times)):
    
    """
        DRY. so wen can create temporary array for the force vector and the repeating part
    """
    tmpRhs = ((1.0-theta) * (A*u[n-1] - F) - theta * F)
    rhs = M * u[n-1] + k * M * v[n-1] - k*k * theta * tmpRhs

    """
        Create the system matrix
    """
    S = M + A * k*k * theta*theta
    
    """
        Apply boundary conditions.
    """
    if n==1:
        solver.assembleDirichletBC(S,
                            solver.parseArgToBoundaries([grid.findBoundaryByMarker(1), 1], grid),
                            rhs)
        
    else:
        solver.assembleDirichletBC(S,
                            solver.parseArgToBoundaries([grid.findBoundaryByMarker(1), 0], grid),
                            rhs)
        solver.assembleDirichletBC(S,
                            solver.parseArgToBoundaries([grid.findBoundaryByMarker(2), 0], grid),
                            rhs)
    
    """
        Solve for u
    """
    u[n] = solver.linsolve(S, rhs)
    
    """
        Solve the second equation for v 
    """
    
    rhs = M * v[n-1] - A * u[n] * k * theta - k * tmpRhs

    """
        Be aware of python's #No copy at all!#. 
        So we need to take a copy of the mass element matrix ourself to keep in
        safe environment.
    """
    S = pg.RSparseMatrix(M)
    #S = M

    if n==1:
        solver.assembleDirichletBC(S,
                            solver.parseArgToBoundaries([grid.findBoundaryByMarker(1), 1], grid),
                            rhs)
        
    else:
        solver.assembleDirichletBC(S,
                            solver.parseArgToBoundaries([grid.findBoundaryByMarker(1), 0], grid),
                            rhs)
        solver.assembleDirichletBC(S,
                            solver.parseArgToBoundaries([grid.findBoundaryByMarker(2), 0], grid),
                            rhs)
    
    v[n] = solver.linsolve(S, rhs)

    e[n] = 0.5 * (pg.dot(M * v[n], v[n]) + pg.dot(A * u[n], u[n]))
    
    print(n, times[n], "Energy", e[n], k/(h*h))

    
import matplotlib.animation as animation
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(e)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

lineFEM, = ax.plot(pg.x(grid.positions()), u[0])
ax.set_ylim(-0.1, 0.1)

def init():
    #lineANA.set_ydata(np.ma.array(x, mask=True))
    #lineFEM.set_xdata(np.ma.array(x, mask=True))
    return lineFEM 

def animate(i):
#    lineANA.set_ydata(uTana[i])  # update the data

    ut = u[i]
    #ut = np.ma.masked_where(abs(u[i]) <1e-6, u[i])

    lineFEM.set_ydata(ut)  # update the data
    return lineFEM

ani = animation.FuncAnimation(fig, animate, np.arange(0, len(times)),
                              init_func=init, interval=1, blit=False)

"""
.. image:: PLOT2RST.current_figure
    :scale: 75

"""
plt.show()