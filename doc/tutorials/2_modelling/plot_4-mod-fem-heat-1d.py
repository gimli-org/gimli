#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""

Heat equation in 1D
-------------------

Assume isotropic and homogeneous heat equation in one dimension:

<<<<<<< Updated upstream
.. math::

    \Delta u(t,x) - check(-) \frac{\partial u(t,x)}{\partial t} & = f(t,x)\\
    u(0,x) & = \sin(\pi x)\in x=\Omega \\
* Handling of time discretization

As showcase we assume the homogeneous heat equation on isotropic and 
homogeneous media in one dimension:
>>>>>>> Stashed changes

We will solve this for :math:`(t,x) \in [0,1]
\text{s} \times \Omega=[0,1]\text{m}`
temporal :math:`k=0.04\text{s}` & spatial discretization :math:`h=0.1\text{m}`

See: :py:mod:`pygimli.viewer`

Discussion: FEM is less suited for problems with piece-wise (element) constant
solutions, because linear shape functions demand twice differentiable solution.
For diffusion and wave equation with partially starting gradients = 0 you can
obtain numeric undulations (acausal overshoots) caused by the shape functions.

<<<<<<< Updated upstream
=======
We will solve for temperature :math:`u(t,x)` on the one dimensional 
domain :math:`\Omega = x = [0, 1]\text{m}` for a time interval 
:math:`t \in [0,1] \text{s}`
>>>>>>> Stashed changes
"""

import pygimli as pg
import pygimli.solver as solver
import matplotlib.pyplot as plt
import numpy as np

<<<<<<< Updated upstream
grid = pg.createGrid(x=np.linspace(0.0, 1.0, 100))
times = np.arange(0, 1.0, 0.04)
=======
x=np.linspace(0.0, 1.0, 11)
grid = pg.createGrid(x)
>>>>>>> Stashed changes

dirichletBC = [[1, 0],  # top
               [2, 0]]  # bottom

probeID = int(grid.nodeCount() / 2)

###############################################################################
<<<<<<< Updated upstream
# For this case we have an analytical solution:
=======
# The time discretization is a simple array

times = np.arange(0, 1, 0.002)

###############################################################################
# We plot the exact solution as reference solution

#plt.plot(times, uAna(times, grid.node(probeID).pos()[0]), label='exact')

###############################################################################
#For the numerical solution we review the main equation in a time discrete view
#for the Laplace operator in one dimension:
#
#.. math::
#
#  \frac{u(t+h,x)-u(t,x)}{h} - \frac{\partial^2 u(t,x)}{\partial x^2} = 0

###############################################################################
# with the time discretization step width :math:`h`.
>>>>>>> Stashed changes
#
# .. math::
#
#     u(t,x) = \e^{-\pi^2 t} \sin(\pi x)
#
#
<<<<<<< Updated upstream

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
=======
#  \mathbf{M}_{ij} = <u_i, u_j>&= \int_{\Omega} u_i u_j \qquad\text{Mass element matrix} \\
#  \mathbf{S} &= \int \nabla u \nabla v \qquad\text{Striffness matrix} 

###############################################################################
# .. warning::
#   TODO We need to explain these matrices in a different tutorial. Clean 
#   this when done

h = min(grid.cellSizes())

S = solver.createStiffnessMatrix(grid)
S_FD = S*(1./h)
M = solver.createMassMatrix(grid)
I = solver.identity(grid.nodeCount())

dirichletBC = [[1, 0],  # top
               [2, 0]]  # bottom

boundUdir = solver.parseArgToBoundaries(dirichletBC, grid)

solver.assembleDirichletBC(S, boundUdir)
solver.assembleDirichletBC(S_FD, boundUdir)
solver.assembleDirichletBC(M, boundUdir)

import scipy
import scipy.sparse.linalg
lam = scipy.sparse.linalg.eigsh(pg.utils.sparseMatrix2csr(S_FD), 
                                k=S.rows()-2)[0]
tau = 2/max(lam) # stable for FD

lam = scipy.sparse.linalg.eigsh(pg.utils.sparseMatrix2csr(S), 
                                k=S.rows()-2)[0]
tau = (2*h)**2/max(lam)  # stable for FE

#tau *= h
#tau = 0.1 / (2/h)**2

times = np.arange(0, 0.5, tau)
print('max lambda:', max(lam), '/' , (2/h)**2, 'h:', h, 'tau:', tau, 'nTau', len(times))


u = np.zeros((len(times), grid.nodeCount()))
u[0] = np.sin(np.pi * pg.x(grid))
u[0, 0] = 0.0
u[0,-1] = 0.0

print('c:', pg.solver.checkCFL(times, grid, 1))
print('dt:', tau, 'dx:', h, 'c:', 1. * tau / h)
print('h:', h, 'tau:', tau, 'tau: < ', 2*h**2, 'tau2: < ', 2/max(lam))

# pg.solver.showSparseMatrix(S, full=True)
# pg.solver.showSparseMatrix(M, full=True)
pg.tic()
SFD = (I - S_FD * tau)
for n in range(1, len(times)):
    #u[n] = u[n-1] - (S * (tau/h)) * u[n-1]
    u[n] = SFD * u[n-1]
pg.toc('FD Solution')


#plt.plot(pg.x(grid), np.log10(u[-1]), '-.')
#plt.plot(times, u[:, probeID], 'x', label='explicit FD')

pg.tic()
ut = pg.RVector(grid.nodeCount(), 0.0)
solve = pg.LinSolver(M)
for n in range(1, len(times)):
    # M * u[n] = M * u[n-1] + h * S * u[n-1]
    b = (M - S * tau) * u[n-1]
    u[n] = solve.solve(b)
pg.toc('FE Solution')

# for i in np.linspace(0, len(u)-1, 10):
#     plt.plot(pg.x(grid), np.log10(u[int(i)]))

# #plt.ylim(-6., 1.0)
# pg.wait()
plt.plot(times, u[:, probeID], label='explicit')

for n in range(1, len(times)):
    # (M + h * S ) u[n] = M * u[n-1]
    b = M * u[n - 1]
    A = M + S * tau

    solve = pg.LinSolver(A)
    u[n] = solve.solve(b)

plt.plot(times, u[:, probeID], label='implicit')

plt.show()
>>>>>>> Stashed changes

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
