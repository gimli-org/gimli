#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import pygimli as pg
import pygimli.solver as solver

import matplotlib.pyplot as plt
import numpy as np

from solverFVM import solveFiniteVolume

dx = 0.005
x = np.arange(0, 0.2 + dx, dx)
x0 = max(x)/2.
dt = 1e-2
t = np.arange(dt, 5.0 + dt, dt)
tempLF = 1e-4

timeProbeID = int(1.6 / dt)
spaceProbeID = int(0.15 / dx)
sourcePosID = int(x0 / dx)
u0 = np.zeros(len(x))
u0[sourcePosID] = 1.

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

uG = np.zeros((len(t), len(x)))
for i, ti in enumerate(t):
    g = solver.greenDiffusion1D(np.hstack((-x[:0:-1], x)), ti, a=tempLF)*dx
    uG[i] = np.convolve(g, u0)[len(x)-1:2*len(x)-1]

ax1.plot(x, uG[timeProbeID],
         label='Analytical at ' + str(t[timeProbeID]) + ' s')
ax2.plot(t, uG[:, spaceProbeID],
         label='Analytical at ' + str(x[spaceProbeID]) + ' m')

ax1.set_xlabel('Distance in m')
ax1.set_ylabel('Temperature difference in K')
ax1.grid()

ax2.set_xlabel('Time in s')
ax2.set_ylabel('Temperature difference in K')
ax2.grid()


"""
    Start modelling code
"""

ut = np.zeros(len(x))  # du/dt(x,0)

v = tempLF  # [(m*m)/s]
c = v * dt / (dx*dx)  # [(m*m)/s * s/(m*m)] = []
print("Courant-Friedrichs-Lewy condition:", "CLF =", c)
# r needs to be <= 0.5

uE = np.zeros((len(t), len(x)))
uE[0] = u0
uE[0][0] = 0.0
uE[0][-1] = 0.0

uI = np.array(uE)
uFEM = np.array(uE)

#L = np.diag(np.ones(len(x)-1), -1) - 2. * np.diag(np.ones(len(x))) + np.diag(np.ones(len(x)-1), 1)
Lg = solver.triDiagToeplitz(len(x), a=2.0, l=-1.0, r=-1.0,
                            start=1, end=len(x)-1)
I = solver.identity(len(x))
# L[1,0]=2 # neummann links
# L[-2,-1]=2 # neummann rechts

grid = pg.createGrid(x=x)
uFEM = np.array(uE)
dirichletBC = [[1, 0.0],  # top
               [2, 0.0]]  # bottom

A = solver.createStiffnessMatrix(grid, np.ones(grid.cellCount()))
M = solver.createMassMatrix(grid, np.ones(grid.cellCount()))

for n in range(1, len(t)):
    """
        Direct time and space discretization of the second derivatives
        # diff in space is 2nd order central
    """
    # FD explicit
#    uE[n] = uE[n-1] + c * (Lg * -uE[n-1])

    # FD implicit
    uI[n] = solver.linsolve(I + Lg * c, uI[n-1])

    # FE implicit
    S = M + A * tempLF * dt
    rhs = M * uFEM[n-1]
    solver.assembleBoundaryConditions(grid, S,
                                      rhs=rhs,
                                      boundArgs=dirichletBC,
                                      assembler=solver.assembleDirichletBC)
    uFEM[n] = solver.linsolve(S, rhs)

#ax1.plot(x, uE[timeProbeID], label='FD explicit')
#ax2.plot(t, uE[:,spaceProbeID], label='FD explicit')

ax1.plot(x, uI[timeProbeID], '.-', label='FD implicit')
ax2.plot(t, uI[:, spaceProbeID], '.-', label='FD implicit')

#ax1.plot(x, uFEM[timeProbeID], label='FE implicit')
#ax2.plot(t, uFEM[:,spaceProbeID], label='FE implicit')

uFVM = solveFiniteVolume(grid, a=tempLF, f=0.0, times=t, u0=u0,
                         uDirichlet=[0.0, 0.0])

print(len(uFVM), len(uFVM[0]))
ax1.plot(pg.x(grid.cellCenters()), uFVM[timeProbeID], '.-',
         label='FV implicit')
ax2.plot(t, uFVM[:, spaceProbeID], '.-', label='FV implicit')

uFE = solver.solvePoisson(grid, times=t, a=tempLF, theta=0.5, u0=u0,
                          uBoundary=dirichletBC)

ax1.plot(x, uFE[timeProbeID], label='FE Crank-Nicolsen')
ax2.plot(t, uFE[:, spaceProbeID], label='FE Crank-Nicolsen')

ax1.legend()
ax2.legend()


"""
.. image:: PLOT2RST.current_figure
    :scale: 75

"""
plt.show()
