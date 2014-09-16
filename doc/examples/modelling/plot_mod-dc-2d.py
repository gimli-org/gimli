#!/usr/bin/env python
# encoding: UTF-8

"""
Geoelectric in 2.5d
-------------------

Let us start with a mathematical formulation ...

.. math::

    \nabla\cdot( \sigma \nabla u ) = -I\delta(\vec{r}-\vec{r}_{\text{s}}) \in R^3

Fourier-Cosine-Transform

.. math::
    \nabla\cdot( \sigma \nabla u ) & = -I\delta(\vec{r}-\vec{r}_{\text{s}}) \in R^2 \\
    \frac{\partial u}{\partial \vec{n}} & = 0 \quad\mathrm{on}\quad\text{Surface} z=0

"""

import pygimli as pg

from pygimli.solver import solvePoisson

from pygimli.viewer import showMesh
from pygimli.mplviewer import *

"""
Maybe this is usefully. The analytical solution for one source location.
"""

def uAnalytical(p, sourcePos, k):
    r1A = (p - sourcePos).abs()
    # Mirror on surface at depth=0
    r2A = (p - pg.RVector3(1.0, -1.0, 1.0) * sourcePos).abs()

    # need rho here!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    if r1A > 1e-12 and r2A > 1e-12:
        return (pg.besselK0(r1A * k) + pg.besselK0(r2A *k)) / (2.0 * np.pi)
    else:
        return 0.

"""

Define the derivative of the analytical solution regarding the outer normal direction
:math:`\vec{n}`. So we can define the value for the Neumann type Boundary
conditions for the boundaries in the subsurface.

"""

def mixedBC(boundary, userData):
    sourcePos = userData['sourcePos']
    k = userData['k']
    r1 = boundary.center() - sourcePos
    # Mirror on surface at depth=0
    r2 = boundary.center() - pg.RVector3(1.0, -1.0, 1.0) * sourcePos
    r1A = r1.abs()
    r2A = r2.abs()

    n = boundary.norm()
    # need rho here !!!!!!!!!!!!!!!!!!!!!!!!!!!1

    if r1A > 1e-12 and r2A > 1e-12:
        return k * ((r1.dot(n)) / r1A * pg.besselK1(r1A * k) +
                    (r2.dot(n)) / r2A * pg.besselK1(r2A * k)) / \
        (pg.besselK0(r1A * k) + pg.besselK0(r2A * k))
    else:
        return 0.


"""

Define function for the current source term
:math:`\delta(x-pos), \int f(x) \delta(x-pos)=f(pos)=N(pos)`
Right hand side entries will be shape functions(pos)

"""

def pointSource(cell, f, userData):
    sourcePos = userData['sourcePos']

    if cell.shape().isInside(sourcePos):
        f.setVal(cell.N(cell.shape().rst(sourcePos)), cell.ids())

"""

"""

grid = pg.createGrid(x=np.linspace(-10.0, 10.0, 50), y=np.linspace(-15.0, .0, 50))

#grid = grid.createH2()
grid = grid.createP2()

sourcePosA = [-5.0, -4.0]
sourcePosB = [ 5.0, -4.0]

neumannBC = [[1, mixedBC], #left boundary
             [2, mixedBC], #right boundary
             [4, mixedBC]] #bottom boundary

"""

"""
k = 1e-3
u = solvePoisson(grid, a=1, b=k*k, f=pointSource,
                 duBoundary=neumannBC,
                 userData={'sourcePos': sourcePosA, 'k': k},
                 verbose=True)

u -= solvePoisson(grid, a=1, b=k*k, f=pointSource,
                  duBoundary=neumannBC,
                  userData={'sourcePos': sourcePosB, 'k': k},
                  verbose=True)


#uAna = pg.RVector(map(lambda p__: uAnalytical(p__, sourcePosA, k), grid.positions()))
#uAna -= pg.RVector(map(lambda p__: uAnalytical(p__, sourcePosB, k), grid.positions()))

#err = (1.0 -u/uAna)*100.0

#print "error min max", min(err), max(err)

ax = showMesh(grid, data=u, filled=True, colorBar=True,
              orientation='horizontal', label='Solution $u$',
              showLater=True)[0]

drawMesh(ax, grid)

"""
Instead of the grid we want to add streamlines to the plot to show the gradients
of the solution.
"""

gridCoarse = pg.createGrid(x=np.linspace(-10.0, 10.0, 30), y=np.linspace(-15.0, .0, 30))
drawStreamLines2(ax, grid, u, coarseMesh=gridCoarse, color='Black')

"""
.. image:: PLOT2RST.current_figure
    :scale: 75

"""

plt.show()
