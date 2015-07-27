#!/ussr/bin/env python
# -*- coding: utf-8 -*-
"""

Modelling
---------

This is the first step for the modelling tutorial where we actually use
finite elements computation.

We will not go in deep detail about the finite elements theory here,
as this can be found in several books, e.g., :cite:`Zienkiewicz1977`

In this modelling tutorial we just want to solve some simple problems
to show how the *M* (Modelling) in *GIMLi* works.

We start with a simple elliptic partial differential equation and
with zero boundary values, but a nonzero right hand side.

.. math::

    \nabla\cdot(A \cdot \nabla u ) + B u + C & =
    0 \quad{\mathrm{in}}\quad\Omega \\
    \alpha u + \beta \frac{\partial u }{\partial n} & =
    0 \quad{\mathrm{on}}\quad\partial\Omega

Let :math:`A=1,\,B=0\,` and :math:`C = 1` we get the simplest Poisson equation:

.. math::

    - \Delta u & = 1 \quad{\mathrm{in}}\quad\Omega\\
             u & = 0 \quad{\mathrm{on}}\quad\partial\Omega\\


Model domain is the unit square: :math:`\Omega=[-1, 1]^2`

We start by importing the pygimli package.
"""

import pygimli as pg

import numpy as np
import matplotlib.pyplot as plt

"""
We create a grid for our modelling domain with equidistant spacing in x
and y direction.
"""

grid = pg.createGrid(x=np.linspace(-1.0, 1.0, 10),
                     y=np.linspace(-1.0, 1.0, 10))

"""
Now we can call the solver :py:mod:`pygimli.solver.solve`  for some
default material values and global homogeneous Dirichlet boundary conditions.
"""

u = pg.solver.solve(grid, f=1.,
                    uB=[grid.findBoundaryByMarker(1, 5), 0.0],
                    verbose=True)

"""
The result is drawn with the function :py:mod:`pygimli.show`. 
"""

ax, cbar = pg.show(grid, data=u, colorBar=True, label='P1 Solution $u$', hold=1)

"""
:py:mod:`pygimli.show` is just a shortcut for various routines that can also be called directly.

"""
pg.mplviewer.drawMesh(ax, grid)

"""
.. image:: PLOT2RST.current_figure
    :scale: 50

We repeat the computation with a spatially (H) refined version of the original
grid.
"""

gridh2 = grid.createH2()

uh = pg.solver.solve(gridh2, f=1.,
                     uB=[gridh2.findBoundaryByMarker(1, 5), 0.0],
                     verbose=True)

ax, cbar = pg.show(gridh2, data=uh, colorBar=True, label='H2 Solution $u$')

pg.mplviewer.drawMesh(ax, gridh2)

"""
.. image:: PLOT2RST.current_figure
    :scale: 50

We do the same using a quadratic (P) refinement, i.e. the same number of nodes.
"""

gridp2 = grid.createP2()

up = pg.solver.solve(gridp2, f=1.,
                     uB=[gridp2.findBoundaryByMarker(1, 5), 0.0],
                     verbose=True)

"""
Fortunately there exist an analytical solution for this example.
"""

def uAna(r):
    x = r[0]
    y = r[1]

    ret = 0
    for k in range(1, 151, 2):
        kp = k*np.pi
        s = np.sin(kp * (1. + x)/2) / (k**3 * np.sinh(kp)) * \
            (np.sinh(kp * (1. + y)/2) + np.sinh(kp * (1. - y)/2))
        ret += s
    return (1. - x**2)/2 - 16./(np.pi**3) * ret


"""
To compare the different results the in detail we interpolate our solution
along a probe line through our domain.

"""

x = np.linspace(-1.0, 1.0, 100)

probe = np.zeros((len(x), 3))
probe[:, 0] = x

uH1 = pg.interpolate(srcMesh=grid, inVec=u, destPos=probe)
uH2 = pg.interpolate(srcMesh=gridh2, inVec=uh, destPos=probe)
uP2 = pg.interpolate(srcMesh=gridp2, inVec=up, destPos=probe)

plt.figure()
plt.plot(x, np.array(list(map(uAna, probe))), 'black', linewidth=2,
         label='analytical')
plt.plot(x, uH1, label='linear (H1)')
plt.plot(x, uH2, label='linear (H2)')
plt.plot(x, uP2, label='quadratic (P2)')

plt.xlim([-0.4, 0.4])
plt.ylim([0.25, 0.3])
plt.legend()


"""
.. image:: PLOT2RST.current_figure
    :scale: 50

"""

plt.show()
