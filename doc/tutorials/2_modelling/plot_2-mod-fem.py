#!/ussr/bin/env python
# -*- coding: utf-8 -*-
# sphinx_gallery_thumbnail_number = 2

r"""
Basics of finite element modelling
----------------------------------

This tutorial covers the first steps into finite element computation
using the *M* (Modelling) in *GIMLi*.

We will not dig into deep details about the theory of finite elements
here, as this can be found in several books, e.g.,
:cite:``Zienkiewicz1977``.

However, there is a little need for theory to understand what it means
to use the finite elements method to solve a modelling problem, so we
start with some basics.

Assuming the Poisson equation for a simple partial differention equation
that needs to be solved for the unknown scalar field
:math:`u(\mathbf{r})` within a modelling domain :math:`\mathbf{r}\in\Omega` 
for a non zero right hand side function :math:`f`.

.. math::
   
   - \Delta u & = f \quad{\mathrm{in}}\quad~\Omega\\
            u & = g \quad{\mathrm{on}}\quad\partial\Omega

The Laplace operator :math:`\Delta = \nabla\cdot\nabla` given by the divergence of the gradient
is the sum of second partial derivatives 
the of the field :math:`u(\mathbf{r})` with respect to the Cartesian coordinates
in 1d space :math:`\mathbf{r} = (x)`, in 2d :math:`\mathbf{r} = (x, y)`, or 3d
space :math:`\mathbf{r} = (x, y, z)`. 
On the boundary :math:`\partial\Omega` of the domain we want 
known values of :math:`u=g` as the so called Dirichlet boundary conditions.

An approximated solution :math:`u_h\approx u` for our numerical problem will probably
only satisfy with a rest :math:`R`: :math:`\Delta u_h + f = R`.
If we choose some weighting functions :math:`w`, we can try to minimize 
these residuum over our modelling domain.

.. math::

    \int_{\Omega} R w &= 0 \\
    \int_{\Omega} - \Delta u w & = \int_{\Omega} f w

It is preferable to eliminate the second derivative in the Laplace operator 
either due to integration by parts or by applying production rule and Gauss's law.
This leads to the so called weak formulation:

.. math::

    \int_{\Omega} \nabla u \nabla w - \int_{\partial \Omega}\mathbf{n}\nabla u w & = \int_{\Omega} f w \\
    \int_{\Omega} \nabla u \nabla w & = \int_{\Omega} f w + \int_{\partial \Omega}\frac{\partial u}{\partial\mathbf{n}} w

to be continued soon ...
"""

import pygimli as pg

import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# We create a grid for our modelling domain with equidistant spacing in x
# and y direction.

grid = pg.createGrid(x=np.linspace(-1.0, 1.0, 10),
                     y=np.linspace(-1.0, 1.0, 10))

###############################################################################
# Now we can call the solver :py:mod:`pygimli.solver.solve`  for some
# default material values and global homogeneous Dirichlet boundary conditions.

u = pg.solver.solve(grid, f=1.,
                    bc={'Dirichlet': [grid.findBoundaryByMarker(1, 5), 0.0]},
                    verbose=True)

###############################################################################
# The result is drawn with the function :py:mod:`pygimli.show`.

ax, cbar = pg.show(grid, data=u, label='P1 Solution $u$')

###############################################################################
# :py:mod:`pygimli.show` is just a shortcut for various routines that can also
# be called directly.
#
pg.mplviewer.drawMesh(ax, grid)

###############################################################################
#
# We repeat the computation with a spatially (H) refined version of the
# original grid.

gridh2 = grid.createH2()

uh = pg.solver.solve(gridh2, f=1.,
                     bc={'Dirichlet': [gridh2.findBoundaryByMarker(1, 5), 0.0]},
                     verbose=True)

ax, cbar = pg.show(gridh2, data=uh, label='H2 Solution $u$')

pg.mplviewer.drawMesh(ax, gridh2)

###############################################################################
#
# We do the same using quadratic (P) refinement, i.e. the same number of nodes.

gridp2 = grid.createP2()

up = pg.solver.solve(gridp2, f=1.,
                     bc={'Dirichlet': [gridp2.findBoundaryByMarker(1, 5), 0.0]},
                     verbose=True)
###############################################################################
# Fortunately there exist an analytical solution for this example.

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


###############################################################################
# To compare the different results the in detail we interpolate our solution
# along a probe line through the domain.
#

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

plt.show()
