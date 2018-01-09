#!/ussr/bin/env python
# -*- coding: utf-8 -*-
r"""

Modelling with Boundary Conditions
----------------------------------

We use the preceding example (Poisson equation on the unit square) but want to
specify different boundary conditions on the four sides.

Again, we first import numpy and pygimli, the solver and post processing
functionality.
"""

import numpy as np
import pygimli as pg

from pygimli.solver import solve
from pygimli.viewer import show
from pygimli.mplviewer import drawStreams

###############################################################################
# We create 21 x 21 node grid to solve on.

grid = pg.createGrid(x=np.linspace(-1.0, 1.0, 21),
                     y=np.linspace(-1.0, 1.0, 21))

###############################################################################
# We start considering inhomogeneous Dirichlet boundary conditions (BC).
#
# There are different ways of specifying BCs. They can be maps from markers to
# values, explicit functions or implicit (lambda) functions.
#
# The boundary 1 (top) and 2 (left) are directly mapped to the values 1 and 2.
# On side 3 (bottom) a lambda function 3+x is used (p is the boundary position
# and p[0] its x coordinate. On side 4 (right) a function uDirichlet is used
# that simply returns 4 in this example but can compute anything as a function
# of the individual boundaries b.


###############################################################################
# Short test: setting single node dirichlet BC
u = solve(grid, f=1., uB=[grid.node(2), 0.])

ax, _ = pg.show(grid, u, label='Solution $u$',)
show(grid, ax=ax)

def uDirichlet(boundary):
    """Return a solution value for coordinate p."""
    return 4.0

dirichletBC = [[1, 1.0],                                     # left
               [grid.findBoundaryByMarker(2), 2.0],          # right
               [grid.findBoundaryByMarker(3),
               lambda boundary: 3.0 + boundary.center()[0]], # top
               [grid.findBoundaryByMarker(4), uDirichlet]]   # bottom

###############################################################################
# The BC are passed using the uBoundary keyword. Note that showMesh returns the
# created figure ax ax while drawMesh plots on it and it can also be used as
# a class with plotting or decoration methods.
u = solve(grid, f=1., bc={'Dirichlet': dirichletBC})

ax = show(grid, data=u, colorBar=True,
          orientation='vertical', label='Solution $u$',
          levels=np.linspace(1.0, 4.0, 17), hold=1)[0]

show(grid, ax=ax)

ax.text(0, 1.01, '$u=3+x$', ha='center')
ax.text(-1.01, 0, '$u=1$', va='center', ha='right', rotation='vertical')
ax.text(0, -1.01, '$u=4$', ha='center', va='top')
ax.text(1.02, 0, '$u=2$', va='center', ha='left',  rotation='vertical')

ax.set_title('$\\nabla\cdot(1\\nabla u)=1$')

ax.set_xlim([-1.1, 1.1])  # some boundary for the text
ax.set_ylim([-1.1, 1.1])

###############################################################################
#
# Alternatively we can define the gradients of the solution on the boundary,
# i.e., Neumann type BC. This is done with another map (marker, value) and
# passed by the keyword duBoundary.
neumannBC = [[1, -0.5],  # left
             [grid.findBoundaryByMarker(4), 2.5]]  # bottom

dirichletBC = [3, 1.0]  # top

u = solve(grid, f=0., duB=neumannBC, uB=dirichletBC)

###############################################################################
# Note that on boundary 4 (right) no BC is explicitly applied leading to
# default (natural) BC that are of homogeneous Neumann type
# :math:`\frac{\partial u}{\partial n}=0`

ax = show(grid, data=u, filled=True, colorBar=True,
          orientation='vertical', label='Solution $u$',
          levels=np.linspace(min(u), max(u), 14), hold=True)[0]

# Instead of the grid we now want to add streamlines to show the gradients of
# the solution (i.e., the flow direction).

drawStreams(ax, grid, u)

ax.text(0.0, 1.01, '$u=1$',
        horizontalalignment='center')  # top -- 3
ax.text(-1.0, 0.0, '$\partial u/\partial n=-0.5$',
        va='center', ha='right', rotation='vertical')  # left -- 1
ax.text(0.0, -1.01, '$\partial u/\partial n=2.5$',
        ha='center', va='top')  # bot -- 4
ax.text(1.01, 0.0, '$\partial u/\partial n=0$',
        va='center', ha='left', rotation='vertical')  # right -- 2

ax.set_title('$\\nabla\cdot(1\\nabla u)=0$')

ax.set_xlim([-1.1, 1.1])
ax.set_ylim([-1.1, 1.1])

pg.wait()
