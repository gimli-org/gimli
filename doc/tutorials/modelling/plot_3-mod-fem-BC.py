#!/ussr/bin/env python
# -*- coding: utf-8 -*-
"""
    
Modelling with Boundary Conditions
----------------------------------

This is the first step for the modelling tutorial where we actually use finite elements to compute something. 

We will not go in deep detail about the Finite Elements Theory here, this can be found in several books, e.g., 
:cite:`Zienkiewicz1977`

We just want to solve some problems to show how the *M* in *GIMLi* works.

We will solve a simple version of Poisson's equation with zero boundary values, but a nonzero right hand side:

Let us start with a mathematical formulation ...

.. math::

    \nabla\cdot( A \cdot \nabla u ) + B u + C = 0
   
.. math::

    - \Delta u & = 1 \quad{\mathrm{in}}\quad\Omega\\
               u & = u_d \quad{\mathrm{on}}\quad\partial\Omega\\
     \frac{\partial u}{\partial \vec{n}} & = u_n \quad{\mathrm{on}}\quad\partial\Omega\\


We will solve this equation on the unit square: :math:`\Omega=[-1, 1]^2`

As usually, the library have to be imported first.

**Define term Dirichlet boundaries**

"""

import pygimli as g

from pygimli.solver import solvePoisson

'''

'''
import numpy as np
import matplotlib.pyplot as plt

from pygimli.viewer import showMesh
from pygimli.mplviewer import *


grid = g.createGrid(x=np.linspace(-1.0, 1.0, 50), y=np.linspace(-1.0, 1.0, 50))

#grid = grid.createP2()

"""
Define a list of pairs for boundaries and potential values for the Dirichlet boundaries. 
Definition of the used boundaries are either by the marker directly or by a given list of boundaries.
The value can be defined directly or by a given function.

"""

def uDirichlet(b):
    '''
        Return a solution value for coordinate p.
    '''
    return 4.0

dirichletBC = [[1, 1.0], # top
               [grid.findBoundaryByMarker(2), 2.0], # left
               [grid.findBoundaryByMarker(3), lambda p: 3.0 + p[0]], # bottom
               [grid.findBoundaryByMarker(4), uDirichlet]] # right

"""

"""
u = solvePoisson(grid, f=1.,
                 uBoundary=dirichletBC)

ax = showMesh(grid, data=u, filled=True, colorBar=True,
              orientation='vertical', label='Solution $u$',
              levels=np.linspace(1.0, 4.0, 17), showLater=True)[0]

drawMesh(ax, grid)

ax.text( 0.0, 1.02, '$u=1$')
ax.text(-1.08, 0.0, '$u=2$', rotation='vertical')
ax.text( 0.0, -1.08, '$u=3+x$')
ax.text( 1.02, 0.0, '$u=4$', rotation='vertical')


ax.set_title('$\\nabla\cdot(1\\nabla u)=1$')

ax.set_xlim([-1.1, 1.1])
ax.set_ylim([-1.1, 1.1])
   
"""
.. image:: PLOT2RST.current_figure
    :scale: 75


We can define how the gradients of the solution have to be on the boundaries, 
i.e., Boundary conditions of Neumann type.

"""
neumannBC = [[2, -0.5], # left
             [grid.findBoundaryByMarker(3), 2.4]] # bottom

dirichletBC = [1, 1.0] # top

"""
On boundary 4, the top side, the default or natural boundary conditions are applied 
:math:`\frac{\partial u}{\partial n}=0`
"""
              
u = solvePoisson(grid, f=0.,
                 duBoundary=neumannBC,
                 uBoundary=dirichletBC)

ax = showMesh(grid, data=u, filled=True, colorBar=True,
              orientation='vertical', label='Solution $u$',
              levels=np.linspace(min(u), max(u), 14), showLater=True)[0]

"""
Instead of the grid we want to add streamlines to the plot to show the gradients
of the solution.
"""

drawStreamLines(ax, grid, u, nx=25, ny=25, color='Black')

#drawMesh(ax, grid)
ax.text(0.0, 1.02, '$u=1$',
        horizontalalignment='center' )
ax.text(-1.08, 0.0, '$\partial u/\partial n=-0.5$',
        verticalalignment='center', rotation='vertical')
ax.text(0.0, -1.08, '$\partial u/\partial n=2.5$',
        horizontalalignment='center')
ax.text(1.02, 0.0, '$\partial u/\partial n=0$',
        verticalalignment='center', rotation='vertical')


ax.set_title('$\\nabla\cdot(1\\nabla u)=0$')

ax.set_xlim([-1.1, 1.1])
ax.set_ylim([-1.1, 1.1])

"""
.. image:: PLOT2RST.current_figure
    :scale: 75

"""

P.show()
