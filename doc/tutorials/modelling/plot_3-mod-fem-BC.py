#!/ussr/bin/env python
# -*- coding: utf-8 -*-
"""
    
Modelling with Boundary Conditions
----------------------------------

We use the preceding example (Poisson equation on the unit square) but want to 
specify different boundary conditions on the four sides.

Again, we first import numpy and pygimli, the solver and post processing 
functionality.
"""

import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.solver import solve
from pygimli.viewer import show
from pygimli.mplviewer import drawStreams

"""
We create a 50x50 node grid to solve on.
"""

grid = pg.createGrid(x=np.linspace(-1.0, 1.0, 21),
                     y=np.linspace(-1.0, 1.0, 21))

"""
We start considering inhomogeneous Dirchlet boundary conditions (BC).
There are different ways of specifying BCs. 
They can be maps from markers to values, explicit functions or implicit (lambda) functions as exemplified in the next example.

The boundary 1 (top) and 2 (left) are directly mapped to the values 1.0 and 2.0.
On side 3 (bottom) a lambda function 3+x is used (p is the boundary position and p[0] its x coordinate.
On side 4 (right) a function uDirichlet is used that simply returns 4.0 in this example but can compute everything as a function of the individual boundaries b.
"""

def uDirichlet(b):
    '''
        Return a solution value for coordinate p.
    '''
    return 4.0

dirichletBC = [[1, 1.0], # left
               [grid.findBoundaryByMarker(2), 2.0], # right
               [grid.findBoundaryByMarker(3), lambda p: 3.0 + p[0]], # top
               [grid.findBoundaryByMarker(4), uDirichlet]] # bottom

"""
The BC are passed using the uBoundary keyword.
Note that showMesh returns the created figure axes ax while drawMesh plots on it and it can also be used as a class with plotting or decoration methods.
"""
u = solve(grid, f=1.,
          uBoundary=dirichletBC)

ax = show(grid, data=u, filled=True, colorBar=True,
          orientation='vertical', label='Solution $u$',
          levels=np.linspace(1.0, 4.0, 17), showLater=True)[0]

show(grid, axes=ax)

ax.text( 0.0, 1.02, '$u=1$')
ax.text(-1.08, 0.0, '$u=2$', rotation='vertical')
ax.text( 0.0, -1.08, '$u=3+x$')
ax.text( 1.02, 0.0, '$u=4$', rotation='vertical')

ax.set_title('$\\nabla\cdot(1\\nabla u)=1$')

ax.set_xlim([-1.1, 1.1]) # some boundary for the text
ax.set_ylim([-1.1, 1.1])

"""
.. image:: PLOT2RST.current_figure
    :scale: 75


Alternatively we can define the gradients of the solution on the boundary, i.e., Neumann type BC.
This is done as another map (marker, right-hand-side value) and passed by the keyword duBoundary.
"""
neumannBC = [[1, -0.5], # left
             [grid.findBoundaryByMarker(4), 2.5]] # bottom

dirichletBC = [3, 1.0] # top

u = solve(grid, f=0.,
          duBoundary=neumannBC,
          uBoundary=dirichletBC)

"""
Note that on boundary 4 (right) no BC is explicitly applied leading to default or natural BC that are of homogeneous Neumann type 
:math:`\frac{\partial u}{\partial n}=0`
"""
              
ax = show(grid, data=u, filled=True, colorBar=True,
          orientation='vertical', label='Solution $u$',
          levels=np.linspace(min(u), max(u), 14), showLater=True)[0]

"""
Instead of the grid we now want to add streamlines to the plot to show the gradients of the solution (i.e., the flow direction).
"""

drawStreams(ax, grid, u)

ax.text(0.0, 1.02, '$u=1$',
        horizontalalignment='center' ) #top -- 3
ax.text(-1.08, 0.0, '$\partial u/\partial n=-0.5$',
        verticalalignment='center', rotation='vertical') #left -- 1
ax.text(0.0, -1.08, '$\partial u/\partial n=2.5$',
        horizontalalignment='center') # bot -- 4
ax.text(1.02, 0.0, '$\partial u/\partial n=0$',
        verticalalignment='center', rotation='vertical') #right -- 2


ax.set_title('$\\nabla\cdot(1\\nabla u)=0$')

ax.set_xlim([-1.1, 1.1])
ax.set_ylim([-1.1, 1.1])

"""
.. image:: PLOT2RST.current_figure
    :scale: 75

"""

plt.show()
