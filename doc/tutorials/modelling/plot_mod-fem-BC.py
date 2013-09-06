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
               u & = 0 \quad{\mathrm{on}}\quad\partial\Omega\\
  

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

grid = grid.createP2()

for b in grid.boundaries():
    # default grid boundary marker is 1
    if b.marker() == 1:
        #if x-component of outward norm direction shows right set marker to 2
        if b.norm()[0] == 1.0:
            b.setMarker(3)
        #if y-component of outward norm direction shows up set marker to 3
        if b.norm()[1] == 1.0:
            b.setMarker(4)
        #if y-component of outward norm direction shows down set marker to 4
        if b.norm()[1] == -1.0:
            b.setMarker(2)

"""
Define a list of pairs for boundaries and potential values for the Dirichlet boundaries. 
Definition of the used boundaries are either by the marker directly or by a given list of boundaries.
The value can be defined directly or by a given function.

"""

def uDirichlet(boundary):
    '''
        Return a solution value for the nodes of the given boundary.
    '''
    return 4.0

dirichletBC = [[1, 1.0], # left
               [grid.findBoundaryByMarker(2), 2.0], # right
               [grid.findBoundaryByMarker(3), 3.0], # upper
               [grid.findBoundaryByMarker(4), uDirichlet]] # lower

"""

"""
            
u = solvePoisson(grid, f=1.,
                 uBoundary=dirichletBC,
                 verbose=True)

ax = showMesh(grid, data=u, filled=True, colorBar=True,
              orientation='vertical', label='Solution $u$',
              levels=np.linspace(1.0, 4.0, 14), showLater=True)[0]

drawMesh(ax, grid)

drawSelectedMeshBoundaries(ax, grid.findBoundaryByMarker(1), 
                           color=(1.0, 0.0, 0.0), linewidth=3)
ax.text(-1.08, 0.0, '1', color=(1.0, 0.0, 0.0))
drawSelectedMeshBoundaries(ax, grid.findBoundaryByMarker(2),
                           color=(0.0, 0.66, 0.0), linewidth=3)
ax.text( 0.0, -1.08, '2', color=(0.0, 0.66, 0.0))
drawSelectedMeshBoundaries(ax, grid.findBoundaryByMarker(3),
                           color=(0.0, 0.0, 1.0), linewidth=3)
ax.text( 1.02, 0.0, '3', color=(0.0, 0.0, 1.0))
drawSelectedMeshBoundaries(ax, grid.findBoundaryByMarker(4),
                           color=(0.7, 0.7, 0.0), linewidth=3)
ax.text( 0.0, 1.02, '4', color=(0.7, 0.7, 0.0))

ax.set_xlim([-1.1, 1.1])
ax.set_ylim([-1.1, 1.1])
   
   
"""
.. image:: PLOT2RST.current_figure

"""

P.show()
