#!/ussr/bin/env python
# -*- coding: utf-8 -*-
"""
    
Modelling
---------

This is the first step for the modelling tutorial where we actually use finite elements computation. 

We will not go in deep detail about the Finite Elements Theory here, as this can be found in several books, e.g., 
:cite:`Zienkiewicz1977`

We just want to solve some problems to show how the *M* (Modelling) in *GIMLi* works.

We will solve a simple Poisson equation with zero boundary values, but a nonzero right hand side:

.. math::

    \nabla\cdot( A \cdot \nabla u ) + B u + C = 0
   
.. math::

    - \Delta u & = 1 \quad{\mathrm{in}}\quad\Omega\\
               u & = 0 \quad{\mathrm{on}}\quad\partial\Omega\\
  

Model domain is the unit square: :math:`\Omega=[-1, 1]^2`

As usually, the library needs to be imported first. 
Additionally we import the Poisson solver plus some plotting routines.
"""

import pygimli as pg
from pygimli.solver import solvePoisson
from pygimli.viewer import showMesh
from pygimli.mplviewer import drawMesh

""" 
We make use of NumPy to create vectors and the plotting library of matplotlib (Matlab-like plotting).

"""

import numpy as np
import matplotlib.pyplot as plt
    
"""
We create a grid with equidistant vectors and call the solver by providing  values for A (the material property) and the values at the four boundary with the markers 1 to 4.
"""

grid = pg.createGrid(x=np.linspace(-1.0, 1.0, 10), y=np.linspace(-1.0, 1.0, 10))

#material?

u = solvePoisson(grid, f=1.,
                 uBoundary=[grid.findBoundaryByMarker(1,5), 0.0],
                 verbose=True)

"""
.. error::

    can we calculate analytical solution?

The result is drawn with the function showMesh.    
"""

ax = showMesh(grid, data=u, filled=True, showLater=True, colorBar=True, orientation='vertical', label='P1 Solution $u$')[0]
drawMesh(ax, grid)

"""
.. image:: PLOT2RST.current_figure
    :scale: 50

We repeat the computation with a spatially (H) refined version of the original grid.
"""

gridh2 = grid.createH2()

uh = solvePoisson(gridh2, f=1.,
                  uBoundary=[gridh2.findBoundaryByMarker(1,5), 0.0],
                  verbose=True)

ax = showMesh(gridh2, data=uh, filled=True, showLater=True,colorBar=True, orientation='vertical', label='P2 Solution $u$')[0]

drawMesh(ax, gridh2)

"""
The same we do using a quadratic (P) refinement.
"""

gridp2 = grid.createP2()

up = solvePoisson(gridp2, f=1.,
                  uBoundary=[gridp2.findBoundaryByMarker(1,5), 0.0],
                  verbose=True)

"""
.. image:: PLOT2RST.current_figure
    :scale: 50

To analyse the solution in detail we use a numpy Nx2 matrix containing the x and y values, and the interpolation function.
"""

x = np.linspace(-1.0, 1.0, 100)

probe = np.zeros((len(x), 3))
probe[:, 0] = x

uH1 = pg.interpolate(mesh=grid, data=u, pos=probe)
uH2 = pg.interpolate(mesh=gridh2, data=uh, pos=probe)
uP2 = pg.interpolate(mesh=gridp2, data=up, pos=probe)

plt.figure()
plt.plot(x, uH1, label='linear (H1)')
plt.plot(x, uH2, label='linear (H2)')
plt.plot(x, uP2, label='quadratic (P2)')

plt.legend()


"""
.. image:: PLOT2RST.current_figure
    :scale: 50

"""

plt.show()
