#!/ussr/bin/env python
# -*- coding: utf-8 -*-
"""
    
Modelling
---------

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
"""

import pygimli as g

from pygimli.solver import solvePoisson

'''

'''
import matplotlib.pyplot as plt
import numpy as np
    
from pygimli.viewer import showMesh
from pygimli.mplviewer import *


grid = g.createGrid(x=np.linspace(-1.0, 1.0, 10), y=np.linspace(-1.0, 1.0, 10))

#material?

u = solvePoisson(grid, f=1.,
                 uBoundary=[grid.findBoundaryByMarker(1,5), 0.0],
                 verbose=True)

"""
.. error::

    can we calculate analytical solution?
    
"""

ax = showMesh(grid, data=u, filled=True, showLater=True,
              colorBar=True, orientation='vertical', label='P1 Solution $u$')[0]
drawMesh(ax, grid)

"""
.. image:: PLOT2RST.current_figure
    :scale: 50

"""

gridp2 = grid.createP2()

up = solvePoisson(gridp2, f=1.,
                  uBoundary=[gridp2.findBoundaryByMarker(1,5), 0.0],
                  verbose=True)

ax = showMesh(gridp2, data=up, filled=True, showLater=True,
              colorBar=True, orientation='vertical', label='P2 Solution $u$')[0]

drawMesh(ax, gridp2)


"""
.. image:: PLOT2RST.current_figure
    :scale: 50
    
"""

xp = np.linspace(-1.0, 1.0, 100)

probe = np.zeros((len(xp), 3))
probe[:, 0] = xp
probe[:, 1] = xp

uProbe1 = g.interpolate(mesh=grid, data=u, pos=probe)
uProbe2 = g.interpolate(mesh=gridp2, data=up, pos=probe)

plt.figure()
plt.plot(xp, uProbe1, label='linear')
plt.plot(xp, uProbe2, label='quadratic')

plt.legend()


"""
.. image:: PLOT2RST.current_figure
    :scale: 50

"""

plt.show()
