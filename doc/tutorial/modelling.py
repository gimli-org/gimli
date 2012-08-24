#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''

*This introductory sentence should state the intent and goal of the tutorial. Keep it brief.
This next block should state any assumptions that you the writer are making. Present them in list form.*

Let us start with ...

.. math::

    \nabla\cdot( A \cdot \nabla u ) + B u + C = 0
    
This is the first tutorial example where we actually use finite elements to compute something. 
We will solve a simple version of Poisson's equation with zero boundary values, but a nonzero right hand side:

.. math::

    - \Delta u & = 1 \quad{\mathrm{in}}\quad\Omega\\
               u & = 0 \quad{\mathrm{on}}\quad\partial\Omega\\
  

\begin{align*} -\Delta u &= 1 \qquad\qquad & \text{in}\ \Omega, \\ u &= 0 \qquad\qquad & \text{on}\ \partial\Omega. \end{align*}

We will solve this equation on the unit square, $\Omega=[0,1]^2$

First, the library must be imported.
To avoid name clashes with other libraries we suggest to import `pygimli` and alias it to a simple abbreviation 'g', e.g., by using
'''

import pygimli as g
from pygimli.solver import solvePoisson

'''
As a result, all :ref:`gimli` objects (classes and functions) can be referred to with a preceding `g.`, e.g., 
printing the version string for gimli.

'''
import numpy as np
    
grid = g.createGrid( np.linspace( 0, 1, 5 ), np.linspace( 0, 1, 5 ) )
#grid.createNeighbourInfos( True )
print grid.boundary( 1 ).leftCell()
print grid.boundary( 1 ).rightCell()

c = grid.cell( 0 )
for n in c.nodes():
    print n.id(), n.pos()
    
print c.shape().isInside( [1. / 3., 1. / 3.], True )
grid.findCell( [1. / 3., 1. / 3.] )
#.pot( [1. / 3., 1. / 3.], u)

for i in range( 4 ):
    u = solvePoisson( grid, f = 1., u0 = grid.findBoundaryByMarker( 1 ), verbose = True )
    
    #print grid.findCell( [1. / 3., 1. / 3.] ).pot( [1. / 3., 1. / 3.], u)
    grid = grid.createH2()
    #grid.createNeighbourInfos( True )

from pygimli.viewer import showMesh
from pygimli.mplviewer import drawMesh

print sum(u)

#grid = g.createGrid( np.linspace( 0, 1, 25 ), np.linspace( 0, 1, 25 ) )
#u = solvePoisson( grid, f = 1., u0 = grid.findBoundaryByMarker( 1 ), verbose = True )

ax = showMesh( grid, data = u, filled = True, showLater = True )
drawMesh( ax, grid )




'''
do we have an analytical solution for this example?
'''


import pylab as P
P.show()
