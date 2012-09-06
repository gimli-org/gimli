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
  

We will solve this equation on the unit square, $\Omega=[-1,1]^2$

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
    
#grid = g.createGrid( np.linspace( -1, 1, 5 ), np.linspace( -1, 1, 5 ) )

##.pot( [1. / 3., 1. / 3.], u)

#u = []

#for i in range( 4 ):
    #grid = grid.createH2()
    #grid.createNeighbourInfos()
    #u = solvePoisson( grid, f = 1., u0 = grid.findBoundaryByMarker( 1 ), verbose = True )
    
    #print grid.findCell( [1. / 3., 1. / 3.] ).pot( [1. / 3., 1. / 3.], u)

    ##grid.createNeighbourInfos( True )

from pygimli.viewer import showMesh
from pygimli.mplviewer import *
import pylab as P

grid = g.createGrid( np.linspace( -1, 1, 20 ), np.linspace( -1, 1, 20 ), np.linspace( -1, 1, 20 ) )
u = solvePoisson( grid, f = 1., uBoundary = [ grid.findBoundaryByMarker( 1 ), 0], verbose = True )

#ax = showMesh( grid, data = u, filled = True, showLater = True, colorBar = True, orientation = 'vertical', label = 'Solution er$u$' )[0]
#drawMesh( ax, grid )

'''
.. error::
    
    do we find an analytical solution for this example?
'''

for b in grid.boundaries():
    if b.marker() == 1:
        if b.norm()[0] == 1:
            b.setMarker( 2 )
        if b.norm()[1] == 1:
            b.setMarker( 3 )
        if b.norm()[1] == -1:
            b.setMarker( 4 )
            
u = solvePoisson( grid, f = 0., uBoundary = [ [ grid.findBoundaryByMarker( 2 ), 1.0 ], 
                                              [ grid.findBoundaryByMarker( 4 ), 0.0 ] ],  verbose = True )

#ax = showMesh( grid, data = u, filled = True, showLater = True, colorBar = True, orientation = 'vertical', label = 'Solution $u$'  )[0]
#drawMesh( ax, grid )
#drawSelectedMeshBoundaries( ax, grid.findBoundaryByMarker( 1 ), color = (1.0, 0.0, 0.0), linewidth=2 )
#drawSelectedMeshBoundaries( ax, grid.findBoundaryByMarker( 2 ), color = (0.0, 1.0, 0.0), linewidth=2 )
#drawSelectedMeshBoundaries( ax, grid.findBoundaryByMarker( 3 ), color = (0.0, 0.0, 1.0), linewidth=2 )
    
grid.addExportData('u', u )
grid.exportVTK( "grid" )
    
P.show()
