#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    
Equation of heat
----------------

Let us start with a mathematical formulation ...

.. math::

    \nabla\cdot( a \cdot \nabla u ) + \frac{\partial u}{\partial t} + C = 0
   
.. math::

    - \Delta u & = 1 \quad{\mathrm{in}}\quad\Omega\\
               u & = 0 \quad{\mathrm{on}}\quad\partial\Omega\\
  
We will solve this equation on the unit square: :math:`\Omega=[-1, 1]^2`
"""

import pygimli as pg
from pygimli.solver import solve

"""
"""
import numpy as np
    
grid = pg.createGrid(x=np.linspace(-1.0, 1.0, 21),
                     y=np.linspace(-1.0, 1.0, 21))

vals = pg.Vector(grid.cellCount(), 1.)
for c in grid.cells():
    if abs(c.center()[0]) < 0.1:
        vals[c.id()] = 1.0
        
#grid = grid.createP2()
times = np.arange(0, 1, 1./50)
#material?
            
bcN = {1: -0.5, # left
       2:  0.5, # right
       }  

def bcD4(p, time=0.):
    return 15. * np.cos(2.0 * np.pi * time)

bcD = {4: bcD4, # top
       } 
               
ax,_ = pg.show(grid, vals)
    
u = solve(grid, a=vals, f=0.5, 
          times=times,
          u0=pg.Vector(grid.nodeCount(), 0.0),
          bc={'Dirichlet':bcD, 'Neumann':bcN},
          dynamic=True, verbose=False, progress=True,
          theta=0.6)

q = [pg.solver.grad(grid, u[i]) for i in range(len(u))]
pg.show(grid, u, flux=q, ax=ax, cMin=-15, cMax=15)
