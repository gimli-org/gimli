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
import matplotlib.pyplot as plt
from matplotlib import animation

    
from pygimli.viewer import show
from pygimli.mplviewer import drawStreams

import time

grid = pg.createGrid(x=np.linspace(-1.0, 1.0, 31),
                     y=np.linspace(-1.0, 1.0, 31))

vals = pg.RVector(grid.cellCount(), 1.)
for c in grid.cells():
    if abs(c.center()[0]) < 0.1:
        vals[c.id()] = 10.0
        

#grid = grid.createP2()
times = np.arange(0, 2, 1./20)
#material?
            
neumannBC = [[1, -0.5], # left
             [2, 0.5]]  # right

dirichletBC = [[3, lambda b, t: 1.0 + np.cos(2.0 * np.pi * t)], # top
              [4, 1.0]] #bottom
               
pg.showLater(1)
ax = show(grid)[0]
ax.figure.canvas.draw()

#plt.ion()
#plt.show()

    
u = solve(grid, a=vals, f=0.5, 
          times=times,
          u0=pg.RVector(grid.nodeCount(), 0.0),
          duBoundary=neumannBC,
          uBoundary=dirichletBC,
          ##plotTimeStep=updateDrawU,
          verbose=False, progress=True)

uMin = min(u.flat)
uMax = max(u.flat)

show(grid, u[0], axes=ax)

"""
.. image:: PLOT2RST.current_figure
    :scale: 75

"""

def gen():
    mesh = gen.mesh
    u = gen.u
    for i, ui in enumerate(u[1:]):
        yield i, len(u), ui, mesh
        
gen.mesh = grid
gen.u = u
    
def animate(data):
    i = data[0]
    imax = data[1]
    ui = data[2]
    mesh = data[3]
    print(i,'/', imax)    
    
    global ax
    
    ax.clear()
    ax = show(mesh, data=ui, showLater=True, axes=ax,
              levels=np.linspace(0, 3, 16))[0]

    if min(ui) != max(ui):
        pass
        #drawStreams(ax, mesh, ui)

anim = animation.FuncAnimation(plt.gcf(), animate,
                               gen,
                               interval=2)

"""
.. animate:: anim fps=10 bitrate=1024 dpi=92
    
"""

#plt.show()

pg.showNow()
#plt.show()
