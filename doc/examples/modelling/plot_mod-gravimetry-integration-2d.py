#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
Gravimetry in 2d
----------------

Simple gravimetric field solution with Integration after :cite:`WonBevis1987`.


"""

import numpy as np
import pygimli as pg
from pygimli.meshtools import createCircle
from pygimli.physics.gravimetry import solveGravimetry
from pygimli.physics.gravimetry import gradUCylinderHoriz, gradGZCylinderHoriz
from pygimli.physics.gravimetry import gradUHalfPlateHoriz, gradGZHalfPlateHoriz

radius  = 2.
depth   = 5.
rho     = 1000.0

x = np.arange(-20, 30, 1)
pnts = np.zeros((len(x),2))
pnts[:, 0] = x
pos=[0, -depth]

def plot(x, a1, ga, gza, a2, g, gz):
    a1.plot(x, ga[:,0],  label=r'Analytical $\frac{\partial u}{\partial x}$')
    a1.plot(x, ga[:,1],  label=r'Analytical $\frac{\partial u}{\partial z}$')

    a1.plot(x, g[:,0], label=r'Won & Bevis: $\frac{\partial u}{\partial x}$',
            marker='o', linewidth=0)
    a1.plot(x, g[:,2], label=r'Won & Bevis: $\frac{\partial u}{\partial z}$',
            marker='o', linewidth=0)

    a2.plot(x, gza[:,0],  label=r'Analytical $\frac{\partial^2 u}{\partial z,x}$')
    a2.plot(x, gza[:,1],  label=r'Analytical $\frac{\partial^2 u}{\partial z,z}$')

    a2.plot(x, gz[:,0], label=r'Won & Bevis: $\frac{\partial^2 u}{\partial z,x}$',
            marker='o', linestyle='')
    a2.plot(x, gz[:,2], label=r'Won & Bevis: $\frac{\partial^2 u}{\partial z,z}$',
            marker='o', linestyle='')
    a1.set_xlabel('$x$-coordinate [m]')
    a1.set_ylabel(r'$\frac{\partial u}{\partial (x,z)}$ [mGal]')
    a1.legend()
    
    a2.set_xlabel('$x$-coordinate [m]')
    a2.legend()


fig = pg.plt.figure()
ax = [fig.add_subplot(2,2,i) for i in range(1,5)]

# Horizontal cylinder

ga = gradUCylinderHoriz(pnts, radius, rho, pos=pos)
gza = gradGZCylinderHoriz(pnts, radius, rho, pos=pos)

circ = createCircle([0, -depth], radius=radius, marker=2, area=0.1, segments=32)
g, gz = solveGravimetry(circ, rho, pnts, complete=True)

plot(x, ax[0], ga, gza, ax[1], g, gz)

# Half plate

thickness = 0.1

###mesh = pg.createGrid(x=[-2,2], y=[-2,2], z=[-3,-7])
mesh = pg.createGrid(x=np.linspace(0, 5000, 2),
                     y=[-depth-thickness/2.0, -depth+thickness/2.0] )

ga = gradUHalfPlateHoriz(pnts, thickness, rho, pos=[0, -depth])
gza = gradGZHalfPlateHoriz(pnts, thickness, rho, pos=[0, -depth])
g, gz = solveGravimetry(mesh, rho, pnts, complete=True)

plot(x, ax[2], ga, gza, ax[3], g, gz)

pg.wait()
