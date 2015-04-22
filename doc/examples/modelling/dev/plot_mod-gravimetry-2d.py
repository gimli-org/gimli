#!/usr/bin/env python

"""
Simple gravimetry field solution
"""

import sys
import time

import matplotlib.pyplot as plt
import numpy as np

import pygimli as pg
from pygimli.viewer import *
from pygimli.solver import *
from pygimli.meshtools import *
from pygimli.polytools import *

from pygimli.physics.gravimetry import *

world = createWorldPolygon(start=[-50, -50], end=[50, 50], area=0.1, marker=1)
circ = createCirclePolygon([0, -5], radius=3, area=0.1, marker=2)

mesh = createMesh([world, circ], quality=34)
mesh = mesh.createP2()
print(mesh)
density=solver.parseMapToCellArray([[1, 0.0], [2, 1000.0]], mesh)

u = solve(mesh, a=1, f=density, uB=[[-2,0], [-1,0]]) * pg.physics.constants.G

pnts = [pg.RVector3(x,0) for x in np.arange(-19, 19.1, 1)]

#dgz, dggz = solveGravimetry(mesh, dDensity=density,
                            #pnts=pnts, complete=1)

#print(min(dgz), max(dgz))

duz = np.zeros(len(pnts))
dux = np.zeros(len(pnts))
for i,p in enumerate(pnts):
    c = mesh.findCell(p)
    g = c.grad(p, u)
    print(c, p, g)
    dux[i] =  g[0] * 4.0 * np.pi # wo kommen die 4 pi her?
    duz[i] = -g[1] * 4.0 * np.pi # wo kommen die 4 pi her?

uI = pg.interpolate(mesh, u, pnts)

plt.plot(pg.x(pnts), dux, label='dux')
plt.plot(pg.x(pnts), duz, label='duz')
plt.plot(pg.x(pnts), pg.sqrt(dux*dux + duz*duz), label='du')

dgz, dggz = solveGravimetry(circ, dDensity=1000,
                            pnts=pnts, complete=1)

    
plt.plot(pg.x(pnts), dgz[:,0], label='dgx')
plt.plot(pg.x(pnts), dgz[:,2], label='dgz')
plt.plot(pg.x(pnts), np.sqrt(dgz[:,0]**2 + dgz[:,1]**2 + dgz[:,2]**2), label='dg')


print(dux/dgz[:,0])
print(duz/dgz[:,2])

#plt.plot(pg.x(pnts), uI, label='ui')

plt.legend()
ax, cbar= pg.show(mesh, data=u, showLater=1)
pg.mplviewer.drawMesh(ax, mesh)

plt.show()

