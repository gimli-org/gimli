#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import pygimli as pg
import pygimli.solver as solver
from pygimli.viewer import showMesh
from pygimli.mplviewer import drawMesh, drawStreams
from pygimli.meshtools import createMesh

import matplotlib.pyplot as plt
import numpy as np

from solverFVM import solveFiniteVolume, createFVPostProzessMesh

# build domain
nSteps = 20
dPhi = (0.6 * np.pi)/nSteps
boundary = []

for i in range(1, nSteps+1):
    boundary.append([np.cos(dPhi*i), np.sin(dPhi*i)])
    
poly = pg.Mesh(2)
nodes = []    
for b in boundary:
    nodes.append(poly.createNode(b))
for b in boundary[::-1]:
    nodes.append(poly.createNode(pg.RVector3(b)*0.1))

for i in range(len(nodes)):
    poly.createEdge(nodes[i], nodes[(i+1)%len(nodes)], 1)

mesh = createMesh(poly, quality=34, area=0.001, smooth=[0,10])
f = pg.RVector(mesh.cellCount(), 10)
a = pg.RVector(mesh.cellCount(), 0.1)

#Start FEM solution
swatch = pg.Stopwatch(True)
u = solver.solvePoisson(mesh, a=a, f=f,
                        uBoundary=[1, lambda p_: np.sin(np.arctan2(p_[1], p_[0]))/p_.abs()])

print('FEM:', swatch.duration(True))

#print(min(u), max(u))
ua = np.array(list(map(lambda p_: np.sin(np.arctan2(p_[1], p_[0]))/p_.abs(), mesh.positions())))
print(min(u), max(u))
ax1,cbar = showMesh(mesh, data=u, nLevs=20,
                    cMin=0, cMax=10, colorBar=True, showLater=True)
drawMesh(ax1, mesh)
#drawStreamLines2(ax1, mesh, data=u)
#ax2,cbar = showMesh(mesh, data=(u+1e-6)/(ua+1e-6), filled=True, colorBar=True, showLater=True)
#showMesh(amesh)


bounds = mesh.findBoundaryByMarker(1,99)

print('---:', swatch.duration(True))

uDirichlet = np.array(list(map(lambda p_: np.sin(np.arctan2(p_.center()[1], 
                                                            p_.center()[0]))/p_.center().abs(),
                               bounds)))
u = solveFiniteVolume(mesh, a=a, f=f, uDirichlet=uDirichlet)
print('FVM:', swatch.duration(True))

ax2,cbar = showMesh(mesh, 
                    data=u, 
                    cMin=0, cMax=10, logScale=False,
                    colorBar=True, showLater=True)
drawMesh(ax2, mesh)

mesh2, boundSortIdx = createFVPostProzessMesh(mesh, bounds)

print('---:', swatch.duration(True))
ax3,cbar = showMesh(mesh2, data=pg.cat(u, uDirichlet[boundSortIdx]),
                    cMin=0, cMax=10, colorBar=True, showLater=True)

drawMesh(ax3, mesh2)



#ax3,cbar = showMesh(mesh,
                    #data=np.array(list(map(lambda p_: np.sin(np.arctan2(p_[1],p_[0]))/p_.abs(), mesh.cellCenter()))),
                    #cMin=0, cMax=10, logScale=False,
                    #showLater=True)
#drawMesh(ax3, mesh)

plt.show()
#drawMesh(ax, grid)