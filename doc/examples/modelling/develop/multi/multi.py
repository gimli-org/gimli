#!/usr/bin/env python

"""
Test multi
"""

import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

import pygimli as pg
from pygimli.viewer import *
from pygimli.solver import *
from pygimli.meshtools import *

from pygimli.physics.seismics import ricker, wiggle, solvePressureWave
from pygimli.physics.gravimetry import solveGravimetry

def createTestWorld(maxArea=0.2, verbose=0):
    boundary = []
    boundary.append([-20.0,   0.0]) #0
    boundary.append([-20.0,  -2.0]) #1
    boundary.append([-20.0, -12.0]) #2
    boundary.append([-20.0, -20.0]) #3
    boundary.append([ 20.0, -20.0]) #4
    boundary.append([ 20.0, -12.0]) #5
    boundary.append([ 20.0,  -2.0]) #6
    boundary.append([ 20.0,   0.0]) #7
       
    poly = pg.Mesh(2)
    nodes = [poly.createNode(b) for b in boundary]

    boundaryMarker = [1, 2, 1, 3, 1, 4, 1, 1]
    for i in range(len(boundary)):
        poly.createEdge(nodes[i], nodes[(i+1)%len(boundary)], boundaryMarker[i]) # dirichlet (inflow)
    
    poly.createEdge(nodes[1], nodes[6], 5)
    poly.createEdge(nodes[2], nodes[5], 6)
        
    boundary = []
    boundary.append([-6.0, -4.0]) #0
    boundary.append([-6.0, -8.0]) #1
    boundary.append([ 6.0, -8.0]) #2
    boundary.append([ 6.0, -4.0]) #3
    
    nodes = [poly.createNode(b) for b in boundary]
    for i in range(len(boundary)):
        poly.createEdge(nodes[i], nodes[(i+1)%len(boundary)], 10) # dirichlet (inflow)
  
    poly.addRegionMarker((-19, -1), 1)
    poly.addRegionMarker((-19, -19), 1)
    poly.addRegionMarker((-19, -10), 2)
    poly.addRegionMarker((  0, -7), 3)
  
    mesh = createMesh(poly, quality=33, area=maxArea, smooth=[1,10], verbose=verbose)
    return mesh

def calcVelocity(mesh):
    velBoundary=[[1, [0.0, 0.0]],
                 [2, [1.0, 0.0]],
                 [4, [1.0, 0.0]],
                 [3, [0.0, 0.0]]
                ]
    preBoundary=[[1, 0.0],
                 [3, 0.0]]

    viscosity = pg.solver.parseMapToCellArray([[1, 2.2], [2, 1], [3, 10.2]], mesh)
    viscosity /= 100.

    solutionName = 'cache-vel-' + str(mesh.nodeCount())
    try:
        #vel = pg.load(solutionName + '.bmat')
        vel = np.load(solutionName + '.bmat.npy')
        
    except Exception as e:
        print(e)
        vel, pres, pCNorm, divVNorm = solveStokes_NEEDNAME(mesh, velBoundary,
                                                           preBoundary,
                                                           viscosity=viscosity,
                                                           maxIter=200,
                                                           verbose=1)
        np.save(solutionName + '.bmat', vel)
    return viscosity, vel

    
pg.showLater(1)

mesh = createTestWorld(maxArea=0.2, verbose=0)
meshC = createTestWorld(maxArea=1, verbose=0)

viscosity, vel = calcVelocity(mesh)

ax, cax = show(mesh, data=viscosity, colorBar=1, label='Viscosity')
show(mesh, axes=ax)

ax, cbar = show(mesh, data=pg.cellDataToPointData(mesh, np.sqrt(vel[:,0]*vel[:,0] + vel[:,1]*vel[:,1])),
                logScale=False, colorBar=True, label='|Velocity| in m/s', cbar='b2r')

pg.mplviewer.drawMeshBoundaries(ax, mesh, fitView=True)
show(mesh, data=vel, axes=ax, coarseMesh=meshC)

plt.pause(0.1)

f = pg.RVector(mesh.cellCount(), 0.0)
sourceCell=mesh.findCell([-18., -6])
f[sourceCell.id()]=0.4

times = np.linspace(0, 30, 30)
Peclet = 50.0
densityBrine = 1.2 # g/cm^3
densityBrine *= 1000 # kg/m^3

uMesh1 = solveFiniteVolume(mesh, a=1./Peclet, f=f, vel=vel, times=times, 
                          uBoundary=[2, 0],
                          scheme='PS', verbose=10)
uMesh2 = solveFiniteVolume(mesh, a=1./Peclet, f=0, vel=vel, times=times, 
                          uBoundary=[2, 0], u0=uMesh1[-1],
                          scheme='PS', verbose=10)

conc = np.vstack((uMesh1, uMesh2))

fig = plt.figure()
axCon = fig.add_subplot(2,1,2)
axGra = fig.add_subplot(2,3,2)
show(mesh, data=uMesh1[1], colorBar=1, label='Concentration', cMin=0, cMax=0.03, axes=axCon)
x, y, triangles, z, zIdx = pg.mplviewer.createTriangles(mesh, uMesh1[0])

dz = np.zeros(2*len(times))

for i in range(1, len(conc)):
    tic = time.time()
    axCon.clear()
    gci = axCon.tripcolor(x, y, triangles, conc[i])
    gci.set_clim(0, 0.03)
        
    axGra.clear()
    dDensity = densityBrine * conc[i]
    dg, dgz = solveGravimetry(mesh, dDensity, pnts=[[0.0, 0.0]]) 
    dz[i] = dg[0][2]
    axGra.plot(dz)
    axGra.set_ylabel('Grav at (0.0, 0.0) in mGal')
    axGra.set_xlabel('Time in s')
    
    print(i, time.time()-tic, "sum:", sum(conc[i]), "dsum:", (sum(conc[i])-sum(conc[i-1])),
          "sourceCell:", sourceCell.size())
    plt.pause(0.001)

pg.showNow()
