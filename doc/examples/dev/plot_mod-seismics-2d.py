# coding=utf-8

"""
Simple 2 layer full waveform pressure wave.
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

def drawWaveField(ax, u):
    ui = u / max(pg.abs(u))
    ui = pg.logDropTol(ui, 1e-2)
    cMax = max(pg.abs(ui))

    drawField(ax, mesh, data=ui, cMin=-cMax, cMax=cMax, cmap='RdBu',
              #interpolate=0, shading='gouraud'
              )

def createMesh2():
    boundary = []
    boundary.append([-20.0,   0.0])
    boundary.append([-20.0, -20.0])
    boundary.append([ 20.0, -20.0])
    boundary.append([ 20.0,   0.0])

    poly = pg.Mesh(2)
    nodes = [poly.createNode(b) for b in boundary]

    poly.createEdge(nodes[0], nodes[1], 1) # dirichlet (inflow)
    poly.createEdge(nodes[1], nodes[2], 3) # hom neumann (outflow)
    poly.createEdge(nodes[2], nodes[3], 2) # hom dirichlet (isolation)
    poly.createEdge(nodes[3], nodes[0], 4) # hom dirichlet (isolation)

    mesh = createMesh(poly, quality=34, area=0.05, smooth=[0,10])
    return mesh

dx = 0.2
x = np.arange(-20, 20., dx)
y = np.arange(-20, 0.0, dx)[::-1]

mesh = pg.createGrid(x=x, y=y)
mesh = createMesh2()
print(mesh)

h = pg.math.median(mesh.boundarySizes())

v1 = 1000
v2 = 3000
tmax = 10.1/v1

z = 2.
f0 = 1000.0 # A low wavelength of 50 Hz

velocities = pg.Vector(mesh.cellCount(), v1)

for c in mesh.cells():
    velocities[c.id()] = v1
    if c.center()[1] < -z:
        velocities[c.id()] = v2

dt = h * 0.5/max(velocities)
times = np.arange(0.0, tmax, dt)
print(mesh, "h:", h, "dt:", dt, "n_t", len(times))

uSource = ricker(f0, times, t0=1./f0)

solutionName = 'uGridBig-'+ str(dt) + '-' + str(h)

nx = len(x)
ny = len(y)
dt = times[1] - times[0]

try:
    print("load:", solutionName)
    u = pg.load(solutionName + '.bmat')
except Exception as e:
    print(e)
    u = solvePressureWave(mesh, velocities, times, sourcePos=(0.0, 0.0),
                          uSource=uSource, verbose=10)
    u.save(solutionName + '.bmat')

print(u)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
Ustep = len(u)/2
field = drawWaveField(ax, u[Ustep, :])

plt.show()

def animate(i):
    if i > 0:
        tic = time.time()
        i = i * 10
        ax.clear()
        uShow = drawWaveField(ax, u[i,:])
        print(i, time.time()-tic, len(u), min(pg.abs(u[i,:])), min(u[i,:]), max(u[i,:]), dt*i)


anim = animation.FuncAnimation(fig, animate,
                               frames=int(len(u)),
                               interval=20)

anim.save(solutionName + ".mp4", writer=None, fps=20, dpi=92, codec=None,
          bitrate=24*320, extra_args=None, metadata=None,
          extra_anim=None, savefig_kwargs=None)


#fig.canvas.draw()
plt.show()
