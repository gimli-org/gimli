#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fast Marching Method test using a two-layer model
-------------------------------------------------

This example shows how the FMM implementation of pyGIMLi works. FMM is an
alternative to the Shortest path method and utilized by the Refraction manager.
We illustrate the effectiveness and compare the solutions."""

###############################################################################
# The first-arrival travel-time :math:`t` from a shot point at the origin is
# described by a Hamilton-Jacobi (HJ) type equation, known as Eikonal equation
# :cite:`podvin1991GJI`
#
# .. math::
#
#    (\nabla t)^2 = s^2
#
# where :math:`s` is the spatially distributed slowness (i.e. 1 / velocity).
# Under the assumption of piece-wise slowness (for individual cells) one can
# solve this equation by the so-called Fast-Marching method as first developed
# for 2D grids by Podvon & Lecompte (1991). Starting from the source location,
# a front is successively increased and filled by the travel times.
#
# This example demonstrates the core features of an FMM Python implementation
# and compares the travel time for a two-layer model with the analytical
# solution (known for straight or sloped layers from the text books) and the
# shortest-path algorithm (Moser, 1991) implemented in the pyGIMLi core.
# :cite:`moser1991G`


import time
from math import asin, tan

import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.mplviewer import drawMesh, drawField, drawStreamLines
from pygimli.physics.traveltime import fastMarch


###############################################################################
# First we provide the analytical solution for a given offset vector x.
def analyticalSolution2Layer(x, zlay, v1, v2):
    """Analytical solution: minimum of direct and critically refracted wave."""
    tdirect = np.abs(x) / v1  # direct wave
    alfa = asin(v1 / v2)  # critically refracted wave angle
    xreflec = tan(alfa) * zlay * 2.  # first critically refracted
    trefrac = (x - xreflec) / v2 + xreflec * v2 / v1**2
    return np.minimum(tdirect, trefrac)

###############################################################################
# The model consists of two boxes, in the first is the source (1)
#
# .. code-block:: none
#
#     0--1--------2
#     |           |
#     6-----------3
#     |           |
#     5-----------4
#
# We create a PLC (piece-wise linear complex) and insert the nodes.
xmin, xmax, zlay = -20., 150., 20.  # model dimensions
plc = pg.Mesh(2)
nodes = []
nodes.append(plc.createNode(xmin, 0., 0.))  # 0
nodes.append(plc.createNode(0.0, 0., 0.))  # 1
nodes.append(plc.createNode(xmax, 0., 0.))  # 2
nodes.append(plc.createNode(xmax, -zlay, 0.))  # 3
nodes.append(plc.createNode(xmax, -zlay * 2, 0.))  # 4
nodes.append(plc.createNode(xmin, -zlay * 2, 0.))  # 5
nodes.append(plc.createNode(xmin, -zlay, 0.))  # 6

###############################################################################
# The nodes are connected from from 0 to 6 and back to 0.
# An additional edge is drawn from 6 to 3. Node/edge markers do not matter.
for i in range(6):
    plc.createEdge(nodes[i], nodes[i + 1])

plc.createEdge(nodes[6], nodes[0])
plc.createEdge(nodes[6], nodes[3])

###############################################################################
# We insert region markers (0 and 1) into the two layers and generate the mesh.
tri = pg.TriangleWrapper(plc)
plc.addRegionMarker(pg.RVector3(0., -zlay + .1), 0, 3.)  # 10m^2 max area
plc.addRegionMarker(pg.RVector3(0., -zlay - .1), 1, 10.)
tri.setSwitches('-pzeAfaq34.6')
mesh = pg.Mesh(2)
tri.generate(mesh)
mesh.createNeighbourInfos()
print(mesh)

###############################################################################
# Next we generate a velocity model from the markers by using a map.
# The values are associated to the markers and stored as attributes.
v = [1000., 3000.]
slomap = pg.stdMapF_F()  # mapping markers to real slowness values
for i, vi in enumerate(v):
    slomap.insert(i, 1. / vi)

mesh.mapCellAttributes(slomap)  # map values to attributes using map

###############################################################################
# We initialize the source position and the travel time vector
# initialize sets and tags and define the initial condition.
source = pg.RVector3(0., 0.)  # does not have to be a node!
times = pg.RVector(mesh.nodeCount(), 0.)
upwind, downwind = set(), set()
upTags, downTags = np.zeros(mesh.nodeCount()), np.zeros(mesh.nodeCount())
cell = mesh.findCell(source)
for i, n in enumerate(cell.nodes()):
    times[n.id()] = cell.attribute() * n.pos().distance(source)
    upTags[n.id()] = 1
for i, n in enumerate(cell.nodes()):
    tmpNodes = pg.commonNodes(n.cellSet())
    for nn in tmpNodes:
        if not upTags[nn.id()] and not downTags[nn.id()]:
            downwind.add(nn)
            downTags[nn.id()] = 1

###############################################################################
# Then we start marching until all fields are filled.
tic = time.time()
while len(downwind) > 0:
    fastMarch(mesh, downwind, times, upTags, downTags)

print(time.time() - tic, "s")

###############################################################################
# First, we plot the traveltime field and streamlines
fig, ax = plt.subplots(figsize=(10, 5))
drawMesh(ax, mesh)
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
drawField(ax, mesh, times, cMap='Spectral', fillContour=True)
drawStreamLines(ax, mesh, -times, nx=50, ny=50)

###############################################################################
# We compare the result with the analytical solution along the x axis
x = np.arange(0., 140., 0.5)
tFMM = pg.interpolate(mesh, times, x, x * 0., x * 0.)
tAna = analyticalSolution2Layer(x, zlay, v[0], v[1])
print("min(dt)={} ms  max(dt)={} ms".format(min(tFMM - tAna) * 1000,
                                            max(tFMM - tAna) * 1000))

###############################################################################
# In order to use the Dijkstra, we extract the surface positions >0
mx = pg.x(mesh.positions())
my = pg.y(mesh.positions())
fi = pg.find((my == 0.0) & (mx >= 0))
px = np.sort(mx(fi))

###############################################################################
# A data container with index arrays named s (shot) and g (geophones) is
# created and filled with the positions and shot/geophone indices.
data = pg.DataContainer()
data.registerSensorIndex('s')
data.registerSensorIndex('g')
for pxi in px:
    data.createSensor(pg.RVector3(pxi, 0.0))

ndata = len(px) - 1
data.resize(ndata)
data.set('s', pg.RVector(ndata, 0))  # only one shot at first sensor
data.set('g', pg.utils.grange(1, ndata, 1))  # all others and geophones
fop = pg.TravelTimeDijkstraModelling(mesh, data)
tDijkstra = fop.response(mesh.cellAttributes())

###############################################################################
# We plot the calculated and measured travel times and relative differences
fig, ax = plt.subplots()
ax.plot(x, tAna*1000, 'r-', label='analytical')
ax.plot(x, tFMM*1000, 'b+-', label='FMM')
ax.plot(px[1:], tDijkstra*1000, 'gx-', label='Dijkstra')
ax.set_xlabel('x [m]')
ax.set_ylabel('t [s]')
ax.grid(True)
ax.legend()
ax2 = ax.twinx()
dtFMM = (tFMM - tAna) * 1000
tAnaD = analyticalSolution2Layer(px[1:], zlay, v[0], v[1])
dtDijkstra = (tDijkstra - tAnaD) * 1000
ax2.set_ylabel('dt [ms]')
if 1:  # relative differences in percent
    ax2.set_ylabel('dt [%]')
    dtFMM /= (tAna * 10)
    dtDijkstra /= (tAnaD * 10)
ax2.plot(x, dtFMM, 'b.-')
ax2.plot(px[1:], dtDijkstra, 'g.-')
pg.wait()

###############################################################################
# Note that the Fast Marching Method is implemented in a modelling class
# (right now Python and very slow but to be replaced by fast C++ soon)
# FMModelling that can be more easily used with the Refraction manager by
# ra.useFMM(True)
