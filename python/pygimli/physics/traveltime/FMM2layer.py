#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Fast Marching test using a two-layer model."""
import time
from math import asin, tan

import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.mplviewer import drawMesh, drawField, drawStreamLines
from pygimli.physics.traveltime import fastMarch

r"""
Solve the Hamilton-Jacobi (HJ) equation, known as the Eikonal equation
.. math::
    |\grad u(x)| & = f(x) \\
    ||\grad t||_2 &= s
    \grad t^2 &= s^2~\cite{SunFomel2000}

where :math:`t` denote traveltime for a spatial distributed slowness :math:`s`

In case f(x) = 1, the solution gives the distance from the boundary
"""


if __name__ == '__main__':
    xmin, xmax, zlay = -20., 150., 20.
    # create PLC (piece-wise linear complex) of two layers
    plc = pg.Mesh(2)
    nodes = []
    #   0-----------1
    #   |           |
    #   5-----------2
    #   |           |
    #   4-----------3
    nodes.append(plc.createNode(xmin, 0., 0.))
    nodes.append(plc.createNode(xmax, 0., 0.))
    nodes.append(plc.createNode(xmax, -zlay, 0.))
    nodes.append(plc.createNode(xmax, -zlay * 2, 0.))
    nodes.append(plc.createNode(xmin, -zlay * 2, 0.))
    nodes.append(plc.createNode(xmin, -zlay, 0.))
    # connect the nodes
    for i in range(5):
        plc.createEdge(nodes[i], nodes[i + 1])

    plc.createEdge(nodes[5], nodes[0])
    plc.createEdge(nodes[5], nodes[2])

    # insert region markers into the two layers and make mesh
    tri = pg.TriangleWrapper(plc)
    plc.addRegionMarker(pg.RVector3(0., -zlay + .1), 0, 3.)  # 10m^2 max area
    plc.addRegionMarker(pg.RVector3(0., -zlay - .1), 1, 10.)
    tri.setSwitches('-pzeAfaq34.6')
    mesh = pg.Mesh(2)
    tri.generate(mesh)
    mesh.createNeighbourInfos()
    print(mesh)

    # make velocity model
    v = [1000., 3000.]
    slomap = pg.stdMapF_F()  # map for mapping real slowness values
    for i, vi in enumerate(v):
        slomap.insert(i, 1. / vi)

    mesh.mapCellAttributes(slomap)  # map values to attributes using map

    # initialize source position and trvel time vector
    source = pg.RVector3(0., 0.)
    times = pg.RVector(mesh.nodeCount(), 0.)

    # initialize sets and tags
    upwind, downwind = set(), set()
    upTags, downTags = np.zeros(mesh.nodeCount()), np.zeros(mesh.nodeCount())

    # define initial condition
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

    # start fast marching
    tic = time.time()
    while len(downwind) > 0:
        fastMarch(mesh, downwind, times, upTags, downTags)

    print(time.time() - tic, "s")

    # compare with analytical solution along the x axis
    x = np.arange(0., 100., 0.5)
    t = pg.interpolate(mesh, times, x, x * 0., x * 0.)
    tdirect = x / v[0]  # direct wave
    alfa = asin(v[0] / v[1])  # critically refracted wave angle
    xreflec = tan(alfa) * zlay * 2.  # first critically refracted
    trefrac = (x - xreflec) / v[1] + xreflec * v[1] / v[0]**2
    tana = np.where(trefrac < tdirect, trefrac, tdirect)  # minimum of both
    print("min(dt)=",
          min(t - tana) * 1000,
          "ms max(dt)=",
          max(t - tana) * 1000,
          "ms")

    # plot traveltime field, a few lines
    fig = plt.figure()
    a = fig.add_subplot(211)
    drawMesh(a, mesh)
    drawField(a, mesh, times, True, 'Spectral')
    drawStreamLines(a, mesh, times, nx=50, ny=50)

    # plot calculated and measured travel times
    a2 = fig.add_subplot(212)
    plt.plot(x, t, 'b.-', x, tdirect, 'r-', x, trefrac, 'g-')
    plt.show()
