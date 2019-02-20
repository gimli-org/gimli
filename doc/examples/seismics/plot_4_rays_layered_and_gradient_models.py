#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Raypaths in layered and gradient models
=======================================

This example performs raytracing for a two-layer and a vertical gradient model
and compares the resulting traveltimes to existing analytical solutions. An
approximation of the raypath is found by finding the shortest-path through a
grid of nodes. The possible angular coverage is small when only corner points of
a cell (primary nodes) are used for this purpose. The angular coverage, and
hence the numerical accuracy of traveltime calculations, can be significantly
improved by a few secondary nodes along the cell edges. Details can be found in
`Giroux & Larouche (2013) <https://doi.org/10.1016/j.cageo.2012.12.005>`_.
"""
# sphinx_gallery_thumbnail_number = 3
from math import asin, tan

import matplotlib.pyplot as plt
import numpy as np
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.mplviewer import drawMesh

################################################################################
# Two-layer model
# ---------------
# We start by building a regular grid.

mesh_layered = mt.createGrid(
    np.arange(-20, 155, step=5), np.linspace(-60, 0, 13))

################################################################################
# We now construct the velocity vector for the two-layer case by iterating over
# the cells. Cells above 25 m depth are assigned :math:`v = 1000` m/s and cells
# below are assigned :math:`v = 3000` m/s.

vel_layered = np.zeros(mesh_layered.cellCount())
for cell in mesh_layered.cells():
    if cell.center().y() < -25:
        vel = 3000.0
    else:
        vel = 1000.0
    vel_layered[cell.id()] = vel

pg.show(mesh_layered, vel_layered, label="Velocity (m/s)")

################################################################################
# We now define the analytical solution. The traveltime at a given offset `x` is
# the minimum of the direct and critically refracted wave, where the latter is
# governed by Snell's law.

def analyticalSolution2Layer(x, zlay=25, v1=1000, v2=3000):
    """Analytical solution for 2 layer case."""
    tdirect = np.abs(x) / v1  # direct wave
    alfa = asin(v1 / v2)  # critically refracted wave angle
    xreflec = tan(alfa) * zlay * 2.  # first critically refracted
    trefrac = (x - xreflec) / v2 + xreflec * v2 / v1**2
    return np.minimum(tdirect, trefrac)

################################################################################
# Vertical gradient model
# -----------------------
# We first create an unstructured mesh:

sensors = np.arange(131, step=10)
plc = mt.createWorld([-20, -60], [150, 0], worldMarker=False)
for pos in sensors:
    plc.createNode([pos, 0.0])
mesh_gradient = mt.createMesh(plc, quality=33, area=3)

################################################################################
# A vertical gradient model, i.e. :math:`v(z) = a + bz`, is defined per cell.

a = 1000
b = 100

vel_gradient = []
for node in mesh_gradient.nodes():
    vel_gradient.append(a + b * abs(node.y()))
vel_gradient = pg.meshtools.nodeDataToCellData(mesh_gradient,
                                               np.array(vel_gradient))
pg.show(mesh_gradient, vel_gradient, label="Velocity (m/s)")

################################################################################
# The traveltime for a gradient velocity model is given by:
#
# .. math::
#
#     v = \left|b^{-1}cosh^{-1}\left(1 + \frac{b^2 x^2}{2a^2}\right)\right|
#

def analyticalSolutionGradient(x, a=1000, b=100):
    """Analytical solution for gradient model."""
    tdirect = np.abs(x) / a  # direct wave
    tmp = 1 + ((b**2 * np.abs(x)**2) / (2 * a**2))
    trefrac = np.abs(b**-1 * np.arccosh(tmp))
    return np.minimum(tdirect, trefrac)

################################################################################
# The loop below calculates the traveltimes and makes the final comparison plot.

fig, ax = plt.subplots(3, 2, figsize=(10, 10), sharex=True)

for j, (case, mesh, vel) in enumerate(
        zip(["layered", "gradient"],
            [mesh_layered, mesh_gradient],
            [vel_layered, vel_gradient])):
    pg.boxprint(case)
    if case is "gradient":
        ana = analyticalSolutionGradient
    elif case is "layered":
        ana = analyticalSolution2Layer
    for boundary in mesh.boundaries():
        boundary.setMarker(0)

    xmin, xmax = mesh.xmin(), mesh.xmax()
    mesh.createNeighbourInfos()

    # In order to use the Dijkstra, we extract the surface positions >0
    mx = pg.x(mesh.positions())
    my = pg.y(mesh.positions())
    fi = pg.find((my == 0.0))
    px = np.sort(mx(fi))

    # A data container with index arrays named s (shot) and g (geophones) is
    # created and filled with the positions and shot/geophone indices.
    data = pg.DataContainer()
    data.registerSensorIndex('s')
    data.registerSensorIndex('g')

    for i, pxi in enumerate(px):
        data.createSensor([pxi, 0.0])
        if pxi == 0:
            source = i

    ndata = len(px) - 1
    data.resize(ndata)
    data.set('s', pg.RVector(ndata, source))  # only one shot at first sensor
    data.set('g', pg.utils.grange(1, ndata, 1))  # all others and geophones

    # We compare the accuracy for 0-5 secondary nodes
    sec_nodes = [0, 1, 3, 5]
    t_all = []
    durations = []
    paths = []

    for n in sec_nodes:
        if n > 0:
            pg.tic()
            mesh2 = mesh.createMeshWithSecondaryNodes(n)
            pg.toc("Mesh generation with %d secondary nodes:" % n)
        else:
            mesh2 = mesh

        # Perform traveltime calculations and log time with pg.tic() and pg.toc()
        pg.tic()
        fop = pg.TravelTimeDijkstraModelling(mesh2, data)
        t_all.append(fop.response(1 / vel))
        durations.append(pg.dur())
        pg.toc("Raytracing with %d secondary nodes:" % n)

        # This is to show single raypaths.
        dij = pg.Dijkstra(fop.createGraph(1 / vel))
        dij.setStartNode(mesh2.findNearestNode([0, 0]))
        paths_per_receiver = []
        for receiver in sensors:
            ni = dij.shortestPathTo(mesh2.findNearestNode([receiver, 0]))
            pos = mesh2.positions(withSecNodes=True)[ni]
            paths_per_receiver.append([pg.x(pos), pg.y(pos)])

        paths.append(paths_per_receiver)

    t_ana = ana(px[1:])

    # Upper subplot
    ax[1, j].plot(px[1:], t_ana * 1000, label="Analytical solution")

    for i, n in enumerate(sec_nodes):
        ax[1, j].plot(
            px[1:], t_all[i] * 1000,
            label="Dijkstra (%d sec nodes, %.2f s)" % (n, durations[i]))

    ax[2, j].plot(px, np.zeros_like(px),
                  label="Zero line")  # to keep color cycle

    for i, n in enumerate(sec_nodes):
        ax[2, j].plot(px[1:], np.abs(t_all[i] - t_ana) * 1000)

    ax[1, j].legend()

    pg.show(mesh, vel, ax=ax[0, j], label="Velocity (m/s)", hold=True,
            logScale=False, cMap="summer_r", coverage=0.7)
    drawMesh(ax[0, j], mesh, color="white", lw=0.21)

    recs = [1, 3, 8, 13]
    cols = ["orangered", "blue", "black"]
    to_plot = [0, 1, 3]
    for i, k in enumerate(to_plot):
        n = sec_nodes[k]
        for r, p in enumerate(recs):
            if r == 0:
                lab = "Raypath with %d sec nodes" % n
            else:
                lab = None

            ax[0, j].plot(paths[i][p][0], paths[i][p][1], cols[i], label=lab)
            ax[0, j].plot(sensors[p], 0.0, "kv", ms=10)

    ax[0, j].plot(0.0, 0.0, "ro", ms=10)
    ax[0, j].set_ylim(mesh.ymin(), 2)

ax[0, 0].set_title("Two-layer model")
ax[0, 1].set_title("Vertical gradient model")
ax[0, 0].legend()
ax[0, 0].set_ylabel("y (m)")
ax[1, 0].set_ylabel("Traveltime (ms)")
ax[2, 0].set_ylabel("Absolute difference to\nanalytical solution (ms)")
ax[2, 0].set_xlabel("x (m)")
fig.tight_layout()
