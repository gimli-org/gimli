#!/usr/bin/env python
# encoding: utf-8
r"""
Building a hybrid mesh in 2D
----------------------------

In some cases, the modelling domain may require flexibility in one region and
equidistant structure in another. In this short example, we demonstrate how to
accomplish this for a two-dimensional mesh consisting of a region with regularly
spaced quadrilaterals and a region with unstructured triangles."""

###############################################################################
# We start by importing numpy, matplotlib and pygimli with its required components.

import numpy as np

import pygimli as pg
from pygimli.meshtools import appendTriangleBoundary, merge2Meshes
from pygimli.mplviewer import drawMesh
from pygimli.viewer import showMesh

###############################################################################
# We continue by building a regular grid and assign the marker 2 to all cells.
xmin, xmax = 0., 50.
zmin, zmax = -50., -25.

dx = 1.0
xreg = np.arange(xmin, xmax + dx, dx, 'float')
zreg = np.arange(zmin, zmax + dx, dx, 'float')

mesh1 = pg.Mesh(2)
mesh1.create2DGrid(xreg, zreg, 0)
for c in mesh1.cells():
    c.setMarker(2)

print(mesh1)

###############################################################################
# Next, we build an unstructured region on top by creating the polygon and
# calling triangle via pygimli's TriangleWrapper.

poly = pg.Mesh(2)  # empty 2d mesh
nStart = poly.createNode(xmin, zmax, 0.0)

nA = nStart
for x in xreg[1:]:
    nB = poly.createNode(x, zmax, 0.0)
    poly.createEdge(nA, nB)
    nA = nB

z2 = 0.
nA = poly.createNode(xmax, z2, 0.0)
poly.createEdge(nB, nA)
nB = poly.createNode(xmin, z2, 0.0)
poly.createEdge(nA, nB)
poly.createEdge(nB, nStart)

tri = pg.TriangleWrapper(poly)
tri.setSwitches('-pzeAfaq31')

###############################################################################
# For more information on the triangle switches and the corresponding settings,
# the reader is referred to `the triangle website
# <http://www.cs.cmu.edu/~quake/triangle.switch.html>`_.
#
# Now we can generate the unstructured mesh.
mesh2 = pg.Mesh(2)
tri.generate(mesh2)

for cell in mesh2.cells():
    cell.setMarker(1)

###############################################################################
# Finally, the grid and the unstructured mesh can be merged to single mesh for
# further modelling.

mesh3 = merge2Meshes(mesh1, mesh2)

###############################################################################
# Of course, you can treat the hybrid mesh like any other mesh and append a
# triangle boundary for example with the function
# :py:func:`pygimli.meshtools.grid.appendTriangleBoundary`.

mesh = appendTriangleBoundary(mesh3, -100., 100., quality=31, smooth=True,
                              marker=3, isSubSurface=True)

ax, cbar = showMesh(mesh, mesh.cellMarkers(), cmap="summer",
                    label="Region marker")

drawMesh(ax, mesh)

ax, _ = showMesh(mesh, mesh.cellMarkers(), logScale=False,
                 label="Region marker")

drawMesh(ax, mesh)
pg.wait()
