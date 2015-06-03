#!/usr/bin/env python
# encoding: utf-8

"""
Building a hybrid mesh in 2-D
-----------------------------

In some cases, the modelling domain may require flexibility in one region and
equidistant structure in another. In this short example, we demonstrate how to
accomplish this for a two-dimensional mesh consisting of a region with regularly
spaced quadrilaterals and a region with unstructured triangles.

We start by importing numpy, matplotlib and pygimli with its required components.
"""
import numpy as np
from matplotlib import pyplot as plt

import pygimli as pg
from pygimli.viewer import showMesh
from pygimli.mplviewer import drawMesh
from pygimli.meshtools import merge2Meshes, appendTriangleBoundary

"""
We continue by building a regular grid and assign the marker 2 to all cells.
"""
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

"""
.. lastcout::

Next, we build an unstructured region on top by creating the polygon and calling
triangle via pygimli's TriangleWrapper.
"""

# append rectangle above, search upper lines
poly = pg.Mesh(2)  # empty 2d mesh
n1 = poly.createNode(xmin, zmax, 0.0)
n0 = pg.Node(n1)  # make a copy for later
for x in xreg[1:]:
    n2 = poly.createNode(x, zmax, 0.0)
    poly.createEdge(n1, n2)
    n1 = n2

z2 = 0.
n2 = poly.createNode(xmax, z2, 0.0)
poly.createEdge(n1, n2)
n1 = poly.createNode(xmin, z2, 0.0)
poly.createEdge(n1, n2)
poly.createEdge(n1, n0)

tri = pg.TriangleWrapper(poly)
tri.setSwitches('-pzeAfaq31')

"""
For more information on the triangle switches and the corresponding settings,
the reader is referred to `the triangle website <http://www.cs.cmu.edu/~quake/triangle.switch.html>`_.

Now we can generate the unstructured mesh.
"""
mesh2 = pg.Mesh(2)
tri.generate(mesh2)

for cell in mesh2.cells():
    cell.setMarker(1)

"""
.. lastcout::

Finally, the grid and the unstructured mesh can be merged to single mesh for further
modelling.
"""
pg.show(mesh1)

# looking for *** Error in `/usr/bin/python3.3': corrupted double-linked list: 0x00000000022a0ac0 ***
mesh3 = merge2Meshes(mesh1, mesh2)
pg.show(mesh3)

"""
.. lastcout::

Of course, you can treat the hybrid mesh like any other mesh and append a triangle
boundary for example with :py:func:`pygimli.meshtools.grid.appendTriangleBoundary`.
"""

mesh = appendTriangleBoundary(mesh3, -100., 100., quality=31,
                              smooth=True, marker=3, isSubSurface=True)

ax, cbar = showMesh(mesh, mesh.cellMarker(), 
                    cmap="summer",
                    label="Region marker", 
                    showLater=True)

drawMesh(ax, mesh)

ax, _ = showMesh(mesh, mesh.cellMarker(),
                 logScale=False,
                 label="Region marker",
                 showLater=True)

drawMesh(ax, mesh)

plt.xlim(40,60)
plt.ylim(-30, -20)
plt.show()


