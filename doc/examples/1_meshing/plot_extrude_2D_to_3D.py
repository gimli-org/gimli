#!/usr/bin/env python
# encoding: utf-8
r"""
Extrude a 2D mesh to 3D
=======================

This example shows how to extrude a 2D mesh to 3D. This can be helpful for
closed laboratory geometries for example. If you are looking for more flexible
ways to create 3D meshes, have a look at TetGen and Gmsh.
"""

import numpy as np

import pygimli as pg
import pygimli.meshtools as mt

###############################################################################
# We start by generating a 2D mesh.
plc = mt.createCircle([-1, -4], radius=1.5, area=0.1, segments=25)
circle = mt.createMesh(plc)
for cell in circle.cells():
    cell.setMarker(cell.id())
###############################################################################
# We now extrude this mesh to 3D given a *z* vector.

z = np.geomspace(1, 5, 5)
cylinder = pg.createMesh3D(circle, z)
pg.show(cylinder, cylinder.cellMarkers(), label="Cell markers")
