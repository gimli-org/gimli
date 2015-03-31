#!/usr/bin/env python
# encoding: utf-8

"""
Mesh interpolation
==================

In this tutorial, we look at the mesh interpolation options in GIMLi. Although,
the example shown here is in 2D, the same routines can be applied when
converting 3D data to a 2D mesh for instance.
"""

import numpy as np
import pygimli as pg
import matplotlib.pyplot as plt

from pygimli.meshtools import createParaMesh2dGrid
from pygimli.mplviewer import drawMesh, drawModel


"""
Create coarse mesh
------------------
"""

n = 5
nc = np.linspace(0, -2, n)
coarse = pg.createMesh2D(nc, nc)
mcx = pg.x(coarse.cellCenter())
mcy = pg.y(coarse.cellCenter())
coarse_data = np.cos(1.5 * mcx) * np.sin(1.5 * mcy)
fig, ax = plt.subplots()
drawModel(ax, coarse, coarse_data)
drawMesh(ax, coarse)

"""
Create fine mesh
----------------
"""

nf = np.linspace(0, -2, 4*n)
fine = pg.createMesh2D(nf, nf)
mfx = pg.x(fine.cellCenter())
mfy = pg.y(fine.cellCenter())
fine_data = np.cos(1.5 * mfx) * np.sin(1.5 * mfy)
fig, ax = plt.subplots()
drawModel(ax, fine, fine_data)
drawMesh(ax, fine)


"""
Interpolate data to different meshes
------------------------------------

We define two functions that take the input mesh, the input data and the output
mesh as parameters and return the input data interpolated to the output mesh.
"""

def nearest_neighbor_interpolation(inmesh, indata, outmesh, nan=99.9):
    """ Nearest neighbor interpolation. """
    outdata = []
    for pos in outmesh.cellCenters():
        cell = inmesh.findCell(pos)
        if cell:
            outdata.append(indata[cell.id()])
        else:
            outdata.append(nan)
    return outdata

def linear_interpolation(inmesh, indata, outmesh):
    """ Linear interpolation using `pg.interpolate()` """
    outdata = pg.RVector() # empty
    pg.interpolate(inmesh, indata, outmesh.cellCenters(), outdata)
    return outdata

"""
Visualization
-------------
"""

meshes = [coarse, fine]
datasets = [coarse_data, fine_data]
ints = [nearest_neighbor_interpolation, linear_interpolation]

fig, axes = plt.subplots(2,2, figsize=(5,5))

# Coarse data to fine mesh
drawModel(axes[0,0], fine, ints[0](coarse, coarse_data, fine), showCbar=False)
drawMesh(axes[0,0], fine)
drawModel(axes[0,1], fine, ints[1](coarse, coarse_data, fine), showCbar=False)
drawMesh(axes[0,1], fine)

# Fine data to coarse mesh
drawModel(axes[1,0], coarse, ints[0](fine, fine_data, coarse), showCbar=False)
drawMesh(axes[1,0], coarse)
drawModel(axes[1,1], coarse, ints[1](fine, fine_data, coarse), showCbar=False)
drawMesh(axes[1,1], coarse)

titles = ["Coarse to fine\nwith nearest neighbors", "Coarse to fine\nwith linear interpolation",
          "Fine to coarse\nwith nearest neighbors", "Fine to coarse\nwith linear interpolation"]

for a, title in zip(axes.flat, titles):
    a.set_title(title + "\n")

fig.tight_layout()
