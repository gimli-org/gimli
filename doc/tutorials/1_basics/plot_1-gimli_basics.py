#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""

GIMLi Basics
============

This is the first tutorial where we demonstrate the general use of
:term:`GIMLi` in Python, i.e., :term:`pyGIMLi`.

The modelling as well as the inversion part of :term:`GIMLi` often requires a
spatial discretization for the domain of interest, the so called
:gimliapi:`GIMLI::Mesh`.
This tutorial shows some basic aspects of handling a mesh.

First, the library needs to be imported.
To avoid name clashes with other libraries we suggest to ``import pygimli`` and
alias it to the simple abbreviation ``pg``: CR
"""
import pygimli as pg


###############################################################################
# Every part of the c++ namespace :gimliapi:`GIMLI` is bound to python and can
# be used with the leading ``pg.``
#
# For instance get the current version for :term:`GIMLi` with:

print(pg.__version__)

###############################################################################
# This yields:
#
# Now we know the name space :gimliapi:`GIMLI` and can create a first mesh.
# A mesh is represented by a collection of nodes, cells and boundaries,
# i.e., geometrical entities.
#
# .. note::
#
#     A regularly spaced mesh consisting of rectangles or hexahedrons is
#     usually called grid. However, a grid is just a special variant of a mesh
#     so GIMLi treat it the same. The only difference is how they are created.
#
# GIMLi provides a collection of tools for mesh import, export and generation.
# A simple grid generation is built-in but we also provide wrappers for
# unstructured mesh generations, e.g., :term:`Triangle`, :term:`Tetgen` and
# :term:`Gmsh`. To create a 2d grid you need to give two arrays/lists for the
# sampling points in x and y direction, respectively, or just numbers.

grid = pg.createGrid(x=[-1.0, 0.0, 1.0, 4.0], y=(-1.0, 0.0, 1.0, 4.0))


###############################################################################
# The returned object ``grid`` is an instance of :gimliapi:`GIMLI::Mesh` and
# provides various methods for modification and io-operations. General
# informations can be simply printed.
#
print(grid)

###############################################################################
# yields:
#
# Or you can access them manually:
#
print('Mesh: Nodes:', grid.nodeCount(),
      'Cells:', grid.cellCount(),
      'Boundaries:', grid.boundaryCount())

###############################################################################
#
# You can iterate through all cells of the general type :gimliapi:`GIMLI::Cell`
# that also provides a lot of methods. Here we list the number of nodes and the
# node ids per cell:

for cell in grid.cells():
    print("Cell", cell.id(), "has", cell.nodeCount(),
          "nodes. Node IDs:", [n.id() for n in cell.nodes()])

print(type(grid.cell(0)))

###############################################################################
# To find the grid generation input arrays ``x`` and ``y``, you can use the
# build-in :gimliapi:`GIMLI::Vector` (pre-defined with value type double as
# ``pg.Vector``), standard python lists or :term:`numpy` arrays,
# which are widely compatible with :term:`GIMLi` vectors.

import numpy as np

grid = pg.createGrid(x=np.linspace(-1.0, 1.0, 10),
                     y=1.0 - np.logspace(np.log10(1.0), np.log10(2.0), 10))

###############################################################################
#
# We can find that this new ``grid`` contains
#

print(grid.cellCount())

###############################################################################
# rectangles of type :gimliapi:`GIMLI::Quadrangle` being derived from the
# base type :gimliapi:`GIMLI::Cell`,
#

print(grid.boundaryCount())

###############################################################################
# edges of type :gimliapi:`GIMLI::Edge`, which are boundaries of the general
# type :gimliapi:`GIMLI::Boundary`.
#
# The mesh can be saved and loaded in our binary mesh format ``.bms``.
# Or exported into ``.vtk`` format for 2D or 3D visualization using
# :term:`Paraview`.
#
# However, we recommend visualizing 2-dimensional content using python scripts
# that provides better exports to graphics files (e.g., png, pdf, svg).
# In :term:`pygimli` we provide some basic post-processing routines using
# the :term:`matplotlib` visualization framework. The main visualization call
# is :py:mod:`pygimli.viewer.show` which is sufficient for the most meshs,
# fields, models and streamline views.

pg.viewer.show(grid)

###############################################################################
# For more control you can also use the appropriate draw methods
# :py:mod:`pygimli.viewer.mpl.drawMesh`.

pg.wait()
