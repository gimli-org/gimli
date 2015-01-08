#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

GIMLi Basics
------------

This is the first tutorial where we demonstrate the general use of :term:`GIMLi` in Python, i.e., :term:`pyGIMLi`.

The modelling as well as the inversion part of :term:`GIMLi` often requires discretization, so handling meshes is an important aspect of this tutorial.

First, the library needs to be imported.
To avoid name clashes with other libraries we suggest to ``import pygimli`` and alias it to a simple abbreviation ``pg``, e.g., by using

"""
import pygimli as pg


"""
Every part of the c++ namespace :gimliapi:`GIMLI` is bind to python and can be now used with the leading ``pg.``

For instance get the current version for :term:`GIMLi` with:
"""

print(pg.__version__)

"""
This yields:

.. lastcout::

Now we know the namespace :gimliapi:`GIMLI` name space and can start creating a first mesh.
A mesh is spatial discretization and is represented by a collection of nodes, cells and boundaries, i.e., geometrical entities.

.. note::

    A regular spaced mesh consisting of rectangles or hexahedrons is usually called grid. 
    However, since a grid is just a kind of a mesh GIMLi does not separate between.
    The only difference is the way of creation.

GIMLi provides a collection of tools for mesh import, export and generation.
A simple grid generation is built-in but we also provide wrappers for unstructured mesh generations, e.g., :term:`Triangle`, :term:`tetgen` and :term:`Gmsh`.
Incorporation of other generators should be straight forward, e.g., :py:mod:`pygimli.meshtools`
"""
import numpy as np
grid = pg.createGrid(x=np.linspace(-1.0, 1.0, 3), y=-np.linspace(-1.0, 1.0, 3))
#grid = pg.createGrid(x=[-1.0, 0.0, 1.0, 4.0], y=[-1.0, 0.0, 1.0, 4.0])


"""
``grid`` is an instance of :gimliapi:`GIMLI::Mesh` and provides various methods for modification and io-operations.

"""
print(grid)

"""
yields:

.. lastcout::

For instance, you can iterate through all elements of the general type :gimliapi:`GIMLI::Cell`, which in turn also provides a lot of methods:

"""

for cell in grid.cells():
    print(cell.id(), type(cell), cell.nodeCount())

"""
.. lastcout::

To define the grid generation input arrays ``x`` and ``y``, you can also use the build-in :gimliapi:`GIMLI::Vector` (pre-defined with value type double as ``pg.RVector``), standard python lists or of course :term:`numpy` arrays, which are widely compatible with :term:`GIMLi` vectors.

"""

import numpy as np


grid = pg.createGrid(x=np.linspace(-1.0, 1.0, 10),
                    y=1.0 - np.logspace(np.log10(1.0), np.log10(2.0), 10))


"""

We can found that this new ``grid`` contains

"""

print((grid.cellCount()))

"""

.. lastcout::

rectangles of type :gimliapi:`GIMLI::Quadrangle` being derived from the base type :gimliapi:`GIMLI::Cell`,

"""

print((grid.boundaryCount()))

"""

.. lastcout::

edges of type :gimliapi:`GIMLI::Edge`, which are boundaries of the general type :gimliapi:`GIMLI::Boundary`.

The mesh can be saved/loaded in/from our binary mesh format ``.bms``. Or exported in ``.vtk`` format for 2D or 3D visualization using
:term:`Paraview`.

However, we recommend visualizing 2-dimensional content using python scripts that provides better exports to graphics files (e.g. png, pdf, svg).
In :term:`pygimli` we provide some basic post-processing routines using the :term:`matplotlib` visualization framework.
See: :py:mod:`pygimli.viewer`

"""


pg.show(grid)

"""
.. image:: PLOT2RST.current_figure

"""
