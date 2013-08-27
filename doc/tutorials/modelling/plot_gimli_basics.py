#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

GIMLi Basics
------------

*This introductory sentence should state the intent and goal of the tutorial. Keep it brief.*

This is the first tutorial where we actually want to show the general use of :term:`GIMLi`, i.e., :term:`pyGIMLi`. 

The modelling as well as the inversion part of :term:`GIMLi` is usually based on discretization, so handling meshes is an important aspect of this tutorial.

In the First, the library have to be imported.
To avoid name clashes with other libraries we suggest to ``import pygimli`` and alias it to a simple abbreviation ``g``, e.g., by using

"""
import pygimli as g

"""
Every part of the namespace :gimliapi:`GIMLI` is bind to python and can be now used with the leading ``g``

For instance get the current version for :term:`GIMLi` with:
"""

print(g.versionStr())

"""
This yield: 

.. lastcout::

Now we know the gimli name space and can start creating a first mesh. A Mesh is a collection of cells and boundaries, i.e, geometrical entities.

.. note::

    A regular spaced mesh consisting of rectangles or hexahedrons is called grid here. However, internal a grid is also a mesh.

Gimli provide a collection of tools for mesh generation. Easy grid generation is build in but we provide also some wrapper for unstructured meshs. e.g. :term:`triangle` and :term:`tetgen`. Incorporation of other generators should be straight forward but is discussed in :ref:`part:examples`.
"""

grid = g.createGrid(x=[-1.0, 0.0, 1.0], y=[-1.0, 0.0, 1.0])

"""
``grid`` is an instance of :gimliapi:`GIMLI::Mesh` and provide a couple methods for modifications and io-operations.

"""
print(grid)

"""
yield:

.. lastcout::

For instance you can iterate through all elements of type :gimliapi:`GIMLI::Cell`, which also provide a lot of methods, with:

"""

for cell in grid.cells():
    print cell.id(), type(cell), cell.nodeCount()

"""
.. lastcout::

To define the grid generation input arrays ``x`` and ``y``, you can also use the build-in :gimliapi:`GIMLI::Vector` (predefined with value type double as ``g.RVector``), standard python lists or of course :term:`numpy` arrays.

"""

import numpy as np

grid = g.createGrid(x=np.linspace(-1.0, 1.0, 10), 
                    y=1.0 - np.logspace(np.log10(1.0), np.log10(2.0), 10))


"""

We can found that this ``grid`` contains,

"""

print(grid.cellCount())

"""

.. lastcout::

rectangles of type :gimliapi:`GIMLI::Quadrangle`, which are of base type :gimliapi:`GIMLI::Cell`.

"""

print(grid.boundaryCount())

"""

.. lastcout::

edges of type :gimliapi:`GIMLI::Edge`, which are generally called boundaries of type :gimliapi:`GIMLI::Boundary`.

The mesh can be saved/loaded in/from our binary mesh format ``.bms``. Or exported in ``.vtk`` format to the visualization via
:term:`Paraview`

For sure we can visualize the most 2-dimensional content within python scripts. 
In :term:`pygimli` we provide some basic post-processing routines using the :term:`Pylab` visualization framework.

"""

from pygimli.viewer import showMesh

showMesh(grid)

"""
.. image:: PLOT2RST.current_figure

"""
