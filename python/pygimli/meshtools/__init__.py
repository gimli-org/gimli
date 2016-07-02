# encoding: utf-8
"""
meshtools
=========

This package provides functions for mesh generation and modification.

.. note::

    Although we discriminate here between grids (structured meshes) and meshes
    (unstructured), both objects are treated the same internally.
"""

from .grid import appendTetrahedronBoundary, appendTriangleBoundary
from .mesh import (createMesh, createParaMesh, createParaMesh2DGrid,
                   merge2Meshes, mergeMeshes, readGmsh, readHydrus2dMesh,
                   readHydrus3dMesh, readTetgen, readTriangle)
from .polytools import createParaDomain2D  # keep for backward compatibility
from .polytools import (createCircle, createLine, createParaMeshPLC,
                        createRectangle, createWorld, mergePLC, readPLC)

#  This is neither functional nor good practice
#  __all__ = [name for name in dir() if '_' not in name]

__all__ = ['appendTriangleBoundary',
           'appendTetrahedronBoundary',
           'createMesh',
           'readGmsh',
           'readTriangle',
           'readTetgen',
           'readHydrus2dMesh',
           'readHydrus3dMesh',
           'mergeMeshes',
           'merge2Meshes',
           'createParaMesh',
           'createParaMesh2DGrid',
           'createRectangle',
           'createWorld',
           'createCircle',
           'createLine',
           'createParaMeshPLC',
           'mergePLC',
           'readPLC',
           'createParaDomain2D'  # keep for backward compatibility
           ]
