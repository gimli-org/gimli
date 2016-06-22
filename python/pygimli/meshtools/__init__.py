# encoding: utf-8

"""
meshtools

This package provides functions for mesh generation and modification.

.. note::
    Although we discriminate here between grids (structured meshes) and meshes
    (unstructured), both objects are treated the same internally.
"""

from .grid import appendTriangleBoundary
from .grid import appendTetrahedronBoundary

from .mesh import createMesh
from .mesh import readGmsh
from .mesh import readTriangle
from .mesh import readTetgen
from .mesh import readHydrus2dMesh
from .mesh import readHydrus3dMesh
from .mesh import mergeMeshes
from .mesh import createParaMesh
from .mesh import createParaMesh2DGrid

from .polytools import createRectangle
from .polytools import createWorld
from .polytools import createCircle
from .polytools import createLine
from .polytools import createParaMeshPLC
from .polytools import mergePLC
from .polytools import readPLC

__all__ = ['appendTriangleBoundary',
           'appendTetrahedronBoundary',
           'createMesh',
           'readGmsh',
           'readTriangle',
           'readTetgen',
           'readHydrus2dMesh',
           'readHydrus3dMesh',
           'mergeMeshes',
           'createParaMesh',
           'createParaMesh2DGrid',
           'createRectangle',
           'createWorld',
           'createCircle',
           'createLine',
           'createParaMeshPLC',
           'mergePLC',
           'readPLC']
