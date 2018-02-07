# encoding: utf-8
"""
Mesh generation and modification.

.. note::

    Although we discriminate here between grids (structured meshes) and meshes
    (unstructured), both objects are treated the same internally.
"""

from .grid import (createGrid, appendTetrahedronBoundary, appendTriangleBoundary)

from .mesh import (createMesh, createParaMesh, createParaMesh2DGrid,
                   merge2Meshes, refineQuad2Tri,
                   mergeMeshes, readGmsh, readHydrus2dMesh,
                   readHydrus3dMesh, readTetgen, readTriangle, convertHDF5Mesh,
                   readHDF5Mesh, readFenicsHDF5Mesh, exportHDF5Mesh,
                   exportFenicsHDF5Mesh)
from .polytools import createParaDomain2D  # keep for backward compatibility
from .polytools import (createCircle, createLine, createParaMeshPLC,
                        createPolygon, createRectangle, createWorld, mergePLC,
                        readPLC, exportPLC, writePLC)
from .quality import (quality)

from .mapping import (nodeDataToCellData,
                      cellDataToNodeData,
                      nodeDataToBoundaryData,
                      cellDataToBoundaryData,
                      fillEmptyToCellArray,
                      tapeMeasureToCoordinates,
                      interpolate,
                      interpolateAlongCurve
                      )

#  This is neither functional nor good practice  #  why?
#  __all__ = [name for name in dir() if '_' not in name]

__all__ = ['appendTriangleBoundary',
           'appendTetrahedronBoundary',
           'createMesh',
           'readGmsh',
           'readTriangle',
           'readTetgen',
           'readHydrus2dMesh',
           'readHydrus3dMesh',
           'readHDF5Mesh',
           'readFenicsHDF5Mesh',
           'refineQuad2Tri',
           'mergeMeshes',
           'merge2Meshes',
           'createParaMesh',
           'createParaMesh2DGrid',
           'createPolygon',
           'createRectangle',
           'createWorld',
           'createCircle',
           'createLine',
           'createParaMeshPLC',
           'convertHDF5Mesh',
           'exportHDF5Mesh',
           'exportFenicsHDF5Mesh',
           'mergePLC',
           'readPLC',
           'writePLC',
           'exportPLC',
           'createParaDomain2D',  # keep for backward compatibility
           'quality'
           ]
