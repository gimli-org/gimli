# encoding: utf-8
"""
Mesh generation and modification.

.. note::

    Although we discriminate here between grids (structured meshes) and meshes
    (unstructured), both objects are treated the same internally.
"""

from ..core import createMesh1D, createMesh1DBlock, createMesh2D, createMesh3D
from .grid import appendTetrahedronBoundary, appendTriangleBoundary, createGrid
from .mapping import (cellDataToBoundaryData, cellDataToNodeData,
                      fillEmptyToCellArray, interpolate, interpolateAlongCurve,
                      nodeDataToBoundaryData, nodeDataToCellData,
                      tapeMeasureToCoordinates)

from .mesh import (convert, convertMeshioMesh, convertHDF5Mesh, createMesh,
                   createParaMesh,
                   createParaMesh2DGrid, exportFenicsHDF5Mesh, exportHDF5Mesh,
                   exportSTL, extrudeMesh, merge2Meshes, mergeMeshes,
                   readFenicsHDF5Mesh, readGmsh, readHDF5Mesh,
                   readHydrus2dMesh, readHydrus3dMesh, readSTL, readTetgen,
                   readTriangle, refineHex2Tet, refineQuad2Tri)
from .polytools import createParaDomain2D  # keep for backward compatibility
from .polytools import (createCircle, createCube, createCylinder, createFacet,
                        createLine, createParaMeshPLC, createPolygon,
                        createRectangle, createWorld, exportPLC, mergePLC,
                        mergePLC3D, readPLC, syscallTetgen, writePLC)
from .quality import quality

#  This is neither functional nor good practice  #  why?
#  __all__ = [name for name in dir() if '_' not in name]

__all__ = [
    'appendTriangleBoundary',
    'appendTetrahedronBoundary',
    'createMesh',
    'readGmsh',
    'readTriangle',
    'readTetgen',
    'readHydrus2dMesh',
    'readHydrus3dMesh',
    'readHDF5Mesh',
    'readFenicsHDF5Mesh',
    'readSTL',
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
