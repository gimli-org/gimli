# -*- coding: utf-8 -*-
"""
Import and extensions of the core Mesh class.
"""
import numpy as np

from math import ceil
from ._pygimli_ import (cat, HexahedronShape, Line,
                        Mesh, MeshEntity, Node, Boundary, RVector,
                        PolygonFace, TetrahedronShape, TriangleFace)
from .logger import deprecated, error, info, warn
from ..meshtools import mergePLC, exportPLC
from .base import isScalar, isArray, isPos, isR3Array, isComplex


def __Mesh_str(self):
    st = "Mesh: Nodes: " + str(self.nodeCount()) + " Cells: " + str(
        self.cellCount()) + " Boundaries: " + str(self.boundaryCount())

    if (self.secondaryNodeCount() > 0):
        st += " secNodes: " + str(self.secondaryNodeCount())

    if len(list(self.dataMap().keys())) > 0:
        st += "\nMesh contains data: "

        uniqueNames = {}
        for d in self.dataMap().keys():

            uName = d
            if '_x' in d:
                uName = uName[0:d.find('_x')]

            if '_y' in d or '_z' in d:
                continue

            if '#' in d:
                uName = uName[0:d.find('#')]


            if not uName in uniqueNames:
                uniqueNames[uName] = []

            uniqueNames[uName].append(d)

        for d, v in uniqueNames.items():
            if len(v) > 1:
                st += d + "[0,...,{}) ".format(len(v))
            else:
                st += d + " "

    return st
Mesh.__repr__ =__Mesh_str


def __addPLCs__(self, other):
    if self.isGeometry() and other.isGeometry():
        return mergePLC([self, other])
    else:
        error("Addition is only supported for PLCs, i.e. meshs without cells.")
Mesh.__add__ = __addPLCs__


def __MeshEntity_str(self):
    """Give mesh entity infos."""
    s = self.__repr__()
    s += '\tID: ' + str(self.id()) + \
         ', Marker: ' + str(self.marker()) + \
         ', Size: ' + str(self.size()) + '\n'

    if isinstance(self, PolygonFace) and len(self.nodes()) > 5:
        s += '\t' + str(self.nodeCount()) + " Nodes.\n"
    else:
        for n in self.nodes():
            s += '\t' + str(n.id()) + " " + str(n.pos()) + "\n"
    return s
MeshEntity.__repr__ =__MeshEntity_str


def __Node_str(self):
    """Give node infos."""
    s = '\tID: ' + str(self.id()) + \
         ', Marker: ' + str(self.marker())
    s += '\t' + str(self.pos()) + '\n'
    return s
Node.__repr__ =__Node_str

def __Mesh_setVal(self, key, val):
    """Index access to the mesh data.

    Multiple arrays via matrix will be saved too.
    """
    # print(key, len(val), isR3Array(val))

    if isR3Array(val):
        return self.addData(key, val)

    if isinstance(val, list) and isinstance(val[0], (RVector, np.ndarray)) or \
        val.ndim == 2:

        #print(val.ndim)
        maxDigit = ceil(np.log10(len(val)))

        for i, v in enumerate(val):
            #print(i, v, maxDigit, '{}#{}'.format(key, str(i).zfill(maxDigit)))

            self.addData('{}#{}'.format(key, str(i).zfill(maxDigit)),
                         np.asarray(v))
    else:
        self.addData(key, val)

    #print('keys', self.dataMap.keys())
Mesh.__setitem__ = __Mesh_setVal


def __Mesh_getVal(self, key):
    """Index access to the mesh data"""
    if self.haveData(key):
        return self.data(key)
    else:

        uniqueNames = {}
        for d in self.dataMap().keys():
            if '_y' in d or '_z' in d:
                continue

            uName = d

            if '_x' in uName:
                uName = uName[0:uName.find('_x')]
                d1 = self.data(d)
                d2 = self.data(d.replace('_x', '_y'))
                d3 = self.data(d.replace('_x', '_z'))
                dat = np.array([d1, d2, d3]).T
            else:
                dat = self.data(d)

            if '#' in uName:
                uName = uName[0:uName.find('#')]

            if not uName in uniqueNames:
                uniqueNames[uName] = []

            uniqueNames[uName].append(dat)

        if key in uniqueNames:
            if len(uniqueNames[key]) == 1:
                return uniqueNames[key][0]
            try:
                return np.array(uniqueNames[key])
            except:
                return uniqueNames[key]

        critical('The mesh does not have the requested data:', key,
              '. Available:', uniqueNames)


Mesh.__getitem__ = __Mesh_getVal


def __MeshBoundingBox__(self):
    bb = self.boundingBox()
    mi = [bb.min()[i] for i in range(self.dim())]
    ma = [bb.max()[i] for i in range(self.dim())]
    return [mi, ma]
Mesh.bb = __MeshBoundingBox__


def __MeshGetCellMarker__(self):
    deprecated(msg='Mesh::cellMarker()', hint='Mesh::cellMarkers()')
    return self.cellMarkers()


def __MeshSetCellMarker__(self, m):
    deprecated(msg='Mesh::setCellMarker()', hint='Mesh::setCellMarkers()')
    return self.setCellMarkers(m)


def __MeshHoleMarkers__(self):
    return self.holeMarker()

Mesh.cellMarker = __MeshGetCellMarker__
Mesh.setCellMarker = __MeshSetCellMarker__
Mesh.holeMarkers = __MeshHoleMarkers__


def __createSecondaryNodes__(self, n=3, verbose=False):
    """Create `n` equally distributed secondary nodes on the mesh boundaries.
    This is useful to increase the accuracy of traveltime calculations.

    Parameters
    ----------
    n : int
        Number of secondary nodes (the default is 3).
    verbose : bool
        Optionally output number of added nodes.

    Returns
    -------
    pg.Mesh
        Copy of the given mesh with secondary nodes.
    """
    self.createNeighborInfos()

    if self.boundary(0).nodeCount() != self.boundary(0).allNodeCount():
        warn("Mesh already contains secondary nodes. Not adding any more.")
    else:
        if self.dim() == 2:
            for b in self.boundaries():
                A = b.node(0).pos()
                B = b.node(1).pos()
                line = Line(A, B)
                for i in range(n):
                    sn = self.createSecondaryNode(line.at((i + 1) / (n + 1)))
                    b.addSecondaryNode(sn)
        elif self.dim() == 3:
            for b in self.boundaries():
                bs = b.shape()
                for sx in range(n):
                    nMax = n
                    if isinstance(b, TriangleFace):
                        nMax = n - sx
                    for sy in range(nMax):
                        if isinstance(b, TriangleFace):
                            pos = bs.xyz([(sx + 1) / (n + 2),
                                          (sy + 1) / (n + 2)])
                        else:
                            pos = bs.xyz([(sx + 1) / (n + 1),
                                          (sy + 1) / (n + 1)])

                        sn = self.createSecondaryNode(pos)
                        b.addSecondaryNode(sn)

            for c in self.cells():
                # add secondary nodes to the edges of 3 Entities

                edges = []
                if isinstance(c.shape(), HexahedronShape):
                    #   7------6
                    #  /|     /|
                    # 4------5 |
                    # | 3----|-2
                    # |/     |/
                    # 0------1
                    edges.append([c.shape().node(0), c.shape().node(1)])
                    edges.append([c.shape().node(1), c.shape().node(2)])
                    edges.append([c.shape().node(2), c.shape().node(3)])
                    edges.append([c.shape().node(3), c.shape().node(0)])

                    edges.append([c.shape().node(0), c.shape().node(4)])
                    edges.append([c.shape().node(1), c.shape().node(5)])
                    edges.append([c.shape().node(2), c.shape().node(6)])
                    edges.append([c.shape().node(3), c.shape().node(7)])

                    edges.append([c.shape().node(4), c.shape().node(5)])
                    edges.append([c.shape().node(5), c.shape().node(6)])
                    edges.append([c.shape().node(6), c.shape().node(7)])
                    edges.append([c.shape().node(7), c.shape().node(4)])
                elif isinstance(c.shape(), TetrahedronShape):
                    edges.append([c.shape().node(0), c.shape().node(1)])
                    edges.append([c.shape().node(0), c.shape().node(2)])
                    edges.append([c.shape().node(0), c.shape().node(3)])

                    edges.append([c.shape().node(1), c.shape().node(2)])
                    edges.append([c.shape().node(2), c.shape().node(3)])
                    edges.append([c.shape().node(3), c.shape().node(1)])
                else:
                    print(c)
                    warn('cell type unknown')

                for e in edges:
                    line = Line(e[0].pos(), e[1].pos())
                    for i in range(n):
                        sn = self.createSecondaryNode(line.at((i+1)/(n+1)),
                                                      tol=1e-6)
                        c.addSecondaryNode(sn)
        else:
            warn("Unknown dimension. Don't know what to do.")

    if verbose:
        info("Added %d secondary nodes." % self.secondaryNodeCount())


def __createMeshWithSecondaryNodes__(self, n=3, verbose=False):
    m = Mesh(self)
    m.createSecondaryNodes(n, verbose)
    return m
Mesh.createSecondaryNodes = __createSecondaryNodes__
Mesh.createMeshWithSecondaryNodes = __createMeshWithSecondaryNodes__


__Mesh_deform__ = Mesh.deform
def __deform__(self, eps, mag=1.0):
    v = None
    dof = self.nodeCount()
    if hasattr(eps, 'ndim') and eps.ndim == 1:
        v = eps
    elif len(eps) == self.dim():
        if len(eps[0]) == dof:
            if self.dim() == 2:
                v = cat(eps[0], eps[1])
            elif self.dim() == 3:
                v = cat(cat(eps[0], eps[1]), eps[2])
            else:
                v = eps[0]
        else:
            print(self)
            print(len(eps), len(eps[0]))
            error('Size of displacement does not match mesh nodes size.')
    elif len(eps) == self.nodeCount() and eps.ndim == 2:
        v = eps.reshape(self.nodeCount() * eps.shape[1], order='F')

    return __Mesh_deform__(self, v, mag)

Mesh.deform = __deform__

Mesh.exportPLC = exportPLC

# just to keep backward compatibility 20191120
Mesh.createNeighbourInfos = Mesh.createNeighborInfos
Mesh.xmin = Mesh.xMin
Mesh.ymin = Mesh.yMin
Mesh.zmin = Mesh.zMin
Mesh.xmax = Mesh.xMax
Mesh.ymax = Mesh.yMax
Mesh.zmax = Mesh.zMax

def __Boundary_outside__(self):
    """Is the boundary is on the outside of the mesh."""
    return self.leftCell() is not None and self.rightCell() is None

Boundary.outside = __Boundary_outside__

def __Mesh_h__(self):
    return np.array([c.shape().h() for c in self.cells()])
Mesh.h = __Mesh_h__
