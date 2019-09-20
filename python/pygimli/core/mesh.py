# -*- coding: utf-8 -*-
"""
Import and extensions of the core Mesh class.
"""

from ._pygimli_ import (HexahedronShape, Line, Mesh, MeshEntity, Node,
                        PolygonFace, TetrahedronShape, TriangleFace)
from .logger import deprecated, error, info, warn
from ..meshtools import mergePLC, exportPLC


def __Mesh_str(self):
    st = "Mesh: Nodes: " + str(self.nodeCount()) + " Cells: " + str(
        self.cellCount()) + " Boundaries: " + str(self.boundaryCount())
    if (self.secondaryNodeCount() > 0):
        st += " secNodes: " + str(self.secondaryNodeCount())

    return st


Mesh.__str__ = __Mesh_str


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

    if isinstance(self, PolygonFace):
        s += '\t' + str(self.nodeCount()) + " Nodes.\n"
    else:
        for n in self.nodes():
            s += '\t' + str(n.id()) + " " + str(n.pos()) + "\n"
    return s


MeshEntity.__str__ = __MeshEntity_str


def __Node_str(self):
    """Give node infos."""
    s = self.__repr__()
    s += '\tID: ' + str(self.id()) + \
         ', Marker: ' + str(self.marker())
    s += '\t' + str(self.pos()) + '\n'
    return s


Node.__str__ = __Node_str

# For Jupyer Notebook use.. checkme
# Node.__repr__ = Node_str
# Mesh.__repr__ = Mesh_str
# MeshEntity.__repr__ = MeshEntity_str


def __Mesh_setVal(self, key, val):
    """Index access to the mesh data"""
    self.addData(key, val)


Mesh.__setitem__ = __Mesh_setVal


def __Mesh_getVal(self, key):
    """Index access to the mesh data"""
    if self.haveData(key):
        return self.data(key)
    else:
        error('The mesh does not have the requested data:', key)


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
    self.createNeighbourInfos()

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
                    nmax = n
                    if isinstance(b, TriangleFace):
                        nmax = n - sx
                    for sy in range(nmax):
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

Mesh.exportPLC = exportPLC
