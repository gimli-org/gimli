"""
Import and extensions of the core Mesh class.
"""

from .._logger import deprecated
from ._pygimli_ import Mesh, MeshEntity, Node, PolygonFace, Line

def Mesh_str(self):
    return ("Mesh: Nodes: " + str(self.nodeCount()) + " Cells: " + str(
        self.cellCount()) + " Boundaries: " + str(self.boundaryCount()))

def MeshEntity_str(self):
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

def Node_str(self):
    """Give node infos."""
    s = self.__repr__()
    s += '\tID: ' + str(self.id()) + \
         ', Marker: ' + str(self.marker())
    s += '\t' + str(self.pos()) + '\n'
    return s

Node.__str__ = Node_str
Mesh.__str__ = Mesh_str
MeshEntity.__str__ = MeshEntity_str

def __MeshGetCellMarker__(self):
    deprecated(msg='Mesh::cellMarker()', hint='Mesh::cellMarkers()')
    return self.cellMarkers()


def __MeshSetCellMarker__(self, m):
    deprecated(msg='Mesh::setCellMarker()', hint='Mesh::setCellMarkers()')
    return self.setCellMarkers(m)

Mesh.cellMarker = __MeshGetCellMarker__
Mesh.setCellMarker = __MeshSetCellMarker__

def createSecondaryNodes(self, n=3):
    """Create `n` equally distributed secondary nodes on boundaries of the mesh.
    This is useful to increase the accuracy of traveltime calculations.

    TODO:
        * Fill on 2D boundaries

    Parameters
    ----------
    n : int
        Number of secondary nodes (the default is 3).

    Returns
    -------
    pg.Mesh
        Mesh with secondary nodes.
    """
    secMesh = Mesh(self)
    secMesh.createNeighbourInfos()
    for b in secMesh.boundaries():
        A = b.node(0).pos()
        B = b.node(1).pos()
        line = Line(A, B)
        for i in range(n):
            secNode = secMesh.createNode(line.at((i + 1) / (n + 1)))
            b.addSecondaryNode(secNode)
    return secMesh

Mesh.createSecondaryNodes = createSecondaryNodes
