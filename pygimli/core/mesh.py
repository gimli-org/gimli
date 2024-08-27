# -*- coding: utf-8 -*-
"""
Import and extensions of the core Mesh class.
"""
import numpy as np
from math import ceil
from .core import (
    cat,
    HexahedronShape,
    Line,
    RSparseMapMatrix,
    Mesh,
    MeshEntity,
    Node,
    Boundary,
    RVector,
    RVector3,
    PolygonFace,
    TetrahedronShape,
    TriangleFace,
)
from .logger import deprecated, error, info, warn, critical
from .base import isScalar, isArray, isPos, isR3Array, isComplex
from ..meshtools import mergePLC, exportPLC, interpolateAlongCurve


def __Mesh_unique_dataKeys(self):
    """Return unique data keys"""
    uniqueNames = {}
    for d in self.dataMap().keys():

        uName = d
        if "_x" in d:
            uName = uName[0 : d.find("_x")]

        if "_y" in d or "_z" in d:
            continue

        if "#" in d:
            uName = uName[0 : d.find("#")]

        if not uName in uniqueNames:
            uniqueNames[uName] = []

        uniqueNames[uName].append(d)
    return uniqueNames


Mesh.dataKeys = __Mesh_unique_dataKeys


def __Mesh_str(self):
    st = (
        "Mesh: Nodes: "
        + str(self.nodeCount())
        + " Cells: "
        + str(self.cellCount())
        + " Boundaries: "
        + str(self.boundaryCount())
    )

    if self.secondaryNodeCount() > 0:
        st += " secNodes: " + str(self.secondaryNodeCount())

    if len(list(self.dataMap().keys())) > 0:
        st += "\nMesh contains data: "

        uniqueNames = self.dataKeys()

        for d, v in uniqueNames.items():
            if len(v) > 1:
                st += d + "[0,...,{})".format(len(v))
            else:
                st += d

            st += ", "

        st = st.rstrip(", ")

    return st


Mesh.__repr__ = __Mesh_str


def __addPLCs__(self, other):
    if isR3Array(other):
        m = Mesh(dim=self.dim(), isGeometry=True)
        [m.createNode(n) for n in other]
        other = m

    if self.isGeometry() and other.isGeometry():
        return mergePLC([self, other])
    else:
        error("Addition is only supported for PLCs, i.e. meshs without cells.")


Mesh.__add__ = __addPLCs__


def __MeshEntity_str(self):
    """Give mesh entity infos."""
    s = self.__repr__()
    s += (
        "\tID: "
        + str(self.id())
        + ", Marker: "
        + str(self.marker())
        + ", Size: "
        + str(self.size())
        + "\n"
    )

    if isinstance(self, PolygonFace) and len(self.nodes()) > 5:
        s += "\t" + str(self.nodeCount()) + " Nodes.\n"
    else:
        for n in self.nodes():
            s += "\t" + str(n.id()) + " " + str(n.pos()) + "\n"
    return s


MeshEntity.__str__ = __MeshEntity_str


def __Node_str(self):
    """Give node infos."""
    s = "\tID: " + str(self.id()) + ", Marker: " + str(self.marker())
    s += "\t" + str(self.pos()) + "\n"
    return s


Node.__repr__ = __Node_str


def __Mesh_setVal(self, key, val):
    """Index access to the mesh data.

    Multiple arrays via matrix will be saved too.
    """
    if isinstance(val, (float, int)):  # single value
        val = RVector(self.cellCount(), val)

    if isR3Array(val):  # vectorial property
        return self.addData(key, val)
    # list of vectors
    if (
        isinstance(val, list)
        and isinstance(val[0], (RVector, np.ndarray))
        or val.ndim == 2
        or val.ndim == 3
    ):

        maxDigit = ceil(np.log10(len(val)))

        for i, v in enumerate(val):
            self.addData("{}#{}".format(key, str(i).zfill(maxDigit)), np.asarray(v))
    else:
        try:
            if key == "marker":
                self.setCellMarkers(val)
            else:
                self.addData(key, val)
        except:
            ValueError(
                "Could not add data with key " + key + " of shape " + str(np.shape(val))
            )


Mesh.__setitem__ = __Mesh_setVal


def __Mesh_getVal(self, key):
    """Index access to the mesh data"""
    if self.haveData(key):
        return self.data(key)
    else:

        uniqueNames = {}
        for d in self.dataMap().keys():
            if "_y" in d or "_z" in d:
                continue

            uName = d

            if "_x" in uName:
                uName = uName[0 : uName.find("_x")]
                d1 = self.data(d)
                d2 = self.data(d.replace("_x", "_y"))
                d3 = self.data(d.replace("_x", "_z"))
                dat = np.array([d1, d2, d3]).T
            else:
                dat = self.data(d)

            if "#" in uName:
                uName = uName[0 : uName.find("#")]

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

        critical(
            NameError,
            "The mesh does not have the requested data:",
            key,
            ". Available:",
            uniqueNames,
        )


Mesh.__getitem__ = __Mesh_getVal


def __MeshBoundingBox__(self):
    bb = self.boundingBox()
    mi = RVector3([bb.min()[i] for i in range(3)])
    ma = RVector3([bb.max()[i] for i in range(3)])
    return [mi, ma]


Mesh.bb = __MeshBoundingBox__


def __MeshGetCellMarker__(self):
    deprecated(msg="Mesh::cellMarker()", hint="Mesh::cellMarkers()")
    return self.cellMarkers()


def __MeshSetCellMarker__(self, m):
    deprecated(msg="Mesh::setCellMarker()", hint="Mesh::setCellMarkers()")
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
                            pos = bs.xyz([(sx + 1) / (n + 2), (sy + 1) / (n + 2)])
                        else:
                            pos = bs.xyz([(sx + 1) / (n + 1), (sy + 1) / (n + 1)])

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
                    warn("cell type unknown")

                for e in edges:
                    line = Line(e[0].pos(), e[1].pos())
                    for i in range(n):
                        sn = self.createSecondaryNode(
                            line.at((i + 1) / (n + 1)), tol=1e-6
                        )
                        c.addSecondaryNode(sn)
        else:
            warn("Unknown dimension. Don't know what to do.")

    if verbose:
        info("Added %d secondary nodes." % self.secondaryNodeCount())


Mesh.createSecondaryNodes = __createSecondaryNodes__


def __createMeshWithSecondaryNodes__(self, n=3, verbose=False):
    m = Mesh(self)
    m.createSecondaryNodes(n, verbose)
    return m


Mesh.createMeshWithSecondaryNodes = __createMeshWithSecondaryNodes__

__Mesh_deform__ = Mesh.deform


def __deform__(self, eps, mag=1.0):
    v = None
    dof = self.nodeCount()
    if hasattr(eps, "ndim") and eps.ndim == 1:
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
            error("Size of displacement does not match mesh nodes size.")
    elif len(eps) == self.nodeCount() and eps.ndim == 2:
        v = eps.reshape(self.nodeCount() * eps.shape[1], order="F")

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


def __Mesh_findPaths__(self, bounds):
    """Find paths of connected boundaries

    Returns
    -------
    List of list of ids of connected nodes
    """
    import pygimli as pg

    scipy = pg.optImport("scipy")
    scipy.sparse = pg.optImport("scipy.sparse")

    S = pg.core.SparseMapMatrix()
    for b in bounds:
        # S[b.shape().node(0).id(), b.shape().node(1).id()] = 1
        # S[b.shape().node(1).id(), b.shape().node(0).id()] = 1
        S.addVal(b.shape().node(1).id(), b.shape().node(0).id(), 2.0)
        S.addVal(b.shape().node(0).id(), b.shape().node(1).id(), 1.0)

    S = scipy.sparse.dok_matrix(pg.utils.toCOO(S))

    paths = []

    def followPath(path, S, rID):
        # print('start', rID)
        row = S[rID]
        while 1:
            cID = list(row.keys())[0][1]
            # print('row', rID, 'col', cID)
            # print('add', cID)
            path.append(cID)
            # try:
            # print('pop-r', (rID, cID))
            del S[(rID, cID)]
            del S[(cID, rID)]

            col = S[:, cID]

            if len(col) == 1:
                rID = list(col.keys())[0][0]
                path.append(rID)
                # print('add', rID)
                # print('pop-c', (rID, cID))
                del S[(rID, cID)]
                del S[(cID, rID)]
                # print('pop-s', (rID, cID))

                row = S[rID]
                if len(row) != 1:
                    break
            else:
                break

    ## first look for single starting
    for i in range(S.shape[0]):
        rID = i
        row = S[rID]

        if len(row) == 1:
            # single starting
            path = []
            paths.append(path)
            # starting node
            path.append(rID)
            followPath(path, S, rID)

    ## remaining are closed
    for i in range(S.shape[0]):
        rID = i
        row = S[rID]

        if len(row) == 2:
            path = []
            paths.append(path)
            # starting node
            path.append(rID)
            followPath(path, S, rID)

    return paths


Mesh.findPaths = __Mesh_findPaths__


def __Mesh_cutBoundary__(self, marker, boundaryMarker=None):
    """Cut the mesh along a given inner boundary.

    Cut the mesh along a given boundary and convert this inner boundary to an outer. There will be new nodes to cut the connection between neighbouring cells. The new boundary can have an optional boundaryMarker.

    Restrictions
    ------------
        * 2D p1
        * one connected path at once
        * starting node need to be on an outer boundary
        * end node needs to be inside the mesh

    TODO
    ----
        * remove restrictions

    Arguments
    ---------
    mesh: :gimliapi:`GIMLI::Mesh`
        2D
    marker: int
        Marker for the boundary to be cut.
    boundaryMarker: None
        If set to None, boundaryMarker set to marker.

    Example
    -------
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> plc = mt.createCircle(nSegments=24)
    >>> l = mt.createLine(start=[0, -1], end=[0, -0.1], boundaryMarker=2)
    >>> mesh = mt.createMesh([plc, l], area=0.1, quality=30)
    >>> fig, axs= pg.plt.subplots(1, 2)
    >>> ax ,_ = pg.show(mesh, boundaryMarkers=True, ax=axs[0])
    >>> oldNodeCount = mesh.nodeCount()
    >>> print(mesh)
    Mesh: Nodes: 43 Cells: 60 Boundaries: 102
    >>> mesh.cutBoundary(marker=2, boundaryMarker=3)
    >>> print(mesh)
    Mesh: Nodes: 46 Cells: 60 Boundaries: 105
    >>> ## just move the new nodes little rightwards to see the cut
    >>> for n in range(oldNodeCount, mesh.nodeCount()):
    ...     mesh.node(n).setPos(mesh.node(n).pos() + [0.1, 0.0])
    >>> ax, _ = pg.show(mesh, data=range(mesh.cellCount()),
    ...                 boundaryMarkers=True, colorBar=False,
    ...                 showMesh=True, boundaryProps={'lw':2}, ax=axs[1])
    >>> for b in mesh.boundaries():
    ...     if b.marker() != 0:
    ...         c = b.center()
    ...         n = b.norm()
    ...         _ = ax.annotate('', xytext=(c[0], c[1]),
    ...                    xy=((c+n/30.)[0], (c+n/30.)[1]),
    ...                    arrowprops=dict(arrowstyle="-|>", lw=1),
    ...                     )
    """
    import pygimli as pg

    if boundaryMarker is None:
        boundaryMarker = marker

    mesh = self

    def replaceNode_(mesh, c, n1, n2, marker, lastC=None):
        if c is None or n1.id() not in c.ids():
            return
        # pg._y('check in cell', c.id(), n1.id(), n2.id())

        toBeReplaced = []
        for i in range(c.boundaryCount()):
            b = c.boundary(i)
            if b is not None and n1.id() in b.ids():
                # pg._y('\tbound: ', b.id(), ':',
                #         b.node(0).id(), b.node(1).id(), "ma:", b.marker())

                if b.marker() != marker:
                    lC = b.leftCell()
                    rC = b.rightCell()
                    # rcS = None
                    # if rC is not None:
                    #     rcS = rC.id()
                    # lcS = None
                    # if lC is not None:
                    #     lcS = lC.id()

                    # pg._y('\tNeigh: {0} : {1}'.format(lcS, rcS))
                    # pg._r("add:", b.id())
                    toBeReplaced.append(b)

                    if lC != c and lC is not None and lC != lastC:
                        # pg._r("follow up left:", lC.id())
                        replaceNode_(mesh, lC, n1, n2, marker, c)

                    if rC != c and rC is not None and rC != lastC:
                        # pg._r("follow up right:", rC.id())
                        replaceNode_(mesh, rC, n1, n2, marker, c)

        for b in toBeReplaced:
            nIds = [n2.id() if n == n1.id() else n for n in b.ids()]
            # pg._r('replace in boundary', b.id(), ':', n1.id(), n2.id())
            b.setNodes(mesh.nodes(nIds))

        # pg._r('replace in cell:', c.id(), ':', n1.id(), n2.id())
        nIds = [n2.id() if n == n1.id() else n for n in c.ids()]
        c.setNodes(mesh.nodes(nIds))

    paths = mesh.findPaths(mesh.findBoundaryByMarker(marker))
    if len(paths) == 0:
        pg.error("did not found path for marker: {0}".format(marker))

    newNodes = []

    if len(paths[0]) == 0:
        pg.error("did not found path for marker: {0}".format(marker))

    ## step 1 . fix direction along the path
    for i in range(len(paths[0]) - 1):
        nA1 = mesh.node(paths[0][i])
        nB1 = mesh.node(paths[0][i + 1])
        b = pg.core.findBoundary(nA1, nB1)

        if b.node(0) != nA1:
            b.swapNorm()

        lC = b.leftCell()
        rC = b.rightCell()
        if rC is None or lC is None:
            pg.error("Path is not inside the mesh")
            return

    ## add new nodes and decouple cells along the path
    rightCells = []
    for i in range(len(paths[0]) - 1):

        nA1 = mesh.node(paths[0][i])
        nB1 = mesh.node(paths[0][i + 1])
        b = pg.core.findBoundary(nA1, nB1)

        lC = b.leftCell()
        rC = b.rightCell()

        # pg._y(b.node(0).id(), b.node(1).id(), 'N', nA1.id(), nB1.id(), ':', lC.id(), rC.id())

        ### only if on outer boundary .. need check!!
        nA2 = mesh.createNode(nA1.pos(), nA1.marker())
        newNodes.append(nA2)

        # nA2 = newNodes[-1]
        if rC is not None:
            b.setRightCell(None)
            rightCells.append(rC)
            replaceNode_(mesh, rC, nA1, nA2, marker=b.marker())

    newNodes.append(mesh.node(paths[0][-1]))
    for i in range(len(newNodes) - 1):
        b = mesh.createBoundary(
            [newNodes[i + 1].id(), newNodes[i].id()], marker=boundaryMarker
        )
        b.setLeftCell(rightCells[i])


Mesh.cutBoundary = __Mesh_cutBoundary__


def __Mesh__align__(self, pnts):
    """Align 2D mesh along 3D coordinates.

    Align a xy-mesh along xyz-coordinates. x and y coordinates of the 2D mesh will be interpolated to x and y of pnts, where depth y from the mesh will become z and preserves its values.

    TODO
    ....
        * handle z coordinate if pnts contain z

    Args
    ....
    mesh: :gimliapi:`GIMLI::Mesh`
        2D mesh, assumed to be aligned along x-axis. Depth is y-axis.
    pnts: [[x,y],] | [[dx, x, y],]
        * `shape[1] == 2`: Points that will be interpreted as xyz coordinates.
        * `shape[1] == 3`: interpreted as dx, x, y. Dx should start with <=0 max dx should be larger than `mesh.xmax() - mesh.xmin()`
    """
    if self.dim() != 2:
        pg.critical("Only 2D meshes can be aligned to 3D coordinates")

    A = None

    pnts = np.asarray(pnts)
    if pnts.ndim == 2:
        if pnts.shape[1] == 2:

            A = np.zeros((pnts.shape[0], 3))

            from ..utils import cumDist

            A[:, 0] = cumDist(pnts[:, 0:2])
            A[:, 1] = pnts[:, 0]
            A[:, 2] = pnts[:, 1]

        elif pnts.shape[1] == 3:
            A = pnts

    if A is None:
        print(pnts)
        pg.critical("Can't, interprete ptns.")

    tn = [n.pos()[0] for n in self.nodes()]
    zn = [n.pos()[1] for n in self.nodes()]

    p = interpolateAlongCurve(A[:, 1:3], tn, tCurve=A[:, 0])

    for i, n in enumerate(self.nodes()):
        n.setPos((p[i][0], p[i][1], zn[i]))

    self.geometryChanged()


Mesh.align = __Mesh__align__


def __Mesh__swapOrientation__(self):
    """Swap orientation from one right-hand-side type into another.

    - xyz (right, up, up) with z pointing upwards
    - NED (North, East, Down) with z pointing downwards
    """
    self.swapCoordinates(0, 1)  # exchange x and y
    self.scale([1, 1, -1])  # revert z


Mesh.swapOrientation = __Mesh__swapOrientation__


def __Mesh__copy(self):
    """Return copy of a mesh."""
    return Mesh(self)


Mesh.copy = __Mesh__copy


def __Mesh__NED__(self):
    """Return NED copy of mesh."""
    newmesh = self.copy()
    newmesh.swapOrientation()
    return newmesh


Mesh.NED = __Mesh__NED__


def __Mesh__midpoint__(self):
    """Return midpoint."""
    return sum(self.bb()) / 2


Mesh.midpoint = __Mesh__midpoint__


def __Mesh__extent__(self, axis=None):
    """Return extent of mesh.

    Parameters
    ----------
    axis : str|int
        axis along to measure the extent
        0|'x' : x direction
        1|'y' : y direction
        2|'z' : z direction
        None|'max' : maximum of x, y, z
        -1|'d' : diagonal

    Example
    -------
    >>> import numpy as np
    >>> import pygimli.meshtools as mt
    >>> m = mt.createGrid(3, 4, 2)
    >>> print(m.extent())
    3.0
    >>> print(m.extent('x'))
    2.0
    >>> print(m.extent(2))
    1.0
    >>> print(np.round(m.extent('d'), 2))
    3.74
    """
    bb = self.bb()
    dist = bb[1] - bb[0]
    if isinstance(axis, str):
        sa0 = axis.lower()[0]
        if sa0 == "m":
            axis = None
        else:
            axis = "dxyz".find(sa0) - 1

    if axis is None:
        return max(np.abs(dist))
    elif axis < 0:
        return dist.abs()
    else:
        return abs(dist[axis])


Mesh.extent = __Mesh__extent__


def __Mesh__populate__(self, prop: str, value):
    """Fill property of mesh with values from map or vector.

    Parameters
    ==========
    prop : str
        property to be filled, i.e. mesh[prop]
    value : array|dict|map
        values to be sorted into field
        * full length (cellCount) array
        * array indexing into cell markers
        * map like [[1, 100], [2, 1000]]
        * dict like {1: 100, 2: 1000}
    Return
    ======
    field : array (length cellCount)
    """
    from pygimli.solver import parseMapToCellArray

    if isinstance(value, dict):
        self[prop] = parseMapToCellArray(value, self)
    elif hasattr(value, "__iter__"):
        if hasattr(value[0], "__iter__"):
            self[prop] = parseMapToCellArray(value, self)
        elif len(value) == self.cellCount():
            self[prop] = value
        else:
            raise Exception("Length mismatch!")

    return self[prop]


Mesh.populate = __Mesh__populate__


def __Mesh__submesh__(self, relation="=", **kwargs):
    """Select submesh according to keywords.

    Parameters
    ----------
    relation : str
        logical operation ["=", ">", "<"]
    kwargs
        key:val pairs satisfying self[key] OP val
    """
    istrue = np.ones(self.cellCount(), dtype=bool)
    for key, val in kwargs.items():
        if key not in self.dataKeys():
            raise IndexError("Property not in mesh: ", key)
        if relation == ">":
            isk = self[key] > val
        elif relation == "<":
            isk = self[key] < val
        elif relation == "!=":
            isk = self[key] != val
        else:
            isk = self[key] == val

        istrue = np.bitwise_and(istrue, isk)

    return self.createSubMesh(self.cells(np.nonzero(istrue)[0]))


Mesh.submesh = __Mesh__submesh__


def __Mesh__innerBoundaryCenters__(self):
    """Center of all inner boundaries (C1-constraints)."""
    return [b.center() for b in self.boundaries() if not b.outside()]


Mesh.innerBoundaryCenters = __Mesh__innerBoundaryCenters__
