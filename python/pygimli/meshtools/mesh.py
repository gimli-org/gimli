# -*- coding: utf-8 -*-
"""General mesh generation and maintenance."""

import os

import numpy as np

import pygimli as pg


def createMesh(poly, quality=32, area=0.0, smooth=None, switches=None,
               verbose=False, **kwargs):
    """Create a mesh for a given PLC or point list.

    The mesh is created by :term:`triangle` or :term:`tetgen` if the
    pgGIMLi support for these mesh generators is installed.
    A PLC needs to contain nodes and boundaries and should be valid
    in the sense that the boundaries are non-intersecting.

    If poly is a list of coordinates, a simple Delaunay mesh with a convex hull
    will be created. Quality and area arguments are ignored for this case to
    create a mesh with one node for each coordinate position.

    Parameters
    ----------
    poly: :gimliapi:`GIMLI::Mesh` or list or ndarray
        * 2D or 3D gimli mesh that contains the PLC.
        * 2D mesh needs edges
        * 3D mesh needs a plc and tetgen as system component
        * List of x y pairs [[x0, y0], ... ,[xN, yN]]
        * ndarray [x_i, y_i]
        * PLC or list of PLCs

    quality: float
        2D triangle quality sets a minimum angle constraint.
        Be careful with values above 34 degrees.
        3D tetgen quality. Be careful with values below 1.12.

    area: float
        Maximum element size (global). 2D maximum triangle size in m*²,
        3D maximum tetrahedral size in m³.

    smooth: tuple
        [smoothing algorithm, number of iterations]
        0: no smoothing
        1: node center
        2: weighted node center

    switches: str
        Set additional triangle command switches.
        https://www.cs.cmu.edu/~quake/triangle.switch.html

    Returns
    -------
    mesh: :gimliapi:`GIMLI::Mesh`

    Examples
    --------
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> rect = mt.createRectangle(start=[0, 0], end=[4, 1])
    >>> ax, _ = pg.show(mt.createMesh(rect, quality=10))
    >>> ax, _ = pg.show(mt.createMesh(rect, quality=33))
    >>> ax, _ = pg.show(mt.createMesh(rect, quality=33, area=0.01))
    >>> pg.wait()
    """
    #  poly == [pg.Mesh, ]
    if isinstance(poly, list):
        if isinstance(poly[0], pg.Mesh):
            return createMesh(
                pg.meshtools.mergePLC(poly), quality, area, smooth, switches,
                verbose)
    # poly == [pos, pos, ]
    if isinstance(poly, list) or \
        isinstance(poly, type(zip)) or \
        isinstance(poly, pg.core.stdVectorRVector3) or \
        isinstance(poly, pg.core.R3Vector) or \
            (isinstance(poly, np.ndarray) and poly.ndim == 2):
        delPLC = pg.Mesh(2)
        for p in poly:
            delPLC.createNode(p[0], p[1], 0.0)
        return createMesh(delPLC, switches='-zeY')

    # poly == Mesh
    if poly.dim() == 2:
        if poly.nodeCount() == 0:
            pg.critical("No nodes in poly to create a valid mesh")

        if switches is None:
            # -D Conforming delaunay
            # -F Uses Steven Fortune's sweepline algorithm
            switches = 'pzeA'

            if area > 0:
                # The str function turns everything smaller
                # than 0.0001 into the scientific notation 1e-5
                # which can not be read by triangle. The following
                # avoids this even for very small numbers
                switches += 'a' + '{:.20f}'.format(area)

            # switches = switches.replace('.', ',')
            switches += 'q' + str(quality)

            # an EXTRA! -a here else it ignores per region area
            switches += 'a'

        if not verbose:
            switches += 'Q'

        pg.verbose(switches)

        tri = pg.core.TriangleWrapper(poly)
        tri.setSwitches(switches)
        mesh = tri.generate()

        if smooth is not None:
            mesh.smooth(nodeMoving=kwargs.pop('node_move', True),
                        edgeSwapping=False,
                        smoothFunction=smooth[0],
                        smoothIteration=smooth[1])
        return mesh

    else:
        # 3d case
        if quality == 32:
            quality = 1.2

        tmp = pg.optImport('tempfile')
        _, namePLC = tmp.mkstemp(suffix='.poly')

        pg.meshtools.exportPLC(poly, namePLC)
        mesh = pg.meshtools.syscallTetgen(namePLC, quality, area,
                                          verbose=verbose, **kwargs)

        try:
            os.remove(namePLC)
        except BaseException as e:
            print(e)
            print("can't remove:", namePLC)


        return mesh


def refineQuad2Tri(mesh, style=1):
    """Refine mesh of quadrangles into a mesh of triangle cells.

        TODO mixed meshes

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        Mesh containing quadrangle cells.

    style: int [1]
        * 1 bisect each quadrangle into 2 triangles
        * 2 cross-sect each quadrangle into 4 triangles

    Returns
    -------
    ret : :gimliapi:`GIMLI::Mesh`
        Mesh containing triangle cells.

    Examples
    --------
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> quads = pg.createGrid(range(10), range(10))
    >>> ax, _ = pg.show(quads)
    >>> ax, _ = pg.show(mt.refineQuad2Tri(quads, style=1))
    >>> ax, _ = pg.show(mt.refineQuad2Tri(quads, style=2))
    >>> pg.wait()
    """
    out = pg.Mesh(2)
    newNode = None

    for n in mesh.nodes():
        out.createNode(n.pos())

    for c in mesh.cells():

        if style == 1:
            out.createCell([c.node(0).id(), c.node(1).id(), c.node(2).id()],
                           c.marker())
            out.createCell([c.node(0).id(), c.node(2).id(), c.node(3).id()],
                           c.marker())

        elif style == 2:
            newNode = out.createNodeWithCheck(c.center())

            for i in range(4):
                out.createCell([c.node(i).id(), c.node((i + 1) % 4).id(),
                                newNode.id()], c.marker())

        for i in range(c.boundaryCount()):
            b = c.boundary(i)
            if b.marker() != 0:
                out.createBoundary([b.node(0).id(), b.node(1).id()],
                                   b.marker())

    out.createNeighborInfos()

    return out


def refineHex2Tet(mesh, style=1):
    """Refine mesh of hexahedra into a mesh of tetrahedra.

    TODO
    ----
        * mixed meshes
        * preserve boundary information

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        Mesh containing hexrahedron cells, e.g., from a grid.

    style: int [1]
        * 1 bisect each hexahedron int 5 tetrahedrons (numerically favorable)
        * 2 bisect each hexahedron int 6 tetrahedrons

    Returns
    -------
    ret : :gimliapi:`GIMLI::Mesh`
        Mesh containing tetrahedrons cells.

    Examples
    --------
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> hex = pg.createGrid(2, 2, 2)
    >>> print(hex)
    Mesh: Nodes: 8 Cells: 1 Boundaries: 6
    >>> tet = mt.refineHex2Tet(hex, style=1)
    >>> print(tet)
    Mesh: Nodes: 8 Cells: 5 Boundaries: 0
    >>> tet = mt.refineHex2Tet(hex, style=2)
    >>> print(tet)
    Mesh: Nodes: 8 Cells: 6 Boundaries: 0
    """
    out = pg.Mesh(3)

    for n in mesh.nodes():
        out.createNode(n.pos())

    HexahedronSplit5TetID = [[1, 4, 5, 6],
                             [3, 6, 7, 4],
                             [1, 0, 4, 3],
                             [1, 2, 3, 6],
                             [1, 4, 6, 3]]

    HexahedronSplit6TetID = [[0, 1, 2, 6],
                             [0, 2, 3, 6],
                             [0, 1, 6, 5],
                             [0, 4, 5, 6],
                             [0, 3, 7, 6],
                             [0, 4, 6, 7]]

    for c in mesh.cells():
        if style == 1:
            for tet in HexahedronSplit5TetID:
                out.createCell([c.node(tet[0]).id(),
                                c.node(tet[1]).id(),
                                c.node(tet[2]).id(),
                                c.node(tet[3]).id()], c.marker())
        elif style == 2:
            for tet in HexahedronSplit6TetID:
                out.createCell([c.node(tet[0]).id(),
                                c.node(tet[1]).id(),
                                c.node(tet[2]).id(),
                                c.node(tet[3]).id()], c.marker())
    return out


def extrudeMesh(mesh, a, **kwargs):
    r"""Extrude mesh to a higher dimension.

    Generates a 2D mesh by extruding a 1D mesh along y-coordinate using quads.
    We assume a 2D mesh here consisting of nodes and edges.
    The marker of nodes are extruded as edges with the same marker.
    The marker of the edges are extruded as cells with same marker.
    Optionally all y-coordinates can be adjusted to become equal at the end

    Generates a three-dimensional mesh by extruding a two-dimensional mesh
    along the z-coordinate transforming triangles into triangular prisms or
    quads into hexahedrons.
    3D cell markers are set from 2D cell marker.
    The boundary marker for the side boundaries are set from edge markers.

    TODO
    ----
        * document and test marker setting from the core

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`
        Input mesh
    a: iterable (float)
        Additional coordinate to extrude into.

    Additional Parameters
    -------------------
    adjustBottom: bool [False]
        Adjust all nodes that bottom of the mesh has a constant depth (only 2D)

    Returns
    -------
    mesh: :gimliapi:`GIMLI::Mesh`
        Returning mesh of +1 dimension

    Examples
    --------
    >>> import numpy as np
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> topo = [[x, 1.0+np.cos(2*np.pi*1/30*x)] for x in range(31)]
    >>> m1 = mt.createPolygon(topo)
    >>> m1.setBoundaryMarkers(range(m1.boundaryCount()))

    >>> m = mt.extrudeMesh(m1, a=-(np.geomspace(1, 5, 8)-1.0))
    >>> _ = pg.show(m, m.cellMarkers(), showMesh=True)
    >>> m = mt.extrudeMesh(m1, a=-(np.geomspace(1, 5, 8)-1.0),
    ...                    adjustBottom=True)
    >>> _ = pg.show(m, m.cellMarkers(), showMesh=True)
    """
    # case 1d or edge list to 2d grid
    if mesh.dim() == 1 or (mesh.dim() == 2 and mesh.cellCount() == 0):

        adjustBack = kwargs.pop('adjustBottom', False)
        m2 = pg.meshtools.createMesh2D(mesh, y=a, **kwargs)
        if adjustBack:
            minY = min(pg.y(mesh)) + min(a)
            scale = dict()
            n0s = dict()
            for n in mesh.nodes():
                n0s[n.pos()[0]] = n.pos()
                scale[n.pos()[0]] = (minY-n.pos()[1]) / min(a)
            for n in m2.nodes():
                xP = n.pos()[0]
                yP = n.pos()[1]
                y0 = n0s[xP][1]
                n.setPos([xP, (yP-y0)*scale[xP] + y0])
        return m2

    if mesh.dim() == 2:
        return pg.meshtools.createMesh3D(mesh, z=a, **kwargs)
    pg.error('Cannot extrude mesh of dimension:', mesh.dim())


def readGmsh(fName, verbose=False, precision=None):
    r"""Read :term:`Gmsh` ASCII file and return instance of GIMLI::Mesh class.

    Parameters
    ----------
    fName : string
        Filename of the file to read (\\*.msh). The file must conform
        to the `MSH ASCII file version 2
        <http://gmsh.info/doc/texinfo/gmsh.html#MSH-ASCII-file-format>`_ format
    verbose : boolean, optional
        Be verbose during import. Default: False
    precision : None|int, optional
        If not None, then round off node coordinates to the provided number of
        digits using numpy.round. This is useful in case that nodes are
        accessed using their coordinates, in which case numerical discrepancies
        can occur.

    Notes
    -----
    Physical groups specified in Gmsh are interpreted as follows:

    - Points with the physical number 99 are interpreted as sensors. Note that
      physical point groups are ordered with respect to the node tag. For
      example, "Physical Point (99) = {50, 34};" and "Physical Point (99) = {34,
      50};" will yield the same mesh. This must be taken into account when
      defining measurement configurations using electrodes defined in GMSH
      using marker 99.
    - ERT only: Points with markers 999 and 1000 are used to mark calibration
      and reference nodes.
    - Physical Lines and Surfaces define boundaries in 2D and 3D, respectively.
        - Physical Number 1: Homogeneous Neumann condition
        - Physical Number 2: Mixed boundary condition
        - Physical Number 3: Homogeneous Dirichlet condition
        - Physical Number 4: Dirichlet condition
    - Physical Surfaces and Volumes define regions in 2D and 3D, respectively.
        - Physical Number 1: No inversion region
        - Physical Number >= 2: Inversion region

    Examples
    --------
    >>> import tempfile, os
    >>> from pygimli.meshtools import readGmsh
    >>> gmsh = '''
    ... $MeshFormat
    ... 2.2 0 8
    ... $EndMeshFormat
    ... $Nodes
    ... 3
    ... 1 0 0 0
    ... 2 0 1 0
    ... 3 1 1 0
    ... $EndNodes
    ... $Elements
    ... 7
    ... 1 15 2 0 1 1
    ... 2 15 2 0 2 2
    ... 3 15 2 0 3 3
    ... 4 1 2 0 1 2 3
    ... 5 1 2 0 2 3 1
    ... 6 1 2 0 3 1 2
    ... 7 2 2 0 5 1 2 3
    ... $EndElements
    ... '''
    >>> fName = tempfile.mktemp()
    >>> with open(fName, "w") as f:
    ...     f.writelines(gmsh)
    >>> mesh = readGmsh(fName)
    >>> print(mesh)
    Mesh: Nodes: 3 Cells: 1 Boundaries: 3
    >>> os.remove(fName)
    """
    assert precision is None or precision >= 0
    inNodes, inElements, ncount = 0, 0, 0
    fid = open(fName)
    if verbose:
        print('Reading %s... \n' % fName)

    for line in fid:

        if line[0] == '$':
            if line.find('Nodes') > 0:
                inNodes = 1
            if line.find('EndNodes') > 0:
                inNodes = 0
            if line.find('Elements') > 0:
                inElements = 1
            if line.find('EndElements') > 0:
                inElements = 0

        else:
            if inNodes == 1:
                if len(line.split()) == 1:
                    nodes = np.zeros((int(line), 3))
                    if verbose:
                        print('  Nodes: %s' % int(line))
                else:
                    node_coordinates = np.array(line.split(), 'float')[1:]
                    if precision is None:
                        nodes[ncount, :] = node_coordinates
                    else:
                        nodes[ncount, :] = np.round(
                            node_coordinates, precision)
                    ncount += 1

            elif inElements == 1:
                if len(line.split()) == 1:
                    if verbose:
                        print('  Entries: %s' % int(line))
                    points, lines, triangles, tets = [], [], [], []

                else:
                    # Element entries follow the following format:
                    # elm-number elm-type number-of-tags < tag > …
                    #   node-number-list

                    # strip elm-number here
                    entry = [int(e_) for e_ in line.split()][1:]

                    if entry[0] == 15:  # Points
                        # point node, marker (1st tag)
                        points.append((entry[-1], entry[2]))
                    elif entry[0] == 1:
                        lines.append((entry[-2], entry[-1], entry[2]))
                    elif entry[0] == 2:
                        triangles.append((entry[-3], entry[-2], entry[-1],
                                          entry[2]))
                    elif entry[0] == 4:
                        tets.append((entry[-4], entry[-3], entry[-2],
                                     entry[-1], entry[2]))
                    elif entry[0] in [3, 6]:
                        pg.error(
                            "Quadrangles and prisms are not supported yet")

    fid.close()
    lines = np.asarray(lines)
    triangles = np.asarray(triangles)
    tets = np.asarray(tets)

    if verbose:
        print('    Points: %s' % len(points))
        print('    Lines: %s' % len(lines))
        print('    Triangles: %s' % len(triangles))
        print('    Tetrahedra: %s \n' % len(tets))
        print('Creating mesh object... \n')

    # check dimension
    if len(tets) == 0:
        dim, bounds, cells = 2, lines, triangles
        zero_dim = np.abs(nodes.sum(0)).argmin()  # identify zero dimension
    else:
        dim, bounds, cells = 3, triangles, tets
    if verbose:
        print('  Dimension: %s-D' % dim)

    # creating instance of GIMLI::Mesh class
    mesh = pg.Mesh(dim)

    # replacing boundary markers (gmsh does not allow negative phys. regions)
    bound_marker = (pg.core.MARKER_BOUND_HOMOGEN_NEUMANN,
                    pg.core.MARKER_BOUND_MIXED,
                    pg.core.MARKER_BOUND_HOMOGEN_DIRICHLET,
                    pg.core.MARKER_BOUND_DIRICHLET)

    if bounds.any():
        for i in range(4):
            bounds[:, dim][bounds[:, dim] == i + 1] = bound_marker[i]

        # account for CEM markers
        bounds[:, dim][bounds[:, dim] >= 10000] *= -1

        if verbose:
            bound_types = np.unique(bounds[:, dim])
            print('  Boundary types: %s ' % len(bound_types) + str(
                tuple(bound_types)))
    else:
        print("WARNING: No boundary conditions found.",
              "Setting Neumann on the outer edges by default.")

    if verbose:
        regions = np.unique(cells[:, dim + 1])
        print('  Regions: %s ' % len(regions) + str(tuple(regions)))

    for node in nodes:
        if dim == 2:
            mesh.createNode(node[0], node[3 - zero_dim], 0)
        else:
            mesh.createNode(node)

    for cell in cells:
        if dim == 2:
            mesh.createTriangle(
                mesh.node(int(cell[0] - 1)), mesh.node(int(cell[1] - 1)),
                mesh.node(int(cell[2] - 1)), marker=int(cell[3]))
        else:
            mesh.createTetrahedron(
                mesh.node(int(cell[0] - 1)), mesh.node(int(cell[1] - 1)),
                mesh.node(int(cell[2] - 1)), mesh.node(int(cell[3] - 1)),
                marker=int(cell[4]))

    mesh.createNeighborInfos()

    # Set Neumann on outer edges by default (can be overwritten by Gmsh info)
    for b in mesh.boundaries():
        if not b.leftCell() or not b.rightCell():
            b.setMarker(pg.core.MARKER_BOUND_HOMOGEN_NEUMANN)

    for bound in bounds:
        if dim == 2:
            mesh.createEdge(
                mesh.node(int(bound[0] - 1)), mesh.node(int(bound[1] - 1)),
                marker=int(bound[2]))
        else:
            mesh.createTriangleFace(
                mesh.node(int(bound[0] - 1)), mesh.node(int(bound[1] - 1)),
                mesh.node(int(bound[2] - 1)), marker=int(bound[3]))

    # assign marker to corresponding nodes (sensors, reference nodes, etc.)
    if points:
        for point in points:
            mesh.node(point[0] - 1).setMarker(-point[1])

    if verbose:
        if points:
            points = np.asarray(points)
            node_types = np.unique(points[:, 1])
            print('  Marked nodes: %s ' % len(points) + str(tuple(node_types)))
        print('\nDone. \n')
        print('  ' + str(mesh))
    return mesh


def readTriangle(fName, verbose=False):
    r"""Read :term:`Triangle` :cite:`Shewchuk96b` mesh.

    Read :term:`Triangle` :cite:`Shewchuk96b` ASCII mesh files and return an
    instance of GIMLI::Mesh class.
    See: ://www.cs.cmu.edu/~quake/triangle.html

    Parameters
    ----------
    fName : string
        Filename of the file to read (\\*.n, \\*.e)

    verbose : boolean, optional
        Be verbose during import.

    """
    raise Exception("implement me!" + fName + str(verbose))
    # os.system('meshconvert -d2 ' + fName)
    # return pg.Mesh(2)


def readTetgen(fName, comment='#', verbose=False, defaultCellMarker=0,
               loadFaces=True, quadratic=False):
    """
    Read and convert a mesh from the basic :term:`Tetgen` output.

    Read :term:`Tetgen` :cite:`Si2004` ASCII files and return instance
    of :gimliapi:`GIMLI::Mesh` class.
    See: http://tetgen.org/

    Parameters
    ----------

    fName: str
        Base name of the tetgen output, without ending. All additional files
        (.node, .ele and .face) need to have the same basename.

    comment: str ('#')
        String consisting of all symbols indicating a comment in the input
        files. Standard for tetgen files is the '#'.

    verbose: boolean (True)
        Enables console output during the import process.

    defaultCellMarker: int (0)
        :term:`Tetgen` files can contain cell markers, but do not have to.
        If no markers are found, the given integer is used.

    loadFaces:
        Optional decision whether the faces of :term:`Tetgen` output (.face)
        are loaded or not. Note that without the -f in during the tetgen call,
        the faces in the .face file will only contain the faces of the original
        input poly file and not all faces. If only a part of the faces are
        imported, a createNeighborInfos call of the mesh will fail.

    quadratic: boolean (False)
        Returns a P2 (quadratic) refined mesh when True (to be removed, as soon
        as direct import of quadratic meshes is possible).

    Returns
    -------
    mesh: :gimliapi:`GIMLI::Mesh`
    """
    mesh = pg.Mesh(3)

    # Part 1/3: Nodes, essential
    with open(fName + '.node', 'r') as node_in:
        node_lines = pg.utils.filterLinesByCommentStr(node_in.readlines(),
                                                      comment)
    node_1 = node_lines[0].split()
    node_count = int(node_1[0])
    assert int(node_1[1]) == 3, 'Wrong dim: {}, should be 3.'.format(node_1[1])
    number_node_attr = int(node_1[2])
    node_markers = int(node_1[3])
    node_attributes = [] * number_node_attr

    for n in range(node_count):
        node_n = node_lines[n + 1].split()
        if node_markers:
            node_marker_n = int(node_n[-1])
        else:
            node_marker_n = n
        mesh.createNode([float(node_n[1]), float(node_n[2]), float(node_n[3])],
                        marker=node_marker_n)
        for m in range(number_node_attr):
            node_attributes[m].append(float(node_n[4 + m]))

    for k in range(number_node_attr):
        if verbose:
            print('Add node data to mesh.')
        mesh.addData('node_data_{}'.format(k + 1), node_attributes[k])

    # Part 2/3: Tetrahedrons, optional
    if os.path.exists(fName + '.ele'):
        if verbose:
            print('Found .ele file. Adding cells and cell marker.')
        with open(fName + '.ele', 'r') as cell_in:
            cell_lines = pg.utils.filterLinesByCommentStr(cell_in.readlines(),
                                                          comment)
        cell_1 = cell_lines[0].split()
        cell_count = int(cell_1[0])

        nodes_per_cell = int(cell_1[1])  # 4 or 10
        if nodes_per_cell == 10:
            quadratic = True
            raise Exception('Cannot import quadratic meshes directly yet.')

        cell_markers = int(cell_1[2])
        for c in range(cell_count):
            cell_n = cell_lines[c + 1].split()
            if cell_markers:
                cell_marker_n = int(cell_n[-1])
            else:
                cell_marker_n = defaultCellMarker
            mesh.createCell([int(ind) for ind in cell_n[1:5]],
                            marker=cell_marker_n)
            # in order to import quadratic meshes directly, i ned the sorting
            # of the node indices
#           mesh.createCell([int(ind) for ind in cell_n[1:nodes_per_cell + 1]],
#                           marker=cell_marker_n)

        # Part 3/3: Boundaries and Marker, optional
    if os.path.exists(fName + '.face') and loadFaces:
        if verbose:
            print('Found .face file. Adding boundaries and boundary marker.')
        with open(fName + '.face', 'r') as face_in:
            face_lines = pg.utils.filterLinesByCommentStr(face_in.readlines(),
                                                          comment)

        face_1 = face_lines[0].split()
        face_count = int(face_1[0])
        face_markers = int(face_1[1])

        for f in range(face_count):
            face_n = face_lines[f + 1].split()
            if face_markers:
                face_marker_n = int(face_n[-1])
            else:
                face_marker_n = 0
            mesh.createBoundary([int(ind) for ind in face_n[1:4]],
                                marker=face_marker_n)
            # quadratic
#            mesh.createBoundary(
#                [int(ind) for ind in face_n[1:len(face_n) - face_markers]],
#                marker=face_marker_n)

    if quadratic:
        mesh = mesh.createP2()

    if verbose:
        print(mesh)

    return mesh


def readHydrus2dMesh(fileName='MESHTRIA.TXT'):
    """
    Import mesh from Hydrus2D.

    Parameters
    ----------
    fName : str, optional
        Filename of Hydrus output file.

    See Also
    --------
    readHydrus3dMesh : Similar routine for three-dimensional meshes.

    References
    ----------
    .. http://www.pc-progress.com/en/Default.aspx?h3d-description
    """
    fid = open(fileName)
    line = fid.readline().split()
    if 'HYDRUS' in line:
        return readHydrusMeshV3(fileName)

    nNodes = int(line[1])
    nCells = int(line[3])
    mesh = pg.Mesh()
    for _ in range(nNodes):
        line = fid.readline().split()
        mesh.createNode(
            pg.RVector3(float(line[1]) / 100., float(line[2]) / 100., 0.))

    for _ in range(3):
        line = fid.readline()

    for _ in range(nCells):
        line = fid.readline().split()
        if len(line) == 4:
            mesh.createTriangle(
                mesh.node(int(line[1]) - 1), mesh.node(int(line[2]) - 1),
                mesh.node(int(line[3]) - 1), 1)
        elif len(line) == 5:
            mesh.createTetrahedron(
                mesh.node(int(line[1]) - 1), mesh.node(int(line[2]) - 1),
                mesh.node(int(line[3]) - 1), mesh.node(int(line[4]) - 1), 1)

    fid.close()
    mesh.createNeighborInfos()
    return mesh


def readHydrus3dMesh(fileName='MESHTRIA.TXT'):
    """Import mesh from Hydrus3D.

    Parameters
    ----------
    fName : str, optional
        Filename of Hydrus output file.

    See Also
    --------
    readHydrus2dMesh : Similar routine for two-dimensional meshes.

    References
    ----------
    .. http://www.pc-progress.com/en/Default.aspx?h3d-description
    """
    f = open(fileName, 'r')

    if 'HYDRUS' in f.readline():
        return readHydrusMeshV3(fileName)

    for i in range(5):
        line1 = f.readline()

    nNodes = int(line1.split()[0])
    nCells = int(line1.split()[1])
    # print(nNodes, nCells)

    line1 = f.readline()
    nodes = []
    dx = 0.01
    mesh = pg.Mesh()
    for _ in range(nNodes):
        pos = f.readline().split()
        p = pg.RVector3(
            float(pos[1]) * dx, float(pos[2]) * dx, float(pos[3]) * dx * (-1.))
        n = mesh.createNode(p)
        nodes.append(n)

    line1 = f.readline()
    line1 = f.readline()
    cells = []
    for _ in range(nCells):
        pos = f.readline().split()
        i, j, k, l = int(pos[1]), int(pos[2]), int(pos[3]), int(pos[4]),
        c = mesh.createTetrahedron(nodes[i - 1], nodes[j - 1], nodes[k - 1],
                                   nodes[l - 1])
        cells.append(c)

    f.close()
    mesh.createNeighborInfos()
    return mesh


def readHydrusMeshV3(fileName):
    """Import mesh from Hydrus3D. File Version 3.

    TODO:
    ----
    * 3D

    Parameters
    ----------
    fileName : str
        Filename of Hydrus output file.

    """
    with open(fileName) as fid:
        lines = fid.readlines()
        mesh = None
        nNodes = 0
        nBound = 0
        nCells = 0
        for i, line in enumerate(lines):
            if 'Mesh Dim' in line:
                dim = int(line.split('\r\n')[0].split(':')[1])
                if dim != 2:
                    pg.error('3D mesh not yet implemented.')

                mesh = pg.Mesh(dim=dim)

            if 'DIMENSIONS' in line:
                nNodes, nEle1, nEle2, nEle3 = lines[i+8].split(
                    '\r\n')[0].split()
                if mesh.dim() == 1:
                    nNodes = int(nNodes)
                    nCells = int(nEle1)
                elif mesh.dim() == 2:
                    nNodes = int(nNodes)
                    nBound = int(nEle1)
                    nCells = int(nEle2)
                else:
                    nNodes = int(nNodes)
                    nBound = int(nEle2)
                    nCells = int(nEle3)

                print(nNodes, nBound, nCells)

            if 'NODAL' in line:
                for l in lines[i+7:i+7+nNodes]:
                    vals = l.split('\r\n')[0].split()
                    if mesh.dim() == 1:
                        mesh.createNode([float(vals[1]), 0.0])
                    elif mesh.dim() == 2:
                        mesh.createNode([float(vals[1]), float(vals[2])])
                    elif mesh.dim() == 3:
                        mesh.createNode([float(vals[1]), float(vals[2]),
                                         float(vals[3])])
                i = i + 7 + nNodes

            if '1D-ELEMENTS' in line:
                for l in lines[i + 8:]:
                    vals = l.split('\r\n')[0].split()
                    if mesh.dim() == 2:
                        try:
                            mesh.createBoundary([int(vals[2])-1,
                                                 int(vals[3])-1],
                                                marker=int(vals[1]))
                        except BaseException as e:
                            if mesh.boundaryCount() != nBound:
                                print(e)
                                pg.error('Something is wrong in the file.'
                                         '{:d} 1D-Elements are announced but'
                                         'only {:d} found'.format(
                                            nBound, mesh.boundaryCount()))
                                i = i + 8 + mesh.boundaryCount()
                            break
                    else:
                        pg.error('Something wrong with mesh dimension')

            if '2D-ELEMENTS' in line:
                for l in lines[i + 10:]:
                    vals = l.split('\r\n')[0].split()
                    if mesh.dim() == 2:
                        try:
                            nn = [int(vals[2])-1, int(vals[3])-1,
                                  int(vals[4])-1]
                            # nn = [int(v)-1 for v in vals[2:5]]
                            if int(vals[5]) > 0:
                                nn.append(int(vals[5])-1)

                            mesh.createCell(nn, marker=int(vals[1]))
                        except BaseException as e:
                            if mesh.cellCount() != nCells:
                                print(e)
                                pg.error('Something is wrong in the file. '
                                    '{0} 2D-Elements announced but only '
                                    '{1} found'.format(nCells, mesh.cellCount()))

                                i = i + 10 + mesh.cellCount()
                                break
                    else:
                        pg.error('3D mesh not yet implemented.')

    return mesh


def readGambitNeutral(fileName, verbose=False):
    r"""Import Gambit Neutral meshes *.neu.

    See. https://www.sharcnet.ca/Software/Gambit/html/users_guide/ug01.htm

    Not fully implemented. If needed, we can improve this importer just send us
    an example file.

    Parameters
    ----------
    fName : string
        Filename of the file to read (\\*.n, \\*.e \\*.f)

    verbose : boolean, optional
        Be verbose during import.
    """
    with open(fileName, 'r') as fi:
        content = fi.readlines()
    fi.close()

    mesh = pg.Mesh(2)

    for i, line in enumerate(content):
        if 'ENDOFSECTION' in line:
            break

        try:
            nVertices = int(content[i - 1].split()[0])
            nElements = int(content[i - 1].split()[1])
        except:  # bare except!
            raise Exception("Cannot interpret GAMBIT Neutral header: " +
                            content[0:i])

    for i, line in enumerate(content):
        if 'NODAL COORDINATES' in line:
            break
    for j in range(nVertices):
        vx = float(content[i + j + 1].split()[1])
        vy = float(content[i + j + 1].split()[2])
        mesh.createNode((vx, vy))

    for i, line in enumerate(content):
        if 'ELEMENTS/CELLS' in line:
            break
    for j in range(nElements):
        nNodes = int(content[i + j + 1].split()[1])
        nodes = []
        for k in range(nNodes):
            nodes.append(int(content[i + 1 + j].split()[3 + k]) - 1)
        mesh.createCell(nodes)

    if verbose:
        print("Gambit neutral file imported: ", mesh)

    mesh.createNeighborInfos()
    return mesh


def convertHDF5Mesh(h5Mesh, group='mesh', indices='cell_indices',
                    pos='coordinates', cells='topology', marker='values',
                    marker_default=0, dimension=3, verbose=True,
                    useFenicsIndices=False):
    """Converts instance of a hdf5 mesh to a :gimliapi:`GIMLI::Mesh`.

    For full documentation please see :py:mod:`pygimli:meshtools:readHDF5Mesh`.
    """
    # open mesh containing group inside hdf file
    inmesh = h5Mesh.get(group)

    # get node positions
    mesh_pos = inmesh[pos][:]

    # get cells
    if useFenicsIndices:
        # case 1/2: using Fenics specific indices
        # TODO why extra indices? cause this is really slow
        # figure out the nature and purpose of the extra Fenics indices
        mesh_indices = inmesh[indices]
        try:
            mesh_cells = inmesh[cells][mesh_indices]
        except IndexError as IE:
            print('Fenics Indices aren\'t just in arbitrary order in range\
(0, cellCount) as expected. Need Fix.')
            raise IE
    else:
        # case 2/2: indices implicit: [0, ... cellCount)
        mesh_cells = inmesh[cells][:]

    # get marker
    try:
        # case 1/3: marker found in hdf file
        mesh_marker = inmesh[marker][:]
    except KeyError:
        if isinstance(marker_default, int):
            # case 2/3: marker not found, given as int
            mesh_marker = np.ones(len(mesh_cells), dtype=int) * marker_default
        else:
            # case 3/3: marker not found, given as array_like
            mesh_marker = marker_default

    # start building pygimli mesh
    mesh = pg.Mesh(dimension)
    for node in mesh_pos:
        mesh.createNode(node)

    for i, cell in enumerate(mesh_cells):
        mesh.createCell(pg.core.IndexArray(cell), marker=int(mesh_marker[i]))

    mesh.createNeighborInfos()

    if verbose:
        print('converted mesh:', mesh)
    return mesh


def readHDF5Mesh(fileName, group='mesh', indices='cell_indices',
                 pos='coordinates', cells='topology', marker='values',
                 marker_default=0, dimension=3, verbose=True,
                 useFenicsIndices=False):
    """Function for loading a mesh from HDF5 file format.

    Returns an instance of :gimliapi:`GIMLI::Mesh` class.
    Default values for keywords are suited for :term:`FEniCS` syntax
    .h5 meshes.

    Requirements: h5py module

    TODO:
        * Fenics hdf5 meshes do not have boundary markers.

    Parameters
    ----------

    fileName: string
        Name of the mesh to be transformed into :term:`pyGIMLi` format.

    group: string ['domains']
        hdf group that contains the mesh information (see other keyword
        arguments). Default is 'domains' for :term:`FEniCS` compatibility.

    indices: string ['cell_indices']
        Key for the part of the hdf file containing the indices of the cells.

    pos: string ['coordinates']
        Key for the part of the hdf file containing the nodepositions.

    cells: string ['topology']
        Key for the part of the hdf file containing the array which defies the
        cells. Usually of shape (cellCount, 3) for 2D meshes or (cellCount, 4)
        for 3D tetrahedra meshes. For each cell the indices of the
        corresponding node indices is given.

    marker: string ['values']
        If marker is part of the hdf data container, the corresponding array
        is used as identifier for the cell markers. If not found, the cell
        markers will be set to marker_default.

    marker_default: int or array [0]
        Default marker if no markers are found in the hdf file. If array, size
        has to match the cellCount of the mesh.

    dimension: int [3]
        Dimension of the input/output mesh, no own check for dimensions yet.
        Fixed on 3 for now.

    Returns
    -------

    mesh:
        :gimliapi:`GIMLI::Mesh`

    """
    h5py = pg.optImport('h5py',
                        requiredFor='import mesh in .h5 data format')
    h5 = h5py.File(fileName, 'r')
    if verbose:
        print('loaded hdf5 mesh:', h5)

    mesh = convertHDF5Mesh(h5, group=group, indices=indices, pos=pos,
                           cells=cells, marker=marker,
                           marker_default=marker_default, dimension=dimension,
                           verbose=verbose, useFenicsIndices=useFenicsIndices)

    h5.close()
    return mesh


def readFenicsHDF5Mesh(fileName, group='mesh', verbose=True):
    """ Reads :term:`FEniCS` mesh from file format .h5 and returns a
    :gimliapi:`GIMLI::Mesh`.
    """
    mesh = readHDF5Mesh(fileName, group=group, indices='cell_indices',
                        pos='coordinates', cells='topology', marker='values',
                        marker_default=0, dimension=3, verbose=verbose,
                        useFenicsIndices=False)
    return mesh


def exportHDF5Mesh(mesh, exportname, group='mesh', indices='cell_indices',
                   pos='coordinates', cells='topology', marker='values'):
    """Writes given :gimliapi:`GIMLI::Mesh` in a hdf5 format file.

    3D tetrahedral meshes only! Boundary markers are ignored.

    Keywords are explained in :py:mod:`pygimli.meshtools.readHDFS`
    """
    h5py = pg.optImport('h5py',
                        requiredFor='export mesh in .h5 data format')

    if not isinstance(mesh, pg.Mesh):
        mesh = pg.Mesh(mesh)

    # prepare output for writing in hdf data container
    pg_pos = mesh.positions()
    mesh_pos = np.array((np.array(pg.x(pg_pos)), np.array(pg.y(pg_pos)),
                         np.array(pg.z(pg_pos)))).T

    mesh_cells = np.zeros((mesh.cellCount(), 4))  # hard coded for tetrahedrons
    for i, cell in enumerate(mesh.cells()):
        mesh_cells[i] = cell.ids()

    mesh_indices = np.arange(0, mesh.cellCount() + 1, 1, dtype=np.int64)
    mesh_markers = np.array(mesh.cellMarkers())

    with h5py.File(exportname, 'w') as out:
        for grp in np.atleast_1d(group):  # can use more than one group
            # writing indices
            idx_name = '{}/{}'.format(grp, indices)
            out.create_dataset(idx_name, data=mesh_indices, dtype=int)
            # writing node positions
            pos_name = '{}/{}'.format(grp, pos)
            out.create_dataset(pos_name, data=mesh_pos, dtype=float)
            # writing cells via indices
            cells_name = '{}/{}'.format(grp, cells)
            out.create_dataset(cells_name, data=mesh_cells, dtype=int)
            # writing marker
            marker_name = '{}/{}'.format(grp, marker)
            out.create_dataset(marker_name, data=mesh_markers, dtype=int)
            out[grp][cells].attrs['celltype'] = np.string_('tetrahedron')
            out[grp][cells].attrs.create('partition', [0])
    return True


def exportFenicsHDF5Mesh(mesh, exportname):
    """Exports Gimli mesh in HDF5 format suitable for Fenics.

    Equivalent to calling the function
    :py:mod:`pygimli.meshtools.exportHDF5Mesh(mesh, exportname, group=['mesh',
    'domains'], indices='cell_indices', pos='coordinates', cells='topology',
    marker='values')`.

    Parameters
    ----------

    mesh: :gimliapi:GIMLI::Mesh`
        Mesh to be saved.

    exportname: string
        Name under which the mesh is saved.
    """
    return exportHDF5Mesh(mesh, exportname, group=['mesh', 'domains'],
                          indices='cell_indices', pos='coordinates',
                          cells='topology', marker='values')


def readEIDORSMesh(fileName, matlabVarname, verbose=False):
    """Reads finite element model in EIDORS format and returns pygimli mesh.

    Parameters
    ----------
    fileName : str
        name of the .mat file containing the EIDORS model
    matlabVarname : str
        variable name of .mat file in MATLAB workspace
    """
    scipy = pg.optImport("scipy", requiredFor="read EIDORS mesh.")

    def toDict(matobj):
        d = {}
        for name in matobj._fieldnames:
            elem = matobj.__dict__[name]
            if isinstance(elem, scipy.io.matlab.mio5_params.mat_struct):
                d[name] = toDict(elem)
            else:
                d[name] = elem
        return d

    def check_keys(d):
        for key in d:
            if isinstance(d[key], scipy.io.matlab.mio5_params.mat_struct):
                d[key] = toDict(d[key])
        return d

    def loadmat(fileName):
        data = scipy.io.loadmat(fileName, struct_as_record=False, squeeze_me=True)
        return check_keys(data)

    def get_nested(data, *args):
        if args and data:
            element = args[0]
            if element:
                value = data.get(element)
                return value if len(args) == 1 else get_nested(value,
                                                               *args[1:])

    matlab_eidors = loadmat(fileName)
    python_eidors = get_nested(matlab_eidors, matlabVarname)
    # if input eidors data is forward model instead of an image
    if 'nodes' in python_eidors.keys():
        nodes = get_nested(python_eidors, "nodes")
        elems = get_nested(python_eidors, "elems")
        boundary_numbers = get_nested(python_eidors, "boundary_numbers")

    # it is an image with elem_data and fwd_model
    else:
        nodes = get_nested(python_eidors, "fwd_model", "nodes")
        elems = get_nested(python_eidors, "fwd_model", "elems")
        boundary_numbers = get_nested(python_eidors, "fwd_model",
                                      "boundary_numbers")

    dim_nodes = np.size(nodes, 1)
    dim_elems = np.size(elems, 1)

    if verbose:
        print('Reading %s... ' % fileName)
        print('found %s  %s-dimensional nodes... ' % (np.size(nodes, 0),
                                                      dim_nodes))

    if dim_elems == 3 and verbose:
        print('found %s triangles... ' % np.size(elems, 0))
    elif dim_elems == 4 and verbose:
        print('found %s tetrahedrons... ' % np.size(elems, 0))

    if 'elem_data' in python_eidors.keys():
        elem_data = get_nested(python_eidors, "elem_data")
        if verbose:
            print('found %s element data... ' % len(elem_data))
    else:
        if verbose:
            print("found no element data... ")

    region_markers = np.unique(boundary_numbers)
    no_of_regions = len(region_markers)

    if boundary_numbers is not None:
        if verbose:
            print('Found %s unique Regions with Markers' % no_of_regions)
            print('%s' % region_markers)

    if boundary_numbers is None:
        if verbose:
            print('Found no Unique Region with Region Markers')

    # create nodes from eidors model
    mesh = pg.Mesh()

    # if two dimensional eidors model
    if (dim_nodes == 2 and dim_elems == 3):
        if verbose:
            print('converting to %d-D pygimli mesh... ' % dim_nodes)
        for i in range(len(nodes)):
            mesh.createNode(pg.RVector3(nodes[i, 0], nodes[i, 1], 0.))
        for i in range(len(elems)):
            mesh.createTriangle(
                mesh.node(int(elems[i, 0]) - 1),
                mesh.node(int(elems[i, 1]) - 1),
                mesh.node(int(elems[i, 2]) - 1), 1)

    # for non-planar 2D models
    if (dim_nodes == 3 and dim_elems == 3):
        if verbose:
            print('converting to pygimli mesh...')
            print('found 3d nodes with 2d elements')
        for i in range(len(nodes)):
            mesh.createNode(pg.RVector3(nodes[i, 0], nodes[i, 1], nodes[i, 2]))
        for i in range(len(elems)):
            mesh.createTriangle(
                mesh.node(int(elems[i, 0]) - 1),
                mesh.node(int(elems[i, 1]) - 1),
                mesh.node(int(elems[i, 2]) - 1), 1)

    # if three dimensional eidors model
    if (dim_nodes == 3 and dim_elems == 4):
        if verbose:
            print('converting to %d-D pygimli mesh... ' % dim_nodes)
        for i in range(len(nodes)):
            mesh.createNode(pg.RVector3(nodes[i, 0], nodes[i, 1], nodes[i, 2]))
        for i in range(len(elems)):
            mesh.createTetrahedron(
                mesh.node(int(elems[i, 0]) - 1),
                mesh.node(int(elems[i, 1]) - 1),
                mesh.node(int(elems[i, 2]) - 1),
                mesh.node(int(elems[i, 3]) - 1), 1)

    if boundary_numbers is not None:
        mesh.setCellMarkers(boundary_numbers.astype(int))

    return mesh


def readSTL(fileName, ascii=True):
    """Read :term:`STL` surface mesh and returns a :gimliapi:`GIMLI::Mesh`.

    Read :term:`STL` surface mesh and returns a :gimliapi:`GIMLI::Mesh`
    of triangle boundary faces. Multiple solids are supported with increasing
    boundary marker.

    TODO: ASCII=False, read binary STL

    Parameters
    ----------

    fileName : str
        name of the .stl file containing the STL surface mesh

    ascii : bool [True]
        STL ASCII format
    """
    mesh = pg.Mesh(dim=3)
    mesh.importSTL(fileName)
    return mesh

    # readPos = lambda s: pg.pos(float(s[0], float(s[1]), float(s[2])))

    # with open(fileName, 'r') as fi:
    #     content = fi.readlines()
    # fi.close()

    # marker = -1

    # for i, line in enumerate(content):
    #     if 'solid' in line:
    #         marker += 1
    #     elif 'facet' in line:
    #         norm = readPos(content[i].split()[2:5])
    #         v1 = readPos(content[i+2].split()[1:4])
    #         v2 = readPos(content[i+2].split()[1:4])
    #         v3 = readPos(content[i+2].split()[1:4])

    #         mesh.createBoundary([v1, v2, v3], marker=marker)
    #         i += 7

    # return mesh


def exportSTL(mesh, fileName, ascii=True):
    """Write :term:`STL` surface mesh and returns a :gimliapi:`GIMLI::Mesh`.

    Export a three dimensional boundary :gimliapi:`GIMLI::Mesh` into a
    :term:`STL` surface mesh. Boundaries with different marker
    will be separated into different STL solids.

    TODO:
        * ASCII=False, write binary STL
        * QuadrangleFace Boundaries
        * p2 Boundaries

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        Mesh to be exported. Only Boundaries of type TriangleFace will be
        exported.

    fileName : str
        name of the .stl file containing the STL surface mesh

    ascii : bool [True]
        STL ASCII format
    """
    marker = pg.unique(pg.sort(mesh.boundaryMarkers()))

    if '.stl' not in fileName:
        fileName = fileName + '.stl'

    fi = open(fileName, 'w')
    for m in marker:
        me = mesh.createSubMesh(mesh.boundaries(mesh.boundaryMarkers() == m))

        fi.write('solid ' + str(m) + '\n')

        for b in me.boundaries():
            n = b.norm()
            fi.write('facet normal %f %f %f\n' % (n[0], n[1], n[2]))
            fi.write('\touter loop\n')
            fi.write('\t\tvertex %f %f %f\n' % (b.node(0).pos()[0],
                                                b.node(0).pos()[1],
                                                b.node(0).pos()[2]))
            fi.write('\t\tvertex %f %f %f\n' % (b.node(1).pos()[0],
                                                b.node(1).pos()[1],
                                                b.node(1).pos()[2]))
            fi.write('\t\tvertex %f %f %f\n' % (b.node(2).pos()[0],
                                                b.node(2).pos()[1],
                                                b.node(2).pos()[2]))
            fi.write('\tendloop\n')
            fi.write('endfacet\n')

        fi.write('endsolid\n')

    fi.close()


def transform2DMeshTo3D(mesh, x, y, z=None):
    """2D mesh into 3D coordinates.

    Transform a 2D mesh into 3D coordinates using a point list (e.g. from GPS)

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`

    x,y: [float]
        array of x/y positions along 2d profile

    z: [float]
        optional height to add (topographical correction on top of flat earth)

    See Also
    --------

    References
    ----------
    """
    # get mesh node positions
    mt, mz = pg.x(mesh.positions()), pg.y(mesh.positions())  # mesh tape and z
    # compute length of reference points along tape
    pt = np.hstack((0., np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))))
    #  interpolate node positions from tape to x/y using tape positions
    mx = np.interp(mt, pt, x)
    my = np.interp(mt, pt, y)
    # compute z offset by interpolating z
    if z is None:
        oz = np.zeros(len(mt))
    else:
        oz = np.interp(mt, pt, z)

    # set the positions in the mesh
    for i, node in enumerate(mesh.nodes()):
        node.setPos(pg.RVector3(mx[i], my[i], mz[i] + oz[i]))


def rot2DGridToWorld(mesh, start, end):
    """Rotate a 2D mesh into 3D world coordinates.

    todo:: Complete Documentation. ...rotate a given 2D grid in...
    """
    mesh.rotate(pg.degToRad(pg.RVector3(-90.0, 0.0, 0.0)))

    src = pg.RVector3(0.0, 0.0, 0.0).norm(pg.RVector3(0.0, 0.0, -10.0),
                                          pg.RVector3(10.0, 0.0, -10.0))
    dest = start.norm(start - pg.RVector3(0.0, 0.0, 10.0), end)

    rot = pg.getRotation(src, dest)
    mesh.transform(rot)
    mesh.translate(start)


def merge2Meshes(m1, m2):
    """Merge two meshes into one new mesh and return the combined mesh.

    Merge two meshes into a new mesh and return the combined mesh.
    Note that there is a duplicate check for all nodes which should reuse
    existing node but NO cells or boundaries.

    Parameters
    ----------
    m1: :gimliapi:`GIMLI::Mesh`
        First mesh.
    m2: :gimliapi:`GIMLI::Mesh`
        Second mesh.

    Returns
    -------
    mesh: :gimliapi:`GIMLI::Mesh`
        Resulting mesh.
    """
#    for c in m1.cells():
#    if c.size() < 1e-4:
#    print(c)
#    exit()

    mesh = pg.Mesh(m1)
    mesh.translate(-m1.node(0).pos())
    m3 = pg.Mesh(m2)
    m3.translate(-m1.node(0).pos())

#    for n in m3.nodes():
#    i = mesh.findNearestNode(n.pos())
#    if mesh.node(i).pos().dist(n.pos()) < 0.5:
#    print("DUP", 1)
#    exit()

    for c in m3.cells():
        mesh.copyCell(c)

    for b in m3.boundaries():
        mesh.copyBoundary(b, tol=1e-6, check=False)

#    if b.id() > 1362:
#    exit()
#    print(mesh.boundary(2905),
#    mesh.boundary(2905).node(0).id(), mesh.boundary(2905).node(0).pos(),
#    mesh.boundary(2905).node(1).id(), mesh.boundary(2905).node(1).pos()
#    )
#    print(mesh.boundary(2906),
#    mesh.boundary(2906).node(0).id(), mesh.boundary(2906).node(0).pos(),
#    mesh.boundary(2906).node(1).id(), mesh.boundary(2906).node(1).pos()
#    )

    for key in list(mesh.dataMap().keys()):
        d = mesh.dataMap()[key]
        d.resize(mesh.cellCount())
        d.setVal(m1.dataMap()[key], 0, m1.cellCount())
        d.setVal(m2.dataMap()[key], m1.cellCount(),
                 m1.cellCount() + m2.cellCount())
        mesh.addData(key, d)

    mesh.translate(m1.node(0).pos())
    return mesh


def mergeMeshes(meshList, verbose=False):
    """Merge several meshes into one new mesh and return the new mesh.

    Merge several meshes into one new mesh and return the new mesh.

    Parameters
    ----------
    meshList : [:gimliapi:`GIMLI::Mesh`, ...] | [str, ...]
        List of at least two meshes (or filenames to meshes) to be merged.

    verbose : bool
        Give some output

    See Also
    --------
    merge2Meshes
    """
    if not isinstance(meshList, list):
        raise Exception("Argument meshList is no list")

    if len(meshList) < 2:
        raise Exception(
            "To few meshes in meshList, at least 2 meshes are needed.")

    if isinstance(meshList[0], str):
        mL = []
        if verbose:
            print("Reading meshes ... ")

        for mFileName in meshList:
            m = pg.Mesh(2)
            m.load(mFileName)
            if verbose:
                print("loaded", mFileName, m)
            mL.append(m)

        meshList = mL
        #meshList = [pg.Mesh(2); m.load(m) ]

    if verbose:
        print("Merging meshes ... ")

    mesh = meshList[0]
    if verbose:
        print(mesh)

    for m in range(1, len(meshList)):
        mesh = merge2Meshes(mesh, meshList[m])
        if verbose:
            print(mesh)

    return mesh


def createParaMesh(*args, **kwargs):
    """Create parameter mesh from list of sensor positions.

    Create parameter mesh from list of sensor positions.

    Parameters
    ----------
    sensors : list of RVector3 objects
        Sensor positions. Must be sorted and unique in positive x direction.
        Depth need to be y-coordinate.
    paraDX : float
        Relative distance for refinement nodes between two electrodes (1=none),
        e.g., 0.5 means 1 additional node between two neighboring electrodes
        e.g., 0.33 means 2 additional equidistant nodes between two electrodes
    paraDepth : float, optional
        Maximum depth for parametric domain, 0 (default) means 0.4 * maximum
        sensor range.
    paraBoundary : float, optional
        Margin for parameter domain in absolute sensor distances. 2 (default).
    paraMaxCellSize: double [0], optional
        Maximum size for parametric size in m*m (0-no constraint)
    boundaryMaxCellSize: double [0], optional
        Maximum cells size in the boundary region in m*m (0-no constraint)
    boundary : float, optional
        Boundary width to be appended for domain prolongation in absolute
        para domain width.
        Values <=0 force the boundary to be 4 times para domain width.

    Returns
    -------
    poly: :gimliapi:`GIMLI::Mesh`
    """
    plc = pg.meshtools.createParaMeshPLC(*args, **kwargs)
    kwargs.pop('paraMaxCellSize', 0)
    kwargs.pop('boundaryMaxCellSize', 0)
    mesh = createMesh(plc, **kwargs, smooth=[2, 10])
    return mesh


def createParaMesh2dGrid(*args, **kwargs):
    """API change here .. use createParaMesh2DGrid instead."""
    pg.deprecated("pls use createParaMesh2DGrid")
    return createParaMesh2DGrid(*args, **kwargs)


def createParaMesh2DGrid(sensors, paraDX=1, paraDZ=1, paraDepth=0, nLayers=11,
                         boundary=-1, paraBoundary=2, **kwargs):
    """Create a grid-style mesh for an inversion parameter mesh.

    Create a grid-style mesh for an inversion parameter mesh.
    Return parameter grid for a given list of sensor positions.
    Uses and forwards arguments to
    :py:mod:`pygimli.meshtools.appendTriangleBoundary`.

    Parameters
    ----------
    sensors : list of RVector3 objects or data container with sensorPositions
        Sensor positions. Must be sorted in positive x direction
    paraDX : float, optional
        Horizontal distance between sensors, relative regarding sensor
        distance. Value must be greater than 0 otherwise 1 is assumed.
    paraDZ : float, optional
        Vertical distance to the first depth layer, relative regarding sensor
        distance. Value must be greater than 0 otherwise 1 is assumed.
    paraDepth : float, optional
        Maximum depth for parametric domain, 0 (default) means 0.4 * maximum
        sensor range.
    nLayers : int, optional [11]
        Number of depth layers.
    boundary : int, optional [-1]
        Boundary width to be appended for domain prolongation in absolute
        para domain width.
        Values lower than 0 force the boundary to be 4 times para domain width.
    paraBoundary : int, optional [2]
        Offset to the parameter domain boundary in absolute sensor spacing.

    Returns
    -------
    mesh: :gimliapi:`GIMLI::Mesh`

    Examples
    --------
    >>> import pygimli as pg
    >>> import matplotlib.pyplot as plt
    >>>
    >>> from pygimli.meshtools import createParaMesh2DGrid
    >>> mesh = createParaMesh2DGrid(sensors=pg.Vector(range(10)),
    ...                             boundary=1, paraDX=1,
    ...                             paraDZ=1, paraDepth=5)
    >>> ax, _ = pg.show(mesh, markers=True, showMesh=True)
    """
    mesh = pg.Mesh(2)

    # maybe separate x y z and sort
    if isinstance(sensors, np.ndarray) or isinstance(sensors, pg.Vector):
        sensors = [pg.RVector3(s, 0) for s in sensors]

    if isinstance(sensors, pg.DataContainer):
        sensors = sensors.sensorPositions()

    sensorX = pg.x(sensors)

    eSpacing = abs(sensorX[1] - sensorX[0])

    xMin = min(sensorX) - paraBoundary * eSpacing
    xMax = max(sensorX) + paraBoundary * eSpacing

    if paraDX == 0:
        paraDX = 1.
    if paraDZ == 0:
        paraDZ = 1.

    dx = paraDX
    dz = paraDZ
    if eSpacing > 0:
        dx = eSpacing * paraDX
        # dz = eSpacing * paraDZ  # not really making sense

    if paraDepth == 0:
        paraDepth = 0.4 * (xMax - xMin)

    # print(xMin, xMax, dx)
    x = pg.utils.grange(xMin, xMax, dx=dx)

    y = -pg.core.increasingRange(dz, paraDepth, nLayers)

    mesh.createGrid(x, y)
    mesh.setCellMarkers([2] * mesh.cellCount())

    paraXLimits = [xMin, xMax]
    #    paraYLimits = [min(y), max(y)]  # not used

    if boundary < 0:
        boundary = abs((paraXLimits[1] - paraXLimits[0]) * 4.0)

    mesh = pg.meshtools.appendTriangleBoundary(
        mesh, xbound=boundary, ybound=boundary, marker=1, **kwargs)

    return mesh


if __name__ == "__main__":
    pass
