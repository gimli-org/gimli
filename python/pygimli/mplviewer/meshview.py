# -*- coding: utf-8 -*-
"""Draw mesh/model/fields with matplotlib."""

import textwrap

import numpy as np

import matplotlib as mpl
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.colors import LogNorm

import pygimli as pg
from pygimli.misc import streamline
from .utils import updateAxes as updateAxes_

from .colorbar import cmapFromName, autolevel


class CellBrowser(object):
    """Interactive cell browser on current or specified ax for a given mesh.

    Cell information can be displayed by mouse picking. Arrow keys up and down
    can be used to scroll through the cells, while ESC closes the cell
    information window.

    Parameters
    ----------
    mesh : 2D pygimli.Mesh instance
        The plotted mesh to browse through.
    ax : mpl axis instance, optional
        Axis instance where the mesh is plotted (default is current axis).

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import pygimli as pg
    >>> from pygimli.mplviewer import drawModel
    >>> from pygimli.mplviewer import CellBrowser
    >>>
    >>> mesh = pg.createGrid(range(5), range(5))
    >>> fig, ax = plt.subplots()
    >>> plc = drawModel(ax, mesh, mesh.cellMarkers())
    >>> browser = CellBrowser(mesh)
    >>> browser.connect()
    ('Interactive cell browser activated on Fig.', ...)
    """

    def __init__(self, mesh, data=None, ax=None):
        """Construct CellBrowser on a specific `mesh`."""
        if ax:
            self.ax = ax
        else:
            self.ax = mpl.pyplot.gca()

        self.fig = self.ax.figure
        self.mesh = mesh
        self.lw = np.ones(mesh.cellCount())
        self.data = data
        self.cellID = None
        self.edgeColors = None
        self.event = None
        self.artist = None
        self.pid = None
        self.kid = None

        if data is not None:
            if mesh.cellCount() != len(data):
                print('data length mismatch mesh.cellCount(): ' +
                      str(len(data)) + "!=" + str(mesh.cellCount()) +
                      ". Mapping data to cellMarkers().")
                self.data = data[mesh.cellMarkers()]

        bbox = dict(boxstyle='round, pad=0.5', fc='w', alpha=0.5)
        arrowprops = dict(arrowstyle='->', connectionstyle='arc3,rad=0.5')
        kwargs = dict(fontproperties='monospace', visible=False,
                      fontsize=mpl.rcParams['font.size'] - 2,
                      weight='bold', xytext=(50, 20), arrowprops=arrowprops,
                      textcoords='offset points', bbox=bbox, va='center')
        self.text = self.ax.annotate(None, xy=(0, 0), **kwargs)

    def connect(self):
        """Connect to matplotlib figure canvas."""
        self.pid = self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.kid = self.fig.canvas.mpl_connect('key_press_event', self.onpress)
        print(("Interactive cell browser activated on Fig.", self.fig.number))

    def disconnect(self):
        """Disconnect from matplotlib figure canvas."""
        self.fig.canvas.mpl_connect(self.pid)
        self.fig.canvas.mpl_connect(self.kid)
        print(("Cell browser disconnected from Figure", self.fig.number))

    def hide(self):
        """Hide info window."""
        self.text.set_visible(False)
        self.artist.set_edgecolors(self.edgeColors)
        self.fig.canvas.draw()

    def highlight(self):
        """Highlight selected cell."""
        if self.edgeColors:
            ec = self.edgeColors.copy()
            ec[self.cellID] = np.ones(ec.shape[1])
            self.artist.set_edgecolors(ec)
            lw = self.lw.copy()
            lw[self.cellID] = 3
            self.artist.set_linewidths(lw)

    def onpick(self, event):
        """Call `self.update()` on mouse pick event."""
        self.event = event
        self.artist = event.artist

        if self.data is None:
            self.data = self.artist.get_array()
            self.edgeColors = self.artist.get_edgecolors()

        if 'mouseevent' in event.__dict__.keys():
            if event.mouseevent.xdata is not None and \
                event.mouseevent.ydata is not None:
                    pos = pg.RVector3(event.mouseevent.xdata, event.mouseevent.ydata)
                    self.cellID = self.mesh.findCell(pos, True).id()
        else:  # variant before (seemed inaccurate)
            self.cellID = event.ind[0]
        self.update()

    def onpress(self, event):
        """Call `self.update()` if up, down, or escape keys are pressed."""
        # print(event, event.key)
        if self.data is None:
            return
        if event.key not in ('up', 'down', 'escape'):
            return
        if event.key is 'up':
            self.cellID += 1
        elif event.key is'down':
            self.cellID -= 1
        else:
            self.hide()
            return
        self.cellID = int(np.clip(self.cellID, 0, self.mesh.cellCount() - 1))
        self.update()

    def update(self):
        """Update the information window."""
        center = self.mesh.cell(self.cellID).center()
        x, y = center.x(), center.y()
        marker = self.mesh.cells()[self.cellID].marker()
        data = self.data[self.cellID]
        header = "Cell %d:\n" % self.cellID
        header += "-" * (len(header) - 1)
        info = "\nx: {:.2f}\n y: {:.2f}\n data: {:.2e}\n marker: {:d}".format(
            x, y, data, marker)
        text = header + textwrap.dedent(info)
        self.text.set_text(text)
        self.text.xy = x, y
        self.text.set_visible(True)
        self.highlight()
        self.fig.canvas.draw()


def drawMesh(ax, mesh, **kwargs):
    """Draw a 2d mesh into a given ax.

    Set the limits of the ax tor the mesh extent.

    Parameters
    ----------

    mesh : :gimliapi:`GIMLI::Mesh`
        The plotted mesh to browse through.
    ax : mpl axis instance, optional
        Axis instance where the mesh is plotted (default is current axis).
    fitView: bool [True]
        Adjust ax limits to mesh bounding box.

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pygimli as pg
    >>> from pygimli.mplviewer import drawMesh
    >>> n = np.linspace(1, 2, 10)
    >>> mesh = pg.createGrid(x=n, y=n)
    >>> fig, ax = plt.subplots()
    >>> drawMesh(ax, mesh)
    >>> plt.show()
    """
    pg.mplviewer.drawMeshBoundaries(ax, mesh, **kwargs)

    if kwargs.pop('fitView', True):
        ax.set_xlim(mesh.xmin(), mesh.xmax())
        ax.set_ylim(mesh.ymin(), mesh.ymax())

    ax.set_aspect('equal')
    updateAxes_(ax)


def drawModel(ax, mesh, data=None,
              cMin=None, cMax=None, logScale=True, cmap=None,
              xlabel=None, ylabel=None, verbose=False, grid=False,
              tri=False,
              **kwargs):
    """Draw a 2d mesh and color the cell by the data.

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        The plotted mesh to browse through.
    ax : mpl axis instance, optional
        Axis instance where the mesh is plotted (default is current axis).
    data : array, optional
        Data to draw. Should either equal numbers of cells or nodes of the
        corresponding `mesh`.
    grid : boolean, optional
        Draw cell boundaries.
    tri : boolean, optional
        use MPL tripcolor (experimental)

    **kwargs:
        * all remaining
            Will be forwarded to the draw functions and matplotlib methods,
            respectively.
    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pygimli as pg
    >>> from pygimli.mplviewer import drawModel
    >>> n = np.linspace(0, -2, 11)
    >>> mesh = pg.createGrid(x=n, y=n)
    >>> mx = pg.x(mesh.cellCenter())
    >>> my = pg.y(mesh.cellCenter())
    >>> data = np.cos(1.5 * mx) * np.sin(1.5 * my)
    >>> fig, ax = plt.subplots()
    >>> drawModel(ax, mesh, data)
    <matplotlib.collections.PolyCollection object at ...>
    """
    if mesh.nodeCount() == 0:
        raise "drawModel: The mesh is empty."

    if tri:
        gci = drawMPLTri(ax, mesh, data, cMin=cMin, cMax=cMax,
                         cmap=cmap, **kwargs)
    else:
        gci = pg.mplviewer.createMeshPatches(ax, mesh, verbose=verbose,
                                             **kwargs)

        cMap = kwargs.pop('cMap', None)
        if cMap is not None:
            cmap = cMap

        if cmap is not None:
            if isinstance(cmap, str):
                gci.set_cmap(cmapFromName(cmap))
            else:
                gci.set_cmap(cmap)

        ax.set_aspect('equal')

        gci.set_antialiased(True)

        if grid:
            gci.set_linewidth(0.3)
            gci.set_edgecolor("0.3")
        else:
            gci.set_linewidth(None)

        if data is None:
            data = pg.RVector(mesh.cellCount())

        if len(data) != mesh.cellCount():

            print(data, mesh)
            print("INFO: drawModel have wrong data length .. " +
                  " indexing data from cellMarkers()")
            viewdata = data(mesh.cellMarkers())
        else:
            viewdata = data

        if min(data) <= 0:
            logScale = False

        pg.mplviewer.setMappableData(gci, viewdata, cMin=cMin, cMax=cMax,
                                     logScale=logScale)

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)

    showMesh = kwargs.pop('showMesh', False)
    if showMesh:
        drawMesh(ax, mesh)

    updateAxes_(ax)
    return gci


def drawSelectedMeshBoundaries(ax, boundaries,
                               color=None, linewidth=1.0):
    """Draw mesh boundaries into a given ax."""
    drawAA = True
    lines = []

    if hasattr(boundaries, '__len__'):
        if len(boundaries) == 0:
            return

    for bound in boundaries:
        lines.append(list(zip([bound.node(0).x(), bound.node(1).x()],
                              [bound.node(0).y(), bound.node(1).y()])))

    lineCollection = mpl.collections.LineCollection(lines, antialiaseds=drawAA)

    if color is None:
        viewdata = [b.marker() for b in boundaries]
        pg.mplviewer.setMappableData(lineCollection,
                                     viewdata,
                                     logScale=False)
    else:
        lineCollection.set_color(color)

    lineCollection.set_linewidth(linewidth)
    ax.add_collection(lineCollection)

    updateAxes_(ax)

    return lineCollection


def drawSelectedMeshBoundariesShadow(ax, boundaries, first='x', second='y',
                                     color=(0.5, 0.5, 0.5, 1.0)):
    """TODO Documentme."""
    polys = []

    for cell in boundaries:
        polys.append(list(zip([getattr(cell.node(0), first)(),
                               getattr(cell.node(1), first)(),
                               getattr(cell.node(2), first)()],
                              [getattr(cell.node(0), second)(),
                               getattr(cell.node(1), second)(),
                               getattr(cell.node(2), second)()])))

    collection = mpl.collections.PolyCollection(polys, antialiaseds=True)

    collection.set_color(color)
    collection.set_edgecolor(color)
    collection.set_linewidth(0.2)
    ax.add_collection(collection)

    updateAxes_(ax)
    return collection


def drawMeshBoundaries(ax, mesh, hideMesh=False, useColorMap=False, **kwargs):
    """Draw mesh on ax with boundary conditions colorized.

    Parameters
    ----------

    hideMesh: bool [False]
        Show only the boundary of the mesh and omit inner edges that
        separate the cells.

    useColorMap: bool[False]
        Apply the default colormap to boundaries with marker values > 0

    **kwargs:
        * fitView : bool [True]
        * linewidth : float [0.3]
            linewidth for edges with marker == 0 if hideMesh is False.

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pygimli as pg
    >>> from pygimli.mplviewer import drawMeshBoundaries
    >>> n = np.linspace(0, -2, 11)
    >>> mesh = pg.createGrid(x=n, y=n)
    >>> for bound in mesh.boundaries():
    ...     if not bound.rightCell():
    ...         bound.setMarker(pg.MARKER_BOUND_MIXED)
    ...     if bound.center().y() == 0:
    ...         bound.setMarker(pg.MARKER_BOUND_HOMOGEN_NEUMANN)
    >>> fig, ax = plt.subplots()
    >>> drawMeshBoundaries(ax, mesh)
    """
    if not mesh:
        raise Exception("drawMeshBoundaries(ax, mesh): invalid mesh")

    if not mesh.dimension() == 2:
        raise Exception("No 2d mesh: dim = ", mesh.dimension())

    if mesh.nodeCount() < 2:
        raise Exception("drawMeshBoundaries(ax, mesh): to few nodes",
                        mesh.nodeCount())

    if kwargs.pop('fitView', True):
        ax.set_xlim(mesh.xmin() - 0.05, mesh.xmax() + 0.05)
        ax.set_ylim(mesh.ymin() - 0.05, mesh.ymax() + 0.05)

#    drawAA = True
#    swatch = pg.Stopwatch(True)
    mesh.createNeighbourInfos()

    if not hideMesh:
        linewidth = kwargs.pop('linewidth', 0.3)
        drawSelectedMeshBoundaries(ax, mesh.findBoundaryByMarker(0),
                                   color=(0.0, 0.0, 0.0, 1.0),
                                   linewidth=linewidth)

    drawSelectedMeshBoundaries(ax, mesh.findBoundaryByMarker(
                                     pg.MARKER_BOUND_HOMOGEN_NEUMANN),
                               color=(0.0, 1.0, 0.0, 1.0), linewidth=1.0)
    drawSelectedMeshBoundaries(ax, mesh.findBoundaryByMarker(
                                     pg.MARKER_BOUND_MIXED),
                               color=(1.0, 0.0, 0.0, 1.0), linewidth=1.0)

    b0 = [b for b in mesh.boundaries() if b.marker() > 0]
    if useColorMap:
        drawSelectedMeshBoundaries(ax, b0, color=None,
                                   linewidth=2)
    else:
        drawSelectedMeshBoundaries(ax, b0, color=(0.0, 0.0, 0.0, 1.0),
                                   linewidth=2)

    b4 = [b for b in mesh.boundaries() if b.marker() < -4]
    drawSelectedMeshBoundaries(ax, b4, color=(0.0, 0.0, 0.0, 1.0),
                               linewidth=1.5)

    if mesh.cellCount() == 0:
        drawPLC(ax, mesh, **kwargs)

    updateAxes_(ax)


def drawPLC(ax, mesh, fillRegion=True, regionMarker=True, boundaryMarker=False, **kwargs):
    """Draw 2D PLC into the given ax.

    Parameters
    ----------

    fillRegion: bool [True]
        Fill the regions with default colormap.

    regionMarker: bool [True]
        show region marker

    boundaryMarker: bool [False]
        show boundary marker

    **kwargs

    Examples
    --------
    """
#    eCircles = []
    cols = []

    if fillRegion and mesh.boundaryCount() > 0:
        tmpMesh = pg.meshtools.createMesh(mesh, quality=20, area=0)
        if tmpMesh.cellCount() == 0:
            pass
        else:
            drawModel(ax=ax, mesh=tmpMesh, data=tmpMesh.cellMarkers(),
                      nLevs=len(pg.unique(pg.sort(tmpMesh.cellMarkers()))),
                      levels=pg.utils.unique(tmpMesh.cellMarkers()),
                      tri=True, alpha=0.5, linewidth=0.0, edgecolors='k',
                      snap=True)

    for n in mesh.nodes():
        col = (0.0, 0.0, 0.0, 0.5)

        if n.marker() == pg.MARKER_NODE_SENSOR:
            col = (0.0, 0.0, 0.0, 1.0)

#        eCircles.append(mpl.patches.Circle((n.pos()[0], n.pos()[1])))
        ms = kwargs.pop('markersize', 5)
        # ax.plot(n.pos()[0], n.pos()[1], 'bo', markersize=ms,
                # color=col)
#        eCircles.append(mpl.patches.Circle((n.pos()[0], n.pos()[1]), 0.1))
        cols.append(col)

    if boundaryMarker:
        for b in mesh.boundaries():
            ax.text(b.center()[0], b.center()[1], str(b.marker()),
                    color='red', verticalalignment='center',
                    horizontalalignment='center')  # 'white'

#    p = mpl.collections.PatchCollection(eCircles, color=cols)
#    ax.add_collection(p)

    if regionMarker:
        for reg in mesh.regionMarker():
            ax.text(reg[0], reg[1],
                    str(reg.marker()) + ": " + str(reg.area()),
                    color='black',
                    verticalalignment='center',
                    horizontalalignment='left')  # 'white'

        for hole in mesh.holeMarker():
            ax.text(hole[0], hole[1], 'H', color='black')

    updateAxes_(ax)


def createMeshPatches(ax, mesh, verbose=True, **kwargs):
    """Utility function to create 2d mesh patches within a given ax."""
    if not mesh:
        print("drawMeshBoundaries(ax, mesh): invalid mesh")
        return

    if mesh.nodeCount() < 2:
        print("drawMeshBoundaries(ax, mesh): to few nodes")
        return

    swatch = pg.Stopwatch(True)

    if kwargs.pop('fitView', True):
        ax.set_xlim(mesh.xmin(), mesh.xmax())
        ax.set_ylim(mesh.ymin(), mesh.ymax())

    polys = []

    for cell in mesh.cells():
        if cell.shape().nodeCount() == 3:
            polys.append(list(zip(
                [cell.node(0).x(), cell.node(1).x(), cell.node(2).x()],
                [cell.node(0).y(), cell.node(1).y(), cell.node(2).y()])))
        elif cell.shape().nodeCount() == 4:
            polys.append(list(zip([cell.node(0).x(), cell.node(1).x(),
                                   cell.node(2).x(), cell.node(3).x()],
                                  [cell.node(0).y(), cell.node(1).y(),
                                   cell.node(2).y(), cell.node(3).y()])))
        else:
            print(("unknown shape to patch: ", cell.shape(),
                   cell.shape().nodeCount()))

    patches = mpl.collections.PolyCollection(polys, antialiaseds=False,
                                             picker=True)  # ,lod=True

#    patches.set_edgecolor(None)
    patches.set_edgecolor('face')
#    patches.set_linewidth(1.001)
    ax.add_collection(patches)

    updateAxes_(ax)

    if verbose:
        print(("plotting time = ", swatch.duration(True)))
    return patches


def createTriangles(mesh, data=None):
    """TODO Documentme."""
    x = pg.x(mesh.positions())
#    x.round(1e-1)
    y = pg.y(mesh.positions())
#    y.round(1e-1)

    triCount = 0

    for c in mesh.cells():
        if c.shape().nodeCount() == 4:
            triCount = triCount + 2
        else:
            triCount = triCount + 1

    triangles = np.zeros((triCount, 3))
    dataIdx = list(range(triCount))

    triCount = 0
    for c in mesh.cells():
        if c.shape().nodeCount() == 4:
            triangles[triCount, 0] = c.node(0).id()
            triangles[triCount, 1] = c.node(1).id()
            triangles[triCount, 2] = c.node(2).id()
            dataIdx[triCount] = c.id()
            triCount = triCount + 1

            triangles[triCount, 0] = c.node(0).id()
            triangles[triCount, 1] = c.node(2).id()
            triangles[triCount, 2] = c.node(3).id()
            dataIdx[triCount] = c.id()
            triCount = triCount + 1
        else:
            triangles[triCount, 0] = c.node(0).id()
            triangles[triCount, 1] = c.node(1).id()
            triangles[triCount, 2] = c.node(2).id()
            dataIdx[triCount] = c.id()
            triCount = triCount + 1

    z = None
    if data is not None:
        if len(data) == mesh.cellCount():
            # strange behavior if we just use these slice
            z = np.array(data[dataIdx])
        else:
            z = np.array(data)

    return x, y, triangles, z, dataIdx


def drawMPLTri(ax, mesh, data=None, cMin=None, cMax=None,
               cmap=None, interpolate=False, omitLines=False, **kwargs):
    """Only for triangle/quadrangle meshes currently."""
    x, y, triangles, z, _ = createTriangles(mesh, data)

    gci = None

    levels = kwargs.pop('levels', [])
    nLevs = kwargs.pop('nLevs', 8)
    if len(levels) == 0:
        levels = autolevel(data, nLevs)

    if interpolate and len(data) == mesh.cellCount():
        z = pg.cellDataToPointData(mesh, data)

    if len(z) == len(triangles):
        shading = kwargs.pop('shading', 'flat')
        if shading == 'gouraud':
            z = pg.cellDataToPointData(mesh, data)

        gci = ax.tripcolor(x, y, triangles, facecolors=z, shading=shading,
                           **kwargs)

    elif len(z) == mesh.nodeCount():
        shading = kwargs.pop('shading', None)

        if shading is not None:
            gci = ax.tripcolor(x, y, triangles, z,
                               shading=shading, **kwargs)
        else:
            gci = ax.tricontourf(x, y, triangles, z, levels=levels,
                                 **kwargs)
            if not omitLines:
                ax.tricontour(x, y, triangles, z, levels=levels,
                              colors=kwargs.pop('colors', ['0.5']),
                              **kwargs)
    else:
        gci = None
        raise Exception("Data size does not fit mesh size: ",
                        len(z), mesh.cellCount(), mesh.nodeCount())

    if gci and cMin and cMax:
        gci.set_clim(cMin, cMax)

    if cmap is not None:
        if cmap == 'b2r':
            gci.set_cmap(cmapFromName('b2r'))
        else:
            gci.set_cmap(cmap)

    ax.set_aspect('equal')

    if kwargs.pop('fitView', True):
        ax.set_xlim(mesh.xmin(), mesh.xmax())
        ax.set_ylim(mesh.ymin(), mesh.ymax())

    updateAxes_(ax)
    return gci


def drawField(ax, mesh, data=None, omitLines=False, cmap=None,
              **kwargs):
    """TODO WRITEME.

        Only for triangle/quadrangle meshes currently

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pygimli as pg
    >>> from pygimli.mplviewer import drawField
    >>> n = np.linspace(0, -2, 11)
    >>> mesh = pg.createGrid(x=n, y=n)
    >>> nx = pg.x(mesh.positions())
    >>> ny = pg.y(mesh.positions())
    >>> data = np.cos(1.5 * nx) * np.sin(1.5 * ny)
    >>> fig, ax = plt.subplots()
    >>> drawField(ax, mesh, data)
    <matplotlib.tri.tricontour.TriContourSet ...>
    """
    cMin = kwargs.pop('cMin', None)
    cMax = kwargs.pop('cMax', None)

    return drawMPLTri(ax, mesh, data, cMin=cMin, cMax=cMax,
                      omitLines=omitLines,
                      cmap=cmap, **kwargs)


def drawStreamLines(ax, mesh, u, nx=25, ny=25, **kwargs):
    """Draw streamlines for the gradients of field values u on a mesh.

    The matplotlib routine streamplot needs equidistant spacings so
    we interpolate first on a grid defined by nx and ny nodes.
    Additionally arguments are piped to streamplot.

    This works only for rectangular regions.
    drawStreamLine is more comfortable and more flexible.
    """
    X, Y = np.meshgrid(np.linspace(mesh.xmin(), mesh.xmax(), nx),
                       np.linspace(mesh.ymin(), mesh.ymax(), ny))

    U = X.copy()
    V = X.copy()

    for i, row in enumerate(X):
        for j in range(len(row)):
            p = [X[i, j], Y[i, j]]
            gr = [0.0, 0.0]
            c = mesh.findCell(p)
            if c:
                gr = c.grad(p, u)

            U[i, j] = -gr[0]
            V[i, j] = -gr[1]

    gci = ax.streamplot(X, Y, U, V, **kwargs)

    updateAxes_(ax)
    return gci
# def drawStreamLines(...)


def drawStreamLine_(ax, mesh, c, data, dataMesh=None,
                    linewidth=1.0, dropTol=0.0, **kwargs):
    """Draw a single streamline.

    Draw a single streamline into a given mesh for given data stating at
    the center of cell c.
    The Streamline will be enlarged until she reached a cell that
    already contains a streamline.

    TODO
        linewidth and color depends on absolute velocity
        or background color saturation

    Parameters
    ----------

    ax : matplotlib.ax
        ax to draw into

    mesh : :gimliapi:`GIMLI::Mesh`
        2d Mesh to draw the streamline

    c : :gimliapi:`GIMLI::Cell`
        start cell

    data : iterable float | [float, float]
        If data is an array (per cell or node) gradients are calculated
        otherwise the data will be interpreted as vector field.

    dataMesh : :gimliapi:`GIMLI::Mesh` [None]

        Optional mesh for the data. If you want high resolution
        data to plot on coarse draw mesh.

    linewidth : float [1.0]

        Streamline linewidth

    dropTol : float [0.0]

        Don't draw stream lines with velocity lower than drop tolerance.

    """
    x, y, v = streamline(mesh, data, startCoord=c.center(),
                         dLengthSteps=5,
                         dataMesh=dataMesh,
                         maxSteps=10000,
                         verbose=False,
                         koords=[0, 1])

    if 'color' not in kwargs:
        kwargs['color'] = 'black'

    lines = None

    if len(x) > 2:


        points = np.array([x, y]).T.reshape(-1, 1, 2)

        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        lwidths = pg.RVector(len(v), linewidth)
        lwidths[pg.find(pg.RVector(v) < dropTol)] = 0.0

        lines = LineCollection(segments, linewidths=lwidths, **kwargs)
        ax.add_collection(lines)

        # probably the limits are wrong without plot call
        #lines = ax.plot(x, y, **kwargs)

        #updateAxes_(ax, lines)

        #ax.plot(x, y, '.-', color='black', **kwargs)
    if len(x) > 3:
        xmid = int(len(x) / 2)
        ymid = int(len(y) / 2)
        dx = x[xmid + 1] - x[xmid]
        dy = y[ymid + 1] - y[ymid]
        c = mesh.findCell([x[xmid], y[ymid]])
        dLength = c.center().dist(c.node(0).pos()) / 4.

        if v[xmid] > dropTol:
            ax.arrow(x[xmid], y[ymid], dx, dy, width=dLength / 15.,
                    head_starts_at_zero=True,
                    **kwargs)

    return lines


def drawStreams(ax, mesh, data, startStream=3, **kwargs):
    """Draw streamlines based on an unstructured mesh.

    Every cell contains only one streamline and every new stream line
    starts in the center of a cell. You can alternatively provide a second mesh
    with coarser mesh to draw streams for.

    Parameters
    ----------

    ax : matplotlib.ax
        ax to draw into

    mesh : :gimliapi:`GIMLI::Mesh`
        2d Mesh to draw the streamline

    data : iterable float | [float, float]
        If data is an array (per cell or node) gradients are calculated
        otherwise the data will be interpreted as vector field.

    startStream : int
        variate the start stream drawing, try values from 1 to 3 what every
        you like more.

    **kwargs: forward to drawStreamLine_

        * coarseMesh

            Instead of draw a stream for every cell in mesh, draw a streamline
            segment for each cell in coarseMesh.

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pygimli as pg
    >>> from pygimli.mplviewer import drawStreams
    >>> n = np.linspace(1, 2, 10)
    >>> mesh = pg.createGrid(x=n, y=n)
    >>> nx = pg.x(mesh.positions())
    >>> ny = pg.y(mesh.positions())
    >>> data = np.cos(1.5 * nx) * np.sin(1.5 * ny)
    >>> fig, ax = plt.subplots()
    >>> drawStreams(ax, mesh, data, color='red')
    >>> drawStreams(ax, mesh, data, dropTol=0.9)
    >>> ax.set_aspect('equal')
    """
    viewMesh = None
    dataMesh = None

    if 'coarseMesh' in kwargs:
        viewMesh = kwargs['coarseMesh']
        dataMesh = mesh
        dataMesh.createNeighbourInfos()
        del kwargs['coarseMesh']
    else:
        viewMesh = mesh

    viewMesh.createNeighbourInfos()

    for c in viewMesh.cells():
        c.setValid(True)

    if startStream == 1:
        # start a stream from each boundary cell
        for y in np.linspace(viewMesh.ymin(), viewMesh.ymax(), 100):
            c = viewMesh.findCell(
                [(viewMesh.xmax() - viewMesh.xmax()) / 2.0, y])
            if c is not None:
                if c.valid():
                    drawStreamLine_(ax, viewMesh, c, data, dataMesh,
                                    **kwargs)

    elif startStream == 2:
        # start a stream from each boundary cell
        for x in np.linspace(viewMesh.xmin(), viewMesh.xmax(), 100):
            c = viewMesh.findCell(
                [x, (viewMesh.ymax() - viewMesh.ymax()) / 2.0])
            if c is not None:
                if c.valid():
                    drawStreamLine_(ax, viewMesh, c, data, dataMesh,
                                    **kwargs)

    elif startStream == 3:
        # start a stream from each boundary cell
        for b in viewMesh.findBoundaryByMarker(1, 99):
            c = b.leftCell()
            if c is None:
                c = b.rightCell()

            if c.valid():
                drawStreamLine_(ax, viewMesh, c, data, dataMesh,
                                **kwargs)

    # start a stream from each unused cell
    for c in viewMesh.cells():
        if c.valid():
            drawStreamLine_(ax, viewMesh, c, data, dataMesh,
                            **kwargs)

    for c in viewMesh.cells():
        c.setValid(True)

    updateAxes_(ax)
# def drawStreamLines2(...)


def drawSensors(ax, sensors, diam=None, koords=None, verbose=False, **kwargs):
    """Draw sensor positions as black dots with a given diameter.

    Parameters
    ----------

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from pygimli.mplviewer import drawSensors
    >>> sensors = np.random.rand(5, 2)
    >>> fig, ax = plt.subplots()
    >>> drawSensors(ax, sensors, diam=0.02, koords=[0, 1])
    >>> ax.set_aspect('equal')
    """
    if koords is None:
        koords = [0, 2]
        if pg.yVari(sensors):
            koords = [0, 1]

    eCircles = []

    if diam is None:
        eSpacing = sensors[0].distance(sensors[1])
        diam = eSpacing / 8.0

    for e in sensors:
        if verbose:
            print(e, diam, e[koords[0]], e[koords[1]])
        eCircles.append(mpl.patches.Circle((e[koords[0]], e[koords[1]]), diam))

    p = mpl.collections.PatchCollection(eCircles,
                                        **kwargs)
    p.set_zorder(100)
    ax.add_collection(p)

    updateAxes_(ax)


def createParameterContraintsLines(mesh, cMat, cWeight=None):
    """TODO Documentme."""
    C = pg.RMatrix()
    if isinstance(cMat, pg.DSparseMapMatrix):
        cMat.save('tmpC.matrix')
        pg.loadMatrixCol(C, 'tmpC.matrix')
    else:
        C = cMat

    paraMarker = mesh.cellMarkers()
    cellList = dict()

    for cID in range(len(paraMarker)):
        if cID not in cellList:
            cellList[cID] = []
        cellList[cID].append(mesh.cell(cID))

    paraCenter = dict()
    for cID, vals in list(cellList.items()):
        p = pg.RVector3(0.0, 0.0, 0.0)
        for c in vals:
            p += c.center()
        p /= float(len(vals))
        paraCenter[cID] = p

    nConstraints = C[0].size()
    start = []
    end = []
#    swatch = pg.Stopwatch(True)  # not used
    for i in range(0, nConstraints / 2):
        # print i
        # if i == 1000: break;
        idL = int(C[1][i * 2])
        idR = int(C[1][i * 2 + 1])
        # leftCells = []
        # rightCells = []
#        for c, index in enumerate(paraMarker):
#            if idL == index:
#                leftCells.append(mesh.cell(c))
#            if idR == index:
#                rightCells.append(mesh.cell(c))

#        p1 = pg.RVector3(0.0,0.0);
#        for c in leftCells:
#            p1 += c.center()
#        p1 /= float(len(leftCells))

#        p2 = pg.RVector3(0.0,0.0);
#        for c in rightCells:
#            p2 += c.center()
#        print cWeight[i]
#        p2 /= float(len(rightCells))
        p1 = paraCenter[idL]
        p2 = paraCenter[idR]

        if cWeight is not None:
            pa = pg.RVector3(p1 + (p2 - p1) / 2.0 * (1.0 - cWeight[i]))
            pb = pg.RVector3(p2 + (p1 - p2) / 2.0 * (1.0 - cWeight[i]))
        else:
            pa = p1
            pb = p2

        start.append(pa)
        end.append(pb)

#    updateAxes_(ax)  # not existing

    return start, end


def drawParameterConstraints(ax, mesh, cMat, cWeight=None):
    """TODO Documentme."""
    start, end = createParameterContraintsLines(mesh, cMat, cWeight)

    lines = []
    colors = []
    linewidths = []
    for i, _ in enumerate(start):
        lines.append(list(zip([start[i].x(), end[i].x()],
                              [start[i].y(), end[i].y()])))

        linewidth = 0.5
        col = (0.0, 0.0, 1.0, 1.0)
        colors.append(col)
        linewidths.append(linewidth)

    linCol = mpl.collections.LineCollection(lines, antialiaseds=True)

    linCol.set_color(colors)
    linCol.set_linewidth(linewidths)
    ax.add_collection(linCol)

    updateAxes_(ax)


def draw1DColumn(ax, x, val, thk, width=30, ztopo=0, cmin=1, cmax=1000,
                 cmap=None, name=None, textoffset=0.0):
    """Draw a 1D column (e.g., from a 1D inversion) on a given ax.

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from pygimli.mplviewer import draw1DColumn
    >>> thk = [1, 2, 3, 4]
    >>> val = thk
    >>> fig, ax = plt.subplots()
    >>> draw1DColumn(ax, 0.5, val, thk, width=0.1, cmin=1, cmax=4, name="VES")
    <matplotlib.collections.PatchCollection object at ...>
    >>> ax.set_ylim(-np.sum(thk), 0)
    (-10, 0)
    """
    z = -np.hstack((0., np.cumsum(thk), np.sum(thk) * 1.5)) + ztopo
    recs = []
    for i in range(len(val)):
        recs.append(Rectangle((x - width / 2., z[i]), width, z[i + 1] - z[i]))

    pp = PatchCollection(recs)
    col = ax.add_collection(pp)

    pp.set_edgecolor(None)
    pp.set_linewidths(0.0)

    if cmap is not None:
        if isinstance(cmap, str):
            pp.set_cmap(cmapFromName(cmap))
        else:
            pp.set_cmap(cmap)

    pp.set_norm(LogNorm(cmin, cmax))
    pp.set_array(np.array(val))
    pp.set_clim(cmin, cmax)
    if name:
        ax.text(x+textoffset, ztopo, name, ha='center', va='bottom')

    updateAxes_(ax)

    return col


def insertUnitAtNextLastTick(ax, unit, xlabel=True, position=-2):
    """Replace the last-but-one tick label by unit symbol."""
    if xlabel:
        labels = ax.get_xticks().tolist()
        labels[position] = unit
        ax.set_xticklabels(labels)
    else:
        labels = ax.get_yticks().tolist()
        labels[position] = unit
        ax.set_yticklabels(labels)


def plotLines(ax, line_filename, linewidth=1.0, step=1):
    """Read lines from file and plot over model."""
    xz = np.loadtxt(line_filename)
    n_points = xz.shape[0]
    if step == 2:
        for i in range(0, n_points, step):
            x = xz[i:i+step, 0]
            z = xz[i:i+step, 1]
            ax.plot(x, z, 'k-', linewidth=linewidth)
    if step == 1:
        ax.plot(xz[:, 0], xz[:, 1], 'k-', linewidth=linewidth)
