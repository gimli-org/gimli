# -*- coding: utf-8 -*-
"""Generic mesh visualization tools."""

import os
import sys
import time
import traceback

import numpy as np
# plt should not be used outside of viewer.mpl
import matplotlib.pyplot as plt

from .. core.logger import renameKwarg

try:
    import pygimli as pg
    from .mpl import drawMesh, drawModel, drawField
    from .showmatrix import showMatrix
    from .mpl import drawSensors
    from .mpl import createColorBar, updateColorBar
    from .mpl import drawStreams, addCoverageAlpha
    from .mpl import CellBrowser
    from .mpl.colorbar import cmapFromName
except ImportError as e:
    print(e)
    traceback.print_exc(file=sys.stdout)
    pg.critical("ERROR: cannot import the library 'pygimli'."
                "Ensure that pygimli is in your PYTHONPATH ")


def show(obj=None, data=None, **kwargs):
    """Mesh and model visualization.

    Syntactic sugar to show a obj with data. Forwards to
    a known visualization for obj. Typical is
    :py:mod:`pygimli.viewer.showMesh` or
    :py:mod:`pygimli.viewer.mayaview.showMesh3D` to show most of the typical 2D
    and 3D content.
    See tutorials and examples for usage hints. An empty show
    call creates an empty ax window.

    Parameters
    ----------
    obj: obj
        obj can be so far.
        * :gimliapi:`GIMLI::Mesh` or list of meshes
        * DataContainer
        * pg.core.Sparse[Map]Matrix

    data: iterable
        Optionally data to visualize. See appropriate show function.

    Keyword Arguments
    -----------------
    **kwargs
        Additional kwargs forward to appropriate show functions.

        * ax : axe [None]
            Matplotlib axes object. Create a new if necessary.
        * fitView : bool [True]
            Scale x and y limits to match the view.

    Returns
    -------
    Return the results from the showMesh* functions. Usually the axe object
    and a colorbar.

    See Also
    --------
    showMesh
    """
    if "axes" in kwargs: # remove me in 1.2 #20200515
        print("Deprecation Warning: Please use keyword `ax` instead of `axes`")
        kwargs['ax'] = kwargs.pop('axes', None)

    ### Empty call just to create a axes
    if obj is None and not 'mesh' in kwargs.keys():
        ax = kwargs.pop('ax', None)

        if ax is None:
            ax = plt.subplots()[1]
        return ax, None

    ### try to interprete obj containes a mesh
    if hasattr(obj, 'mesh'):
        return pg.show(obj.mesh, obj, **kwargs)

    ### try to interprete obj as ERT Data
    if isinstance(obj, pg.DataContainerERT):
        from pygimli.physics.ert import showERTData
        return showERTData(obj, vals=kwargs.pop('vals', data), **kwargs)

    ### try to interprete obj as matrices
    if isinstance(obj, pg.core.MatrixBase) or \
        (isinstance(obj, np.ndarray) and obj.ndim == 2):
        return showMatrix(obj, **kwargs)

    try:
        from scipy.sparse import spmatrix
        if isinstance(obj, spmatrix):
            return showMatrix(obj, **kwargs)
    except ImportError:
        pass

    ### try to interprete obj as mesh or list of meshes
    mesh = kwargs.pop('mesh', obj)

    if isinstance(mesh, list):
        ax = kwargs.pop('ax', None)
        fitView = kwargs.pop('fitView', ax is None)

        ax, cBar = show(mesh[0], data, hold=1, ax=ax, fitView=fitView, **kwargs)
        xMin = mesh[0].xMin()
        xMax = mesh[0].xMax()
        yMin = mesh[0].yMin()
        yMax = mesh[0].yMax()

        for m in mesh[1:]:
            ax, cBar = show(m, data, ax=ax, hold=1, fitView=False, **kwargs)
            xMin = min(xMin, m.xMin())
            xMax = max(xMax, m.xMax())
            yMin = min(yMin, m.yMin())
            yMax = max(yMax, m.yMax())

#        ax.relim()
#        ax.autoscale_view(tight=True)
        if fitView is not False:
            ax.set_xlim([xMin, xMax])
            ax.set_ylim([yMin, yMax])
        #        print(ax.get_data_interval())
        return ax, cBar

    if isinstance(mesh, pg.Mesh):
        if mesh.dim() == 2:
            if pg.zero(pg.y(mesh)):
                pg.info("swap z<->y coordinates for visualization.")
                meshSwap = pg.Mesh(mesh)
                for n in meshSwap.nodes():
                    n.pos()[1] = n.pos()[2]
                return showMesh(meshSwap, data, **kwargs)

            return showMesh(mesh, data, **kwargs)
        elif mesh.dim() == 3:

            from .vistaview import showMesh3D
            return showMesh3D(mesh, data, **kwargs)
        else:
            pg.error("ERROR: Mesh not valid.", mesh)

    pg.error("Can't interprete obj: {0} to show.".format(obj))
    return None, None


def showMesh(mesh, data=None, hold=False, block=False, colorBar=None,
             label=None, coverage=None, ax=None, savefig=None,
             showMesh=False, showBoundary=None,
             markers=False, **kwargs):
    """2D Mesh visualization.

    Create an axis object and plot a 2D mesh with given node or cell data.
    Returns the axis and the color bar. The type of data determines the
    appropriate draw method.

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`
        2D or 3D GIMLi mesh
    data: iterable [None]
        Optionally data to visualize.

        . None (draw mesh only)
            forward to :py:mod:`pygimli.viewer.mpl.drawMesh`
            or if no cells are given:
            forward to :py:mod:`pygimli.viewer.mpl.drawPLC`

        . [[marker, value], ...]
            List of Cellvalues per cell marker
            forward to :py:mod:`pygimli.viewer.mpl.drawModel`

        . float per cell -- model, patch
            forward to :py:mod:`pygimli.viewer.mpl.drawModel`

        . float per node -- scalar field
            forward to :py:mod:`pygimli.viewer.mpl.drawField`

        . iterable of type [float, float] -- vector field
            forward to :py:mod:`pygimli.viewer.mpl.drawStreams`

        . pg.core.R3Vector -- vector field
            forward to :py:mod:`pygimli.viewer.mpl.drawStreams`

        . pg.core.stdVectorRVector3 -- sensor positions
            forward to :py:mod:`pygimli.viewer.mpl.drawSensors`
    hold: bool [false]
        Set interactive plot mode for matplotlib.
        If this is set to false [default] your script will open
        a window with the figure and draw your content.
        If set to true nothing happens until you either force another show with
        hold=False, you call plt.show() or pg.wait().
        If you want show with stopping your script set block = True.
    block: bool [false]
        Force show drawing your content and block the script until you
        close the current figure.
    colorBar: bool [None], Colorbar
        Create and show a colorbar. If colorBar is a valid colorbar then only
        its values will be updated.
    label: str
        Set colorbar label. If set colorbar is toggled to True. [None]
    coverage: iterable [None]
        Weight data by the given coverage array and fadeout the color.
    ax: matplotlib.Axes [None]
        Instead of creating a new and empty ax, just draw into the given one.
        Useful to combine multiple plots into one figure.
    savefig: string
        Filename for a direct save to disc.
        The matplotlib pdf-output is a little bit big so we try
        an epstopdf if the .eps suffix is found in savefig
    showMesh: bool [False]
        Shows the mesh itself additional.
    showBoundary: bool [None]
        Shows all boundary with marker != 0. A value None means automatic
        True for cell data and False for node data.
    marker: bool [False]
        Show mesh and boundary marker.

    Keyword Arguments
    -----------------
    **kwargs:
        * xlabel: str [None]
            Add label to the x axis
        * ylabel: str [None]
            Add label to the y axis
        fitView: bool
            Fit the axes limits to the view object. Default is True if ax is None else is set to False.
        All remaining will be forwarded to the draw functions
        and matplotlib methods, respectively.

    Examples
    --------
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> world = mt.createWorld(start=[-10, 0], end=[10, -10],
    ...                        layers=[-3, -7], worldMarker=False)
    >>> mesh = mt.createMesh(world, quality=32, area=0.2, smooth=[1, 10])
    >>> _ = pg.viewer.showMesh(mesh, markers=True)

    Returns
    -------
    ax : matplotlib.axes

    cBar : matplotlib.colorbar
    """
    renameKwarg('cmap', 'cMap', kwargs)

    cMap = kwargs.pop('cMap', 'viridis')
    nCols = None
    cBarOrientation = kwargs.pop('orientation', 'horizontal')

    if ax is None:
        ax = plt.subplots()[1]

    # adjust limits only when axis is empty
    if (ax.lines or ax.collections or ax.patches):
        fitViewDefault = False
    else:
        fitViewDefault = True


    # plt.subplots() resets locale setting to system default .. this went
    # horrible wrong for german 'decimal_point': ','
    pg.checkAndFixLocaleDecimal_point(verbose=False)

    if block:
        hold = True

    lastHoldStatus = pg.viewer.mpl.utils.holdAxes__
    if not lastHoldStatus or hold:
        pg.viewer.mpl.hold(val=1)
        hold = True

    gci = None
    validData = False

    if markers:
        kwargs["boundaryMarker"] = True
        if mesh.cellCount() > 0:
            uniquemarkers, uniqueidx = np.unique(
                np.array(mesh.cellMarkers()), return_inverse=True)
            label = "Cell markers"
            cMap = plt.cm.get_cmap("Set3", len(uniquemarkers))
            kwargs["logScale"] = False
            kwargs["cMin"] = -0.5
            kwargs["cMax"] = len(uniquemarkers) - 0.5
            data = np.arange(len(uniquemarkers))[uniqueidx]

    if data is None:
        showMesh = True
        mesh.createNeighborInfos()
        if showBoundary is None:
            showBoundary = True
    elif isinstance(data, pg.core.stdVectorRVector3):
        drawSensors(ax, data, **kwargs)
    elif isinstance(data, pg.core.R3Vector):
        drawStreams(ax, mesh, data, **kwargs)
    else:
        ### data=[[marker, val], ....]
        if isinstance(data, list) and \
            isinstance(data[0], list) and isinstance(data[0][0], int):
            data = pg.solver.parseMapToCellArray(data, mesh)

        if hasattr(data[0], '__len__') and not \
            isinstance(data, np.ma.core.MaskedArray):

            if len(data) == 2:  # [u,v] x N
                data = np.array(data).T

            if data.shape[1] == 2:
                drawStreams(ax, mesh, data, **kwargs)

            elif data.shape[1] == 3:  # probably N x [u,v,w]
                # if sum(data[:, 0]) != sum(data[:, 1]):
                # drawStreams(ax, mesh, data, **kwargs)
                drawStreams(ax, mesh, data[:, 0:2], **kwargs)
            else:
                pg.warn("No valid stream data:", data.shape, data.ndim)
                showMesh = True
        # elif min(data) == max(data):  # or pg.core.haveInfNaN(data):
        #     pg.warn("No valid data: ", min(data), max(data), pg.core.haveInfNaN(data))
        #     showMesh = True
        else:
            validData = True
            if bool(colorBar) is not False:
                colorBar = True

            try:
                if len(data) == mesh.cellCount():
                    kwargs['nCols'] = kwargs.pop('nCols', 256)
                    if label is None:
                        label = ""

                    gci = drawModel(ax, mesh, data, **kwargs)
                    if showBoundary is None:
                        showBoundary = True

                elif len(data) == mesh.nodeCount():
                    kwargs['nLevs'] = kwargs.pop('nLevs', 5)
                    kwargs['nCols'] = kwargs.pop('nCols', kwargs['nLevs']-1)
                    if label is None:
                        label = ""

                    gci = drawField(ax, mesh, data, **kwargs)
                else:
                    pg.error("Data size invalid")
                    print("Data: ", len(data), min(data), max(data), pg.core.haveInfNaN(data))
                    print("Mesh: ", mesh)
                    validData = False
                    drawMesh(ax, mesh)

                if cMap is not None and gci is not None:
                    gci.set_cmap(cmapFromName(cMap))
                    #gci.cmap.set_under('k')

            except BaseException as e:
                pg.error("Exception occurred: ", e)

    if mesh.cellCount() == 0:
        showMesh = False
        if mesh.boundaryCount() == 0:
            pg.viewer.mpl.drawPLC(ax, mesh, showNodes=True,
                                 fillRegion=False, showBoundary=False,
                                 **kwargs)
            showBoundary = False
            #ax.plot(pg.x(mesh), pg.y(mesh), '.', color='black')
        else:
            pg.viewer.mpl.drawPLC(ax, mesh, **kwargs)

    if showMesh:
        if gci is not None and hasattr(gci, 'set_antialiased'):
            gci.set_antialiased(True)
            gci.set_linewidth(0.3)
            gci.set_edgecolor("0.1")
        else:
            pg.viewer.mpl.drawSelectedMeshBoundaries(ax, mesh.boundaries(),
                                            color=kwargs.pop('color', "0.1"), linewidth=0.3)
            #drawMesh(ax, mesh, **kwargs)

    if showBoundary == True or showBoundary == 1:
        b = mesh.boundaries(mesh.boundaryMarkers() != 0)
        pg.viewer.mpl.drawSelectedMeshBoundaries(ax, b,
                                                color=(0.0, 0.0, 0.0, 1.0),
                                                linewidth=1.4)

    fitView = kwargs.pop('fitView', fitViewDefault)
    if fitView:
        ax.set_xlim(mesh.xMin(), mesh.xMax())
        ax.set_ylim(mesh.yMin(), mesh.yMax())
        ax.set_aspect('equal')

    cBar = None

    if label is not None and colorBar is None:
        colorBar = True

    if colorBar and validData:
        labels = ['cMin', 'cMax', 'nCols', 'nLevs', 'logScale', 'levels']
        subkwargs = {key: kwargs[key] for key in labels if key in kwargs}

        subkwargs['label'] = label
        subkwargs['cMap'] = cMap
        subkwargs['orientation'] = cBarOrientation

        if bool(colorBar):
            cBar = createColorBar(gci,
                                  size=kwargs.pop('size', 0.2),
                                  pad=kwargs.pop('pad', None),
                                  **subkwargs
                                  )
        elif colorBar is not False:
            cBar = updateColorBar(colorBar, **subkwargs)

        if markers:
            ticks = np.arange(len(uniquemarkers))
            cBar.set_ticks(ticks)
            labels = []
            for marker in uniquemarkers:
                labels.append(str((marker)))
            cBar.set_ticklabels(labels)

    if coverage is not None:
        if len(data) == mesh.cellCount():
            addCoverageAlpha(gci, coverage,
                             dropThreshold=kwargs.pop('dropThreshold', 0.4))
        else:
            pg.error('show, coverage wrong length, toImplement')
            # addCoverageAlpha(gci, pg.core.cellDataToPointData(mesh, coverage))

    if not hold or block is not False and plt.get_backend().lower() != "agg":
        if data is not None:
            if len(data) == mesh.cellCount():
                CellBrowser(mesh, data, ax=ax)

        plt.show(block=block)
        try:
            plt.pause(0.01)
        except BaseException as _:
            pass

    if hold:
        pg.viewer.mpl.hold(val=lastHoldStatus)

    if savefig:
        print('saving: ' + savefig + ' ...')

        if '.' not in savefig:
            savefig += '.pdf'

        ax.figure.savefig(savefig, bbox_inches='tight')
        # rc params savefig.format=pdf

        if '.eps' in savefig:
            try:
                print("trying eps2pdf ... ")
                os.system('epstopdf ' + savefig)
            except BaseException:
                pass
        print('.. done')

    return ax, cBar


def showBoundaryNorm(mesh, normMap=None, **kwargs):
    """Show mesh boundaries normals.

    Show the mesh and draw a black line along the normal direction of all
    boundaries. If you provide a boundary marker vs. norm direction map,
    then only these norms are drawn.

    Parameters
    ----------

    mesh : :gimliapi:`GIMLI::Mesh`
        2D or 3D GIMLi mesh

    normMap : list
        list of [boundary marker, [norm]] pairs. e.g. [[1, [0.0,1.0]], ... ]

    **kwargs :
        Will be forwarded to the draw functions and matplotlib methods,
        respectively.

    Returns
    -------
    ax : matplotlib.ax
    """
    ax = kwargs.pop('ax', None)

    col = kwargs.pop('color', 'Black')

    if normMap:
        for pair in normMap:
            bounds = mesh.findBoundaryByMarker(pair[0])

            for b in bounds:
                c1 = b.center()

                if (pair[1][0] != 0) or (pair[1][1] != 0):
                    ax.arrow(c1[0], c1[1], pair[1][0], pair[1][1],
                             head_width=0.1, head_length=0.3, color=col,
                             **kwargs)
                else:
                    ax.plot(c1[0], c1[1], 'o', color=col)
        return

    ax = show(mesh, hold=True, ax=ax)[0]
    for b in mesh.boundaries():
        c1 = b.center()
        c2 = c1 + b.norm()
        ax.plot([c1[0], c2[0]], [c1[1], c2[1]], color=col, **kwargs)

    time.sleep(0.05)

    return ax
