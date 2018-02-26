# -*- coding: utf-8 -*-
"""Generic mesh visualization tools."""

import os
import sys
import time
import traceback

# plt should not be used outside of mplviewer
import matplotlib.pyplot as plt
import numpy as np

try:
    import pygimli as pg
    from pygimli.mplviewer import drawMesh, drawModel, drawField
    from pygimli.mplviewer import drawSensors
    from pygimli.mplviewer import createColorBar, updateColorBar
    from pygimli.mplviewer import drawStreams, addCoverageAlpha
    from pygimli.mplviewer import CellBrowser
    from pygimli.mplviewer.colorbar import cmapFromName
except ImportError as e:
    print(e)
    traceback.print_exc(file=sys.stdout)
    raise Exception('''ERROR: cannot import the library 'pygimli'.
        Ensure that pygimli is in your PYTHONPATH ''')


def show(mesh=None, data=None, **kwargs):
    """Mesh and model visualization.

    Syntactic sugar to show a mesh with data. Forwards to
    :py:mod:`pygimli.viewer.showMesh` or
    :py:mod:`pygimli.viewer.mayaview.showMesh3D` to show most of the typical 2D
    and 3D content. See tutorials and examples for usage hints. An empty show
    call create an empty ax window.

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh` or list of meshes
        2D or 3D GIMLi mesh

    **kwargs :
        Will be forwarded to the appropriate show functions.

    Returns
    -------
    Return the results from the showMesh* functions.

    See Also
    --------
    showMesh
    """
    if "axes" in kwargs:
        print("Deprecation Warning: Please use keyword `ax` instead of `axes`")
        kwargs['ax'] = kwargs.pop('axes', None)

    if isinstance(mesh, list):
        ax = kwargs.pop('ax', None)
        ax, cbar = show(mesh[0], data, hold=1, ax=ax, **kwargs)
        xmin = mesh[0].xmin()
        xmax = mesh[0].xmax()
        ymin = mesh[0].ymin()
        ymax = mesh[0].ymax()

        for m in mesh[1:]:
            ax, cbar = show(m, data, ax=ax, hold=1, fitView=False, **kwargs)
            xmin = min(xmin, m.xmin())
            xmax = max(xmax, m.xmax())
            ymin = min(ymin, m.ymin())
            ymax = max(ymax, m.ymax())


#        ax.relim()
#        ax.autoscale_view(tight=True)
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
        #        print(ax.get_data_interval())
        return ax, cbar

    if isinstance(mesh, pg.Mesh):
        if mesh.dim() == 2:
            return showMesh(mesh, data, **kwargs)
        elif mesh.dim() == 3:

            from .mayaview import showMesh3D

            return showMesh3D(mesh, data, **kwargs)
        else:
            print("ERROR: Mesh not valid.")

    ax = kwargs.pop('ax', None)

    if ax is None:
        ax = plt.subplots()[1]

    return ax, None


def showMesh(mesh, data=None, hold=False, block=False, colorBar=None,
             label=None, coverage=None, ax=None, savefig=None, showMesh=False,
             markers=False, **kwargs):
    """2D Mesh visualization.

    Create an axis object and plot a 2D mesh with given node or cell data.
    Returns the axis and the color bar. The type of data determine the
    appropriate draw method.

    Parameters
    ----------

    mesh : :gimliapi:`GIMLI::Mesh`
        2D or 3D GIMLi mesh

    data : iterable [None]
        Optionally data to visualize.

        . None (draw mesh only)
            forward to :py:mod:`pygimli.mplviewer.drawMesh`
            or if no cells are given:
            forward to :py:mod:`pygimli.mplviewer.drawPLC`

        . [[marker, value], ...]
            List of Cellvalues per cell marker
            forward to :py:mod:`pygimli.mplviewer.drawModel`

        . float per cell -- model, patch
            forward to :py:mod:`pygimli.mplviewer.drawModel`

        . float per node -- scalar field
            forward to :py:mod:`pygimli.mplviewer.drawField`

        . iterable of type [float, float] -- vector field
            forward to :py:mod:`pygimli.mplviewer.drawStreams`

        . pg.R3Vector -- vector field
            forward to :py:mod:`pygimli.mplviewer.drawStreams`

        . pg.stdVectorRVector3 -- sensor positions
            forward to :py:mod:`pygimli.mplviewer.drawSensors`


    hold : bool [false]
        Set interactive plot mode for matplotlib.
        If this is set to false [default] your script will open
        a window with the figure and draw your content.
        If set to true nothing happens until you either force another show with
        hold=False, you call plt.show() or pg.wait().
        If you want show with stopping your script set block = True.

    block : bool [false]
        Force show drawing your content and block the script until you
        close the current figure.

    colorBar : bool [None], Colorbar
        Create and show a colorbar. If colorBar is a valid colorbar then only
        his values will be updated.

    label : str
        Set colorbar label. If set colorbar is toggled to True. [None]

    coverage : iterable [None]
        Weight data by the given coverage array and fadeout the color.

    ax : matplotlib.Axes [None]
        Instead of create a new and empty ax, just draw into the a given.
        Useful to combine draws.

    savefig: string
        Filename for a direct save to disc.
        The matplotlib pdf-output is a little bit big so we try
        an epstopdf if the .eps suffix is found in savefig

    showMesh : bool [False]
        Shows the the mesh itself aditional.

    marker : bool [False]
        Show mesh and boundary marker.

    **kwargs :
        * xlabel : str [None]
            Add label to the x axis

        * ylabel : str [None]
            Add label to the y axis

        * all remaining
            Will be forwarded to the draw functions and matplotlib methods,
            respectively.

    Returns
    -------
    ax : matplotlib.axes

    colobar : matplotlib.colorbar
    """
    if ax is None:
        ax = plt.subplots()[1]

    # print('1*'*50)
    # print(locale.localeconv())

    # plt.subplots() resets locale setting to system default .. this went
    # horrible wrong for german 'decimal_point': ','
    pg.checkAndFixLocaleDecimal_point(verbose=False)

    # print('2*'*50)
    # print(locale.localeconv())

    if block:
        hold = True

    if hold:
        lastHoldStatus = pg.mplviewer.utils.holdAxes__
        pg.mplviewer.hold(val=1)

    gci = None
    validData = False

    if markers:
        kwargs["boundaryMarker"] = True
        if mesh.cellCount() > 0:
            uniquemarkers, uniqueidx = np.unique(
                np.array(mesh.cellMarkers()), return_inverse=True)
            label = "Cell markers"
            kwargs["cMap"] = plt.cm.get_cmap("Set3", len(uniquemarkers))
            kwargs["logScale"] = False
            kwargs["cMin"] = -0.5
            kwargs["cMax"] = len(uniquemarkers) - 0.5
            kwargs["grid"] = True
            data = np.arange(len(uniquemarkers))[uniqueidx]
    if data is None:
        drawMesh(ax, mesh, **kwargs)

    elif isinstance(data, pg.stdVectorRVector3):
        drawSensors(ax, data, **kwargs)
    elif isinstance(data, pg.R3Vector):
        drawStreams(ax, mesh, data, **kwargs)
    else:
        #print('-----------------------------')
        #print(data, type(data))
        #print('-----------------------------')

        ### data=[[marker, val], ....]
        if type(data) is list and \
            type(data[0]) is list and type(data[0][0]) is int:
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
                print("No valid stream data:", data.shape, data.ndim)
                drawMesh(ax, mesh, **kwargs)
        elif min(data) == max(data):  # or pg.haveInfNaN(data):
            print("No valid data: ", min(data), max(data), pg.haveInfNaN(data))
            drawMesh(ax, mesh, **kwargs)
        else:
            validData = True
            try:
                if len(data) == mesh.cellCount():
                    gci = drawModel(ax, mesh, data, **kwargs)
                elif len(data) == mesh.nodeCount():
                    gci = drawField(ax, mesh, data, **kwargs)

                cmap = kwargs.pop('cmap', None)
                cMap = kwargs.pop('cMap', None)
                if cMap is not None:
                    cmap = cMap

                if cmap is not None:
                    gci.set_cmap(cmapFromName(cmap))

            except BaseException as e:
                print("Exception occured: " + e)
                print("Data: ", min(data), max(data), pg.haveInfNaN(data))
                print("Mesh: ", mesh)
                drawMesh(ax, mesh, **kwargs)

    ax.set_aspect('equal')

    cbar = None

    if label is not None and colorBar is None:
        colorBar = True

    if colorBar and validData:
        # , **kwargs) # causes problems!
        labels = ['cMin', 'cMax', 'nLevs', 'orientation', 'pad']
        subkwargs = {key: kwargs[key] for key in labels if key in kwargs}

        if colorBar is True or colorBar is 1:
            cbar = createColorBar(gci, label=label, **subkwargs)
        elif colorBar is not False:
            cbar = updateColorBar(colorBar, gci, label=label, **subkwargs)

        if markers:
            ticks = np.arange(len(uniquemarkers))
            cbar.set_ticks(ticks)
            labels = []
            for marker in uniquemarkers:
                labels.append(str((marker)))
            cbar.set_ticklabels(labels)

    if coverage is not None:
        if len(data) == mesh.cellCount():
            addCoverageAlpha(gci, coverage)
        else:
            raise BaseException('toImplement')
            # addCoverageAlpha(gci, pg.cellDataToPointData(mesh, coverage))

    if not hold or block is not False:
        if data is not None:
            if len(data) == mesh.cellCount():
                cb = CellBrowser(mesh, data, ax=ax)
                cb.connect()

        plt.show(block=block)
        try:
            plt.pause(0.01)
        except BaseException as _:

            pass

    if showMesh:
        pg.show(mesh, ax=ax)

    if hold:
        pg.mplviewer.hold(val=lastHoldStatus)

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
            except BaseException as _:
                pass
        print('..done')

    return ax, cbar


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
