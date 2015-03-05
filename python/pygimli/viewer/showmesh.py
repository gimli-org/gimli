# -*- coding: utf-8 -*-

"""
Generic mesh visualization tools. 
"""

try:
    import pygimli as pg
    from pygimli.mplviewer import drawMesh, drawModel, drawField
    from pygimli.mplviewer import drawSensors, showLater
    from pygimli.mplviewer import createColorbar, drawStreams, addCoverageAlpha

except ImportError:
    raise Exception('''ERROR: cannot import the library 'pygimli'.
        Ensure that pygimli is in your PYTHONPATH ''')

import matplotlib.pyplot as plt
import numpy as np


def show(mesh, *args, **kwargs):
    """
    Mesh and model visualization.
    
    Syntactic sugar to show a mesh with data.   
    Forwards to 
    :py:mod:`pygimli.viewer.showmesh.showMesh` or
    :py:mod:`pygimli.viewer.mayaview.showMesh3D` to show most of the
    possible 2D and 3D content.
    See tutorials and examples for usage hints.

    Parameters
    ----------

    mesh : :gimliapi:`GIMLI::Mesh`
        2D or 3D GIMLi mesh

    *args, **kwargs :
        Will be forwarded to the show functions.

    Returns
    -------

    Return the results from the show functions.
    """

    if isinstance(mesh, pg.Mesh):
        if mesh.dimension() == 2:
            return showMesh(mesh, *args, **kwargs)
        elif mesh.dimension() == 3:

            from .mayaview import showMesh3D

            return showMesh3D(mesh, *args, **kwargs)
        else:
            print("ERROR: Mesh not valid.")
    plt.pause(0.001)

def showMesh(mesh, data=None, showLater=False, colorBar=False, coverage=None,
             axes=None, **kwargs):
    """
    2D Mesh visualization.
    
    Create an axes and plot node or cell values for the given 2d mesh. 
    Returns the axes and the color bar.

    Parameters
    ----------

    mesh : :gimliapi:`GIMLI::Mesh`
        2D or 3D GIMLi mesh

    data : iterable [None]
        Optionally data to visualize.

        . None (draw mesh only)
            forward to :py:mod:`pygimli.mplviewer.meshview.drawMesh`

        . float per cell -- model, patch
            forward to :py:mod:`pygimli.mplviewer.meshview.drawModel`

        . float per node -- scalar field
            forward to :py:mod:`pygimli.mplviewer.meshview.drawField`

        . iterable of type [float, float] -- vector field
            forward to :py:mod:`pygimli.mplviewer.meshview.drawStreams`

        . pg.stdVectorRVector3 -- sensor positions
            forward to :py:mod:`pygimli.mplviewer.meshview.drawSensors`

    showLater : bool [false]
        Set interactive plot mode for matplotlib.
        If this is set to false [default] your script will stop to open
        a window with the figure and halted until you close this windows.
        You can set pg.showLater(1) to change this default behavior and
        pg.showNow() to force all pending figures to draw.

    colorBar : bool [false]
        Create and show a colorbar.

    coverage : iterable [None]
        Weight data by the given coverage array and fadeout the color.

    axes : matplotlib.Axes [None]
        Instead of create a new and empty axes, just draw into the a given.
        Useful to combine draws.

    *args, **kwargs :
        Will be forwarded to the draw functions and matplotlib methods,
        respectively.

    Returns
    -------
    axes : matplotlib.axes

    colobar : matplotlib.colobar
    """
    ax = axes

    if ax is None:
        fig, ax = plt.subplots()

    gci = None
    cbar = None
    validData = False
    
    
    if data is None:
        drawMesh(ax, mesh)
    elif type(data) == pg.stdVectorRVector3:
        drawSensors(ax, data)
    else:
        if hasattr(data[0], '__len__') and type(data) != np.ma.core.MaskedArray:
               
            if sum(data[:, 0]) != sum(data[:, 1]):
                drawStreams(ax, mesh, data, **kwargs)
            else:
                print("No valid stream data:",  data)
                drawMesh(ax, mesh)

        elif (min(data) == max(data)) or pg.haveInfNaN(data):
            
            print("No valid data: ",  min(data), max(data), pg.haveInfNaN(data))
            drawMesh(ax, mesh)
        else:
            validData = True
            if len(data) == mesh.cellCount():
                gci = drawModel(ax, mesh, data, **kwargs)
            elif len(data) == mesh.nodeCount():
                gci = drawField(ax, mesh, data, **kwargs)

    ax.set_aspect('equal')

    label=kwargs.pop('label', None)

    if colorBar and validData:
        cbar = createColorbar(gci, label=label, **kwargs)  # , *args, **kwargs) # causes problems!

    if coverage is not None:
        if len(data) == mesh.cellCount():
            addCoverageAlpha(gci, coverage)
        else:
            raise('toImplement')
            addCoverageAlpha(gci, pg.cellDataToPointData(mesh, coverage))

    if not showLater:
        plt.show()

    # fig.canvas.draw()
    return ax, cbar
# def showMesh(...)


def showBoundaryNorm(mesh, *args, **kwargs):
    """
    Show mesh boundaries normals.
    
    Show the mesh and draw a black line along the normal direction of all
    boundaries.

    Parameters
    ----------

    mesh : :gimliapi:`GIMLI::Mesh`
        2D or 3D GIMLi mesh

    Returns
    -------
    axes : matplotlib.axes
    """
    ax = kwargs.pop('axes', None)
    ax = show(mesh, showLater=True, axes=ax)[0]

    for b in mesh.boundaries():
        c1 = b.center()
        c2 = c1 + b.norm()
        ax.plot([c1[0], c2[0]],
                [c1[1], c2[1]], color='Black', *args, **kwargs)

    if not showLater:
        plt.show()

    return ax
