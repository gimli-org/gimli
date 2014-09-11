# -*- coding: utf-8 -*-

try:
    import pygimli as g
    from pygimli.mplviewer import drawMesh, drawModel, drawField, createColorbar, drawStreamLines2
    
except ImportError:
    raise Exception('''ERROR: cannot import the library 'pygimli'. Ensure that pygimli is in your PYTHONPATH ''')

import matplotlib.pyplot as plt
import numpy as np

def show(mesh, *args, **kwargs):
    """Syntactic sugar."""
    if isinstance(mesh, g.Mesh):
        if mesh.dimension() == 2:
            return showMesh(mesh, *args, **kwargs)
        elif mesh.dimension() == 3:
            
            from .mayaview import showMesh3D

            return showMesh3D(mesh, **kwargs)
        else:
            print("ERROR: Mesh not valid.")


def showMesh(mesh, data=None, showLater=False, colorBar=False, axes=None,
             *args, **kwargs):
    """
    Syntactic sugar, short-cut to create axes and plot node or cell values
    return axes, cbar

    Parameters
    ----------
    """

    ret = []

    ax = axes
    
    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        
    gci = None
    cbar = None
    validData = False

    if data is None:
        drawMesh(ax, mesh)
    else:
        #print(data[0], type(data[0]))
        if hasattr(data[0], '__len__') and type(data) != np.ma.core.MaskedArray:
            if sum(data[:,0]) != sum(data[:,1]):
                drawStreamLines2(ax, mesh, data)
            else:
                print("No valid stream data:",  data)
                drawMesh(ax, mesh)

        elif min(data) == max(data):
            print(("No valid data",  min(data), max(data)))
            drawMesh(ax, mesh)
        else:
            validData = True
            if len(data) == mesh.cellCount():
                gci = drawModel(ax, mesh, data, *args, **kwargs)
            elif len(data) == mesh.nodeCount():
                gci = drawField(ax, mesh, data, *args, **kwargs)

    ax.set_aspect('equal')

    if colorBar and validData:
        cbar = createColorbar(gci, *args, **kwargs)

    if not showLater:
        plt.show()

    #fig.show()
    #fig.canvas.draw()
    return ax, cbar
#def showMesh(...)

def showBoundaryNorm(mesh, *args, **kwargs):
    """"""

    ax = showMesh(mesh, showLater=True)[0]

    for b in mesh.boundaries():
        c1 = b.center()
        c2 = c1 + b.norm()
        ax.plot([c1[0], c2[0]],
                [c1[1], c2[1]], color='Black')

    plt.show()

    return ax











