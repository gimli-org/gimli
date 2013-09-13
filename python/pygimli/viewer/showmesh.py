# -*- coding: utf-8 -*-

try:
    import pygimli as g
    from pygimli.mplviewer import drawMesh, drawModel, drawField, createColorbar
except ImportError:
    import sys
    sys.stderr.write('''ERROR: cannot import the library 'pygimli'. Ensure that pygimli is in your PYTHONPATH ''')
    sys.exit(1)

import matplotlib.pyplot as plt
import numpy as np

def showMesh(mesh, data=None, showLater=False, colorBar=False, axis=None,
             *args, **kwargs):
    """
    Syntactic sugar, short-cut to create axes and plot node or cell values
    return axes, cbar
        
    Parameters
    ----------
    """

    ret = []

    a = axis
    if a == None:
        fig = plt.figure()
        a = fig.add_subplot(1,1,1)

    gci = None
    cbar = None
    validData = False

    if data is None:
        drawMesh(a, mesh)
    else:
        if min(data) == max(data):
            print("No valid data",  min(data), max(data))
            drawMesh(a, mesh)
        else:
            validData = True
            if len(data) == mesh.cellCount():
                gci = drawModel(a, mesh, data, *args, **kwargs)
            elif len(data) == mesh.nodeCount():
                gci = drawField(a, mesh, data, *args, **kwargs)

    a.set_aspect('equal')

    if colorBar and validData:
        cbar = createColorbar(gci, *args, **kwargs)

    if not showLater:
        plt.show()

    #fig.show()
    #fig.canvas.draw()
    return a, cbar
#def showMesh(...)

