# -*- coding: utf-8 -*-
"""Plot 3D mesh."""

import numpy as np
import os
import sys

import pygimli as pg

vtki = pg.optImport('vtki', requiredFor="Proper visualization in 3D")
if vtki is None:
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    callback = 'showMesh3DFallback'
else:
    callback = 'showMesh3DVTKI'


PyQt5 = pg.optImport('PyQt5', requiredFor="Make use of pyGIMLi 3D viewer")
if PyQt5 is None:
    use_gui = False
else:
    from .view3d import Show3D
    from PyQt5 import Qt
    use_gui = True


def showMesh3D(mesh, data, **kwargs):
    """
    Calling the defined function to show the 3D object.
    """
    return globals()[callback](mesh, data, **kwargs)


def showMesh3DFallback(mesh, data, **kwargs):
    """
    Plot the 3D object sketchy.
    """
    fig = plt.figure()
    ax = Axes3D(fig)

    if len(mesh.positions()) < 1e4:
        for pos in mesh.positions():
            ax.scatter(pos[0], pos[1], pos[2], 'ko')
            ax.set_title('Fallback')

    plt.show()


def showMesh3DVTKI(mesh, data=None, **kwargs):
    """
    Make use of the actual 3D visualization tool kit

    Parameter
    ---------
    data: pg.Vector or np.ndarray
        Dictionary of cell values, sorted by key. The values need to be
        numpy arrays.

    Note
    ----
    Not having PyQt5 installed results in displaying the first key
    (and values) from the dictionary.
    """
    # FIXME: maybe
    # temporary VTK write & read, may be replaced with direct VTK object.
    tmp = "/tmp/gimli_3d_view_%s.vtk" % os.getpid()
    mesh.exportVTK(tmp)

    # open with vtki
    grid = vtki.read(tmp)

    # add saved data from within the pg.mesh itself
    for k, v in mesh.dataMap():
        grid.cell_arrays[k] = np.asarray(v)
    # add given data from argument
    if data is not None:
        grid.cell_arrays['data'] = np.asarray(data)
    
    params = grid.cell_arrays

    cMap = kwargs.pop('cMap', 'viridis')
    notebook = kwargs.pop('notebook', False)
    if use_gui and notebook is False:
        # app = Qt.QApplication()
        app = Qt.QApplication(sys.argv)
        s3d = Show3D()
        s3d.addMesh(grid, cMap=cMap)
        # FIXME: using the qt interface somehow seems to delete the content of
        # cell-arrays
        if data is not None:
            s3d.addDataToMesh(params)
        app.exec_()

    # elif notebook is True:
    #     tool = vtki.OrthogonalSlicer(grid)
    #     # Get the plotter for adding more datasets:
    #     p = tool.plotter
    #     p.show()

    else:
        plotter = vtki.Plotter(notebook=notebook)
        # add the x, y, z arrows
        plotter.add_bounds_axes()
        plotter.add_axes()

        plotter.add_mesh(grid, cmap=cmap)
        plotter.plot()

