# -*- coding: utf-8 -*-
"""Plot 3D mesh."""

import numpy as np
import os
import sys

import pygimli as pg

vista = pg.optImport('vista', requiredFor="proper visualization in 3D")
if vista is None:
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    callback = 'showMesh3DFallback'
else:
    callback = 'showMesh3DVista'


PyQt5 = pg.optImport('PyQt5', requiredFor="pyGIMLi 3D viewer")
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

    if mesh.nodeCount() < 1e4:
        x = pg.x(mesh.positions())
        y = pg.y(mesh.positions())
        z = pg.z(mesh.positions())
        ax.scatter(x, y, z, 'ko')
    ax.set_title('Fallback, install vista for proper 3D visualization')

    return ax


def showMesh3DVista(mesh, data=None, **kwargs):
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

    # open with vista
    grid = vista.read(tmp)

    # add saved data from within the pg.mesh itself
    for k, v in mesh.dataMap():
        grid.cell_arrays[k] = np.asarray(v)
    # add given data from argument
    if data is not None:
        grid.cell_arrays['data'] = np.asarray(data)
        cMap = kwargs.pop('cMap', 'viridis')
        opacity = 1
    else:
        opacity = 0
        cMap = None

    params = grid.cell_arrays

    notebook = kwargs.pop('notebook', False)
    if notebook:
        vista.set_plot_theme('document')

    if use_gui and notebook is False:
        # app = Qt.QApplication()
        app = Qt.QApplication(sys.argv)
        s3d = Show3D(tmp)
        s3d.addMesh(grid, cMap=cMap, show_edges=True, opacity=opacity)
        app.exec_()

    # elif notebook is True:
    #     tool = vista.OrthogonalSlicer(grid)
    #     # Get the plotter for adding more datasets:
    #     p = tool.plotter
    #     p.show()

    else:
        plotter = vista.Plotter(notebook=notebook)
        # add the x, y, z arrows
        plotter.show_bounds()
        plotter.add_axes()

        plotter.add_mesh(grid, cmap=cMap, show_edges=True, opacity=0.1, color='grey')
        plotter.plot()
