# -*- coding: utf-8 -*-
"""Plot 3D mesh."""

import numpy as np
import os
import sys

import pygimli as pg

pyvista = pg.optImport('pyvista', requiredFor="proper visualization in 3D")
if pyvista is None:
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
    ax.set_title('Fallback, install pyvista for proper 3D visualization')

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

    # open with pyvista
    grid = pyvista.read(tmp)

    hold = kwargs.pop("hold", False)
    cMap = kwargs.pop('cMap', 'viridis')

    # add given data from argument
    add_args = {}
    if data is not None:
        label = kwargs.pop("label", "data")
        if len(data) == mesh.cellCount():
            grid.cell_arrays[label] = np.asarray(data)
        elif len(data) == mesh.nodeCount():
            grid.point_arrays[label] = np.asarray(data)
        grid.set_active_scalar(label)
        add_args["opacity"] = 1
    else:
        add_args["opacity"] = 0.1
        add_args["show_scalar_bar"] = False

    notebook = kwargs.pop('notebook', False)
    if notebook:
        pyvista.set_plot_theme('document')

    if use_gui and notebook is False:
        # add saved data from within the pg.mesh itself
        for k, v in mesh.dataMap():
            grid.cell_arrays[k] = np.asarray(v)
        # app = Qt.QApplication()
        app = Qt.QApplication(sys.argv)
        s3d = Show3D(tmp)
        s3d.addMesh(grid, cMap=cMap)
        app.exec_()

    # elif notebook is True:
    #     tool = pyvista.OrthogonalSlicer(grid)
    #     # Get the plotter for adding more datasets:
    #     p = tool.plotter
    #     p.show()

    else:
        plotter = pyvista.Plotter(notebook=notebook)
        plotter.show_bounds()
        plotter.add_axes()
        plotter.add_mesh(grid, cmap=cMap, show_edges=True, **add_args)
        if data is not None:
            plotter.mesh.set_active_scalar(label)
        if not hold:
            plotter.plot()
        return plotter
