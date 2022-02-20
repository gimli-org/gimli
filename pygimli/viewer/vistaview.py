# -*- coding: utf-8 -*-
"""Plot 3D mesh."""

import sys
import matplotlib.pyplot as plt
import pygimli as pg

PyQt5 = pg.optImport('PyQt5', requiredFor="use pyGIMLi's 3D viewer")
pyvista = pg.optImport('pyvista', requiredFor="properly visualize 3D data")
panel = pg.optImport('panel', requiredFor='pyvista jupyter backend')

print(pyvista.__version__)
if pyvista is None:
    view3Dcallback = 'showMesh3DFallback'
else:
    view3Dcallback = 'showMesh3DVista'
    vers_users = pyvista.__version__
    vers_userf = float(pyvista.__version__[::-1].replace('.', '', 1)[::-1])
    vers_needs = '0.33.0'
    vers_needf = 0.330
    if vers_userf < vers_needf:
        pg.warn("Please consider updating PyVista to at least {}".format(
            vers_needs))
    from pygimli.viewer.pv import drawModel

if panel is None:
    pg.warn("Please install panel to plot 3D models in Jupyter Notebooks - conda install panel")
else:
    pyvista.set_jupyter_backend('panel')

#if pythreejs is None:
 #   pg.warn("Please install pythreejs")
#else:
 #   pyvista.set_jupyter_backend('pythreejs')

#if ipygany is None:
 #   pg.warn("Please install ipygany")
#else:
 #   pyvista.set_jupyter_backend('ipygany')    
# True for Jupyter notebooks and sphinx-builds
#_backend = plt.get_backend().lower()
#inline = "inline" in _backend or _backend == "agg"
#jupyter_backend='ipygany'
if PyQt5 is None:
    inline = True
else:
    from .pv.show3d import Show3D
    from PyQt5 import Qt
    inline = False


def showMesh3D(mesh, data, **kwargs):
    """
    Calling the defined function to show the 3D object.
    """
    if pg.rc['view3D'] == 'fallback':
        return showMesh3DFallback(mesh, data, **kwargs)

    return globals()[view3Dcallback](mesh, data, **kwargs)


def showMesh3DFallback(mesh, data, **kwargs):
    """
    Plot the 3D object sketchy.
    """
    ax = kwargs.pop('ax', None)

    from mpl_toolkits.mplot3d import Axes3D

    if ax is None or not isinstance(ax, Axes3D):
        fig = plt.figure()
        ax = fig.gca(projection='3d', proj_type='persp')
        #ax = fig.gca(projection='3d', proj_type='ortho')

    if mesh.boundaryCount() > 0:
        x, y, tri, z, dataIndex = pg.viewer.mpl.createTriangles(mesh)
        ax.plot_trisurf(x, y, tri, z)
    else:
        if mesh.nodeCount() < 1e4:
            x = pg.x(mesh.positions())
            y = pg.y(mesh.positions())
            z = pg.z(mesh.positions())
            ax.scatter(x, y, z, 'ko')
    ax.set_title('Fallback, install pyvista for proper 3D visualization')

    return ax, None


def showMesh3DVista(mesh, data=None, **kwargs):
    """
    Make use of the actual 3D visualization tool kit

    Parameters
    ----------
    data: pg.Vector or np.ndarray
        Dictionary of cell values, sorted by key. The values need to be
        numpy arrays.

    Returns
    -------
    plotter: pyvista.Plotter
        The plotter from pyvista.
    gui: Show3D [None]
        The small gui based on pyvista. Note that this is returned as 'None'
        if gui is passed as 'False'.

    Note
    ----
    Not having PyQt5 installed results in displaying the first key
    (and values) from the dictionary.
    """
    hold = kwargs.pop('hold', False)
    cmap = kwargs.pop('cmap', 'viridis')
    notebook = kwargs.pop('notebook', inline)
    gui = kwargs.pop('gui', not notebook)

    # add given data from argument
    if gui:
        app = Qt.QApplication(sys.argv)
        s3d = Show3D(app)
        s3d.addMesh(mesh, data, cmap=cmap, **kwargs)
        if not hold:
            s3d.wait()
        return s3d.plotter, s3d  # plotter, gui

    else:
        if notebook:
            plotter = drawModel(None, mesh, data, notebook=notebook, cmap=cmap,
                            **kwargs)
        if not hold:
            plotter.show()
        return plotter, None
