# -*- coding: utf-8 -*-
"""Plot 3D mesh."""

import sys
import matplotlib.pyplot as plt
import pygimli as pg


PyQt5 = pg.optImport('PyQt5', requiredFor="use pyGIMLi's 3D viewer")
pyvista = pg.optImport('pyvista', requiredFor="properly visualize 3D data")

if pyvista is None:
    view3Dcallback = 'showMesh3DFallback'
    pg.rc['view3D'] = 'fallback'
else:
    view3Dcallback = 'showMesh3DVista'
    vers_users = pyvista.__version__
    vers_userf = float(pyvista.__version__[::-1].replace('.', '', 1)[::-1])
    vers_needs = '0.23.2'
    vers_needf = 0.232
    if vers_userf < vers_needf:
        pg.warn("Please consider updating PyVista to at least {}".format(
            vers_needs))
    pg.debug("Using pyvista: {}".format(vers_users))
    from pg.viewer.pv import drawModel

# True for Jupyter notebooks and sphinx-builds
_inlineBackend_ = not pg.viewer.isInteractive()

if PyQt5 is None or _inlineBackend_:
    _inlineBackend_ = True
else:
    from .pv.show3d import Show3D
    from PyQt5 import Qt
    _inlineBackend_ = False


def showMesh3D(mesh, data, **kwargs):
    """
    Calling the defined function to show the 3D object.
    """
    if pg.rc['view3D'] == 'fallback' or pg.rc['view3D'] == None:

        kwargs.pop('gui', None)
        kwargs['alpha'] = kwargs.pop('opacity', None)
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
        ax = fig.add_subplot(111, projection='3d', 
                             proj_type='persp'
                             #proj_type='ortho'
                             )

    if mesh.boundaryCount() > 0:
        x, y, tri, z, dataIndex = pg.viewer.mpl.createTriangles(mesh)
        kwargs.pop('notebook', False)
        kwargs.pop('hold', False)
        ax.plot_trisurf(x, y, tri, z, **kwargs)
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
    notebook = kwargs.pop('notebook', _inlineBackend_)
    gui = kwargs.pop('gui', False)

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
            pyvista.set_plot_theme('document')

        try:
            plotter = drawModel(None, mesh, data, notebook=notebook, cmap=cmap,
                                **kwargs)
        except Exception as e:
            print(e)
            pg.error("fix pyvista bindings")

        if not hold:
            plotter.show()
        return plotter, None
