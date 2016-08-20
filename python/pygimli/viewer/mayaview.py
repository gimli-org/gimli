# -*- coding: utf-8 -*-
"""Plot 3D mesh."""

import sys
import os

from matplotlib import pyplot as plt

showMesh3DFunct = 'showMesh3DMayvi'

try:
    from mayavi import mlab
except ImportError:
    error_msg = """Visualization in 3D requires Mayavi.\n""" + \
                """Try 'pip install mayavi' depending on your system.\n""" + \
                """Fallback to matplotlib \n"""
    sys.stderr.write(error_msg)
    showMesh3DFunct = 'showMesh3DFallback'


def showMesh3D(mesh, interactive=True):
    """TODO DOCUMENTME."""
    return globals()[showMesh3DFunct](mesh, interactive)


def showMesh3DFallback(mesh, interactive=True):
    """TODO DOCUMENTME."""
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = Axes3D(fig)

    if len(mesh.positions()) < 1e4:
        for pos in mesh.positions():
            ax.scatter(pos[0], pos[1], pos[2], 'ko')
            text = ("Proper visualization in 3D requires Mayavi.\n"
                    """Try 'pip install mayavi' depending on your system. """)
            ax.set_title(text + str(interactive))

    plt.show()


def showMesh3DMayvi(mesh, interactive=True):
    """Proof of concept for mayavi binding.

    Parameters
    ----------
    mesh : pygimli.Mesh
    interactive : bool
    """
    # should avoid opening of mayavi window when building documentation
    if not interactive:
        mlab.options.offscreen = True

    fig = mlab.figure(bgcolor=(1, 1, 1), size=(400, 400))

    # temporary VTK write & read, may be replaced with direct VTK object.
    tmp = "/tmp/gimli_3d_view_%s.vtk" % os.getpid()
    mesh.exportVTK(tmp)
    src = mlab.pipeline.open(tmp, figure=fig)
    os.remove(tmp)

    surf = mlab.pipeline.surface(src, figure=fig, opacity=0.5)
    edges = mlab.pipeline.extract_edges(surf, figure=fig)
    mlab.pipeline.surface(edges, color=(0, 0, 0), figure=fig)
    # mlab.pipeline.image_plane_widget(surf, colormap='gray',
    #                                  plane_orientation='x_axes')

    if interactive:
        mlab.show()
    else:
        arr = mlab.screenshot(figure=fig, antialiased=True)
        plt.imshow(arr)
        plt.axis('off')
        plt.show()
