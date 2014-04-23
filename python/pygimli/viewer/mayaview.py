# -*- coding: utf-8 -*-
import sys
import os

from matplotlib import pyplot as plt

try:
    from mayavi import mlab
except ImportError:
    error_msg  = """Visualization in 3D requires Mayavi.\n""" + \
                 """Try 'pip install mayavi' depending on your system.\n"""
    sys.stderr.write(error_msg)
    raise Exception("Visualization in 3D requires Mayavi.")

def showMesh3D(mesh, interactive=True):
    """
    Proof of concept for mayavi binding.

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
    #mlab.pipeline.image_plane_widget(surf, colormap='gray',
                                     #plane_orientation='x_axes')

    if interactive:
        mlab.show()
    else:
        arr = mlab.screenshot(figure=fig, antialiased=True)
        plt.imshow(arr)
        plt.axis('off')
        plt.show()
