import os
from mayavi import mlab
from matplotlib import pyplot as plt

def showMesh3D(mesh, interactive=True):
    """
    Proof of concept for mayavi binding.

    Parameters
    ----------
    mesh : pygimli.Mesh
    interactive : bool
    """
    fig = mlab.figure(bgcolor=(1,1,1), size=(400,400))

    tmp = "/tmp/gimli_3d_view_%s.vtk" % os.getpid()

    mesh.exportVTK(tmp)
    src = mlab.pipeline.open(tmp, figure=fig)
    os.remove(tmp)

    if not interactive:
        mlab.options.offscreen = True

    surf = mlab.pipeline.surface(src, figure=fig, opacity=0.5)
    edges = mlab.pipeline.extract_edges(surf, figure=fig)
    mlab.pipeline.surface(edges, color=(0,0,0), figure=fig)
    if interactive:
        mlab.show()
    else:
        arr = mlab.screenshot(figure=fig, antialiased=True)
        plt.imshow(arr)
        plt.axis('off')
        plt.show()
