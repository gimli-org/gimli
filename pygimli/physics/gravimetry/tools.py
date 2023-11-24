"""Tools for magnetics."""
import numpy as np
import pygimli as pg

def depthWeighting(mesh, cell=False, z0=25, height=0, power=1.5, normalize=True):
    """Return Li&Oldenburg like depth weighting of boundaries or cells."""
    if cell:
        z = np.abs(pg.z(mesh.cellCenter()))
    else:
        z = np.abs(pg.z(mesh.innerBoundaryCenters()))

    weight = 1 / ((z+height)/z0 + 1)**power
    if normalize:
        weight /= np.median(weight)  # normalize that maximum is 1

    return weight
