import time
import os
import numpy as np
import pygimli as pg


def createResultFolder(subfolder):
    now = time.localtime()
    results = str(now.tm_year) + str(now.tm_mon).zfill(2) + \
        str(now.tm_mday).zfill(2) + '-' + str(now.tm_hour).zfill(2) + '.' + \
        str(now.tm_min).zfill(2)

    return createfolders(['./', results, subfolder])


def createfolders(foldername_list):
    """

    """
    path = ''

    for s in foldername_list:
        path = path + s + '/'

    try:
        os.makedirs(path)
    except OSError as e:
        if os.path.exists(path):
            print('Path "{}" already exists.'.format(path))
        else:
            print('Unable to create path "{}".'.format(path))
            raise(e)

    return path


def getSavePath(folder=None, subfolder=''):
    if folder is None:
        path = createResultFolder(subfolder)
    else:
        path = createfolders([folder, subfolder])
    return path


def createGradientModel2D(data, mesh, VTop, VBot):
    """
    Create 2D velocity gradient model.

    Creates a smooth, linear, starting model that takes the slope
    of the topography into account. This is done by fitting a straight line
    and using the distance to that as the depth value.
    Known as "The Marcus method"


    Parameters
    ----------
    data : pygimli DataContainer
        The topography list is in here.
    mesh : pygimli.Mesh
        The parametric mesh used for the inversion
    VTop : float
        The velocity at the surface of the mesh
    VBot : float
        The velocity at the bottom of the mesh

    Returns
    -------
    model : pygimli Vector, length M
        A numpy array with slowness values that can be used to start
        the inversion.
    """

    p = np.polyfit(pg.x(data.sensorPositions()), pg.y(data.sensorPositions()),
                   deg=1)  # slope-intercept form
    n = np.asarray([-p[0], 1.0])  # normal vector
    nLen = np.sqrt(np.dot(n, n))

    x = pg.x(mesh.cellCenters())
    z = pg.y(mesh.cellCenters())
    pos = np.column_stack((x, z))
    d = np.array([np.abs(np.dot(pos[i, :], n) - p[1]) / nLen
                  for i in range(pos.shape[0])])

    return np.interp(d, [min(d), max(d)], [1.0 / VTop, 1.0 / VBot])
