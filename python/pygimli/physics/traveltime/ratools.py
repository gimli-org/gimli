#!/usr/bin/env python
# -*- coding: utf-8 -*-
""""WRITEME"""

import numpy as np
import pygimli as pg


def createRAData(sensors):
    """
    Create a refraction data container.

    Default data container for shot and geophon at every sensor position.
    Predefined sensor indices's 's' as shot position and 'g' as
    geophon position.

    Parameters
    ----------
    sensors : ndarray | R3Vector
        Geophon and shot positions (same)

    Returns
    -------
    data : DataContainer
        Data container with predefined sensor indieces 's' and 'g' for

    """

    data = pg.DataContainer()
    data.registerSensorIndex('s')
    data.registerSensorIndex('g')

    data.setSensorPositions(sensors)
    S, G = [], []
    for s in range(data.sensorCount()):
        for g in range(data.sensorCount()):
            if s is not g:
                S.append(s)
                G.append(g)

    data.resize(len(S))
    data.set('s', S)
    data.set('g', G)
    data.set('valid', np.abs(np.sign(data('g') - data('s'))))

    return data


def createGradientModel2D(data, mesh, vTop, vBot):
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
    vTop : float
        The velocity at the surface of the mesh
    vBot : float
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

    return np.interp(d, [min(d), max(d)], [1.0 / vTop, 1.0 / vBot])
