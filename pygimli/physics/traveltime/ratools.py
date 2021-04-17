#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tools for traveltime/refraction tomography."""

import numpy as np
import pygimli as pg


def createCrossholeData(sensors):
    """
    Create crosshole scheme assuming two boreholes with an equal number of sensors.

    Parameters
    ----------
    sensors : array (Nx2)
        Array with position of sensors.

    Returns
    -------
    scheme : DataContainer
        Data container with `sensors` predefined sensor indices 's' and 'g' for shot and receiver numbers.
    """
    from itertools import product

    if len(sensors) % 2 > 0:
        pg.error(
            "createCrossholeData is only defined for an equal number of sensors in two boreholes."
        )
    sensors = np.sort(sensors, axis=0)
    n = len(sensors) // 2
    numbers = np.arange(n)
    rays = list(product(numbers, numbers + n))

    # Empty container
    scheme = pg.DataContainer()

    # Add sensors
    for sen in sensors:
        scheme.createSensor(sen)

    # Add measurements
    rays = np.array(rays)
    scheme.resize(len(rays))
    scheme["s"] = rays[:, 0]
    scheme["g"] = rays[:, 1]
    scheme["valid"] = np.ones(len(rays))
    scheme.registerSensorIndex("s")
    scheme.registerSensorIndex("g")
    return scheme


def shotReceiverDistances(data, full=False):
    """Return vector of all distances (in m) between shot and receiver.
    for each 's' and 'g' in data.

    Parameters
    ----------
    data : pg.DataContainerERT

    full : bool [False]
        Get distances between shot and receiver position when full is True or
        only form x coordinate if full is False

    Returns
    -------
    dists :  array
        Array of distances

    """
    if full:
        pos = data.sensors()
        s, g = data.id("s"), data.id("g")
        off = [pos[s[i]].distance(pos[g[i]]) for i in range(data.size())]
        return np.absolute(off)
    else:
        px = pg.x(data)
        gx = np.array([px[g] for g in data.id("g")])
        sx = np.array([px[s] for s in data.id("s")])
        return np.absolute(gx - sx)


def createRAData(sensors, shotDistance=1):
    """Create a refraction data container.

    Default data container for shot and geophon at every sensor position.
    Predefined sensor indices's 's' as shot position and 'g' as
    geophon position.

    Parameters
    ----------
    sensors: ndarray | R3Vector
        Geophon and shot positions (same)
    shotDistances: int [1]
        Distance between shot indices.
        
    Returns
    -------
    data : DataContainer
        Data container with predefined sensor indices 's' and 'g' for shot and receiver numbers.
    """
    data = pg.DataContainer()
    data.registerSensorIndex("s")
    data.registerSensorIndex("g")

    if isinstance(sensors, np.ndarray):
        if len(sensors.shape) == 1:
            for x in sensors:
                data.createSensor([x, 0.0, 0.0])
        else:
            data.setSensorPositions(sensors)

    else:
        data.setSensorPositions(sensors)

    S, G = [], []
    for s in range(0, data.sensorCount(), shotDistance):
        for g in range(data.sensorCount()):
            if s is not g:
                S.append(s)
                G.append(g)

    data.resize(len(S))
    data.set("s", S)
    data.set("g", G)
    data.set("valid", np.abs(np.sign(data("g") - data("s"))))

    return data


def createGradientModel2D(data, mesh, vTop, vBot):
    """Create 2D velocity gradient model.

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
    p = np.polyfit(pg.x(data), pg.y(data), deg=1)  # slope-intercept form
    n = np.asarray([-p[0], 1.0])  # normal vector
    nLen = np.sqrt(np.dot(n, n))

    x = pg.x(mesh.cellCenters())
    z = pg.y(mesh.cellCenters())
    pos = np.column_stack((x, z))

    d = np.array(
        [np.abs(np.dot(pos[i, :], n) - p[1]) / nLen for i in range(pos.shape[0])]
    )

    return 1.0 / np.interp(d, [min(d), max(d)], [vTop, vBot])
