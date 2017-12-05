#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Full waveform seismics and utilities
"""
import time
import numpy as np

import pygimli as pg
import pygimli.solver


def ricker(t, f, t0=0.0):
    """
    Create Ricker wavelet.

    Create a Ricker wavelet with a desired frequency and signal length.

    Parameters
    ----------
    t : array [float]
        Time base definition

    f : float
        Frequency of the wavelet in Hz

    t0 : float
        Offset time. Use 1/f to move the wavelet to start nearly from zero.

    Returns
    -------
    y : array_like
        Signal

    Examples
    --------
    Create 100 Hz Wavelet inside 1000 Hz sampled signal of length 0.1s.

    >>> from pygimli.physics.seismics import ricker
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> sampleFrequenz = 1000 # Hz
    >>> t = np.arange(0, 0.1, 1./sampleFrequenz)
    >>> r = ricker(100., t, 1./100)
    >>> lines = plt.plot(t,r,'-x')
    >>> plt.show()

    """
    fact = (np.pi**2) * (f**2) * ((t - t0)**2)
    y = (1.0 - 2.0 * fact) * np.exp(-fact)
    return y


def drawWiggle(ax, x, t, xoffset=0.0,
               posColor='red', negColor='blue', alpha=0.5, **kwargs):
    """
    Draw signal in wiggle style into a given ax.

    Parameters
    ----------
    ax : matplotlib ax
        To plot into

    x : array [float]
        Signal.

    t : array
        Time base for x

    xoffset : float
        Move wiggle plot along x axis

    posColor : str
        Need to be convertible to matplotlib color. Fill positive areas with.

    negColor : str
        Need to be convertible to matplotlib color. Fill negative areas with.

    alpha : float
        Opacity for fill area.

    **kwargs : dict()
        Will be forwarded to matplotlib.axes.fill

    Examples
    --------
    >>> from pygimli.physics.seismics import ricker, drawWiggle
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> t = np.arange(0, 0.02, 1./5000)
    >>> r = ricker(t, 100., 1./100)
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(1,1,1)
    >>> drawWiggle(ax, r, t, xoffset=0, posColor='red', negColor='blue',
    ...            alpha=0.2)
    >>> drawWiggle(ax, r, t, xoffset=1)
    >>> drawWiggle(ax, r, t, xoffset=2, posColor='black', negColor='white',
    ...            alpha=1.0)
    >>> ax.invert_yaxis()
    >>> plt.show()
    """
    wiggle, = ax.plot(x + xoffset, t, color='black', **kwargs)

    if len(t) > 1:
        tracefill = np.array(x)
        tracefill[0] = 0.0
        tracefill[-1] = 0.0
        tracefill[np.nonzero(x > x[0])] = 0
        fill, = ax.fill(tracefill + xoffset, t, color=negColor,
                        alpha=alpha, linewidth=0)

        tracefill = np.array(x)
        tracefill[0] = 0.0
        tracefill[-1] = 0.0
        tracefill[np.nonzero(x < x[0])] = 0
        fill, = ax.fill(tracefill + xoffset, t, color=posColor,
                        alpha=alpha, linewidth=0)


def drawSeismogramm(ax, mesh, u, dt, ids=None, pos=None, i=None):
    r"""Extract and show time series from wave field

    Parameters
    ----------

    ids: list
        List of node ids for the given mesh.
    pos : list
        List of positions for the given mesh. We will look for the nearest node.

    """
    ax.set_xlim(mesh.xmin(), mesh.xmax())
    ax.set_ylim(0., dt*len(u)*1000)
    ax.set_aspect(1)
    ax.set_ylabel('Time (ms)')
    ax.set_xlabel('Distance (m)')

    if i is None:
        i = len(u)-1

    t = np.linspace(0, i*dt*1000, i+1)

    if ids is None and pos is not None:
        ids = []
        for p in pos:
           ids.append(mesh.findNearestNode(p))

    xDist = mesh.node(0).pos().distance(mesh.node(1).pos())
    for iw, n in enumerate(ids):
        pos = mesh.node(n).pos()
        ax.plot(pos[0], 0.05, '^', color='black')
        print(pos)
        trace = pg.cat(pg.RVector(0), u[:(i+1), n])
#        print(i+1, n)
#        print(trace, (max(pg.abs(trace))))

#        if max(pg.abs(trace)) > 1e-8:

        trace *= np.exp(1/1000 * t)
        trace /= (max(pg.abs(trace)))
        trace *= 10

        drawWiggle(ax, trace, t=t, xoffset=pos[0])
    ax.invert_yaxis()


def solvePressureWave(mesh, velocities, times, sourcePos, uSource, verbose):
    r"""
    Solve pressure wave equation.

    Solve pressure wave for a given source function

    .. math::
        \frac{\partial^2 u}{\partial t^2} & = \diverg(a\grad u) + f\\
        finalize equation


    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        Mesh to solve on

    velocities : array
        velocities for each cell of the mesh

    time : array
        Time base definition

    sourcePos : RVector3
        Source position

    uSource : array
        u(t, sourcePos) source movement of length(times)
        Usually a Ricker wavelet of the desired seismic signal frequency.

    Returns
    -------
    u : RMatrix

        Return

    Examples
    --------
    See TODO write example
    """
    A = pg.RSparseMatrix()
    M = pg.RSparseMatrix()

#    F = pg.RVector(mesh.nodeCount(), 0.0)
    rhs = pg.RVector(mesh.nodeCount(), 0.0)
    u = pg.RMatrix(len(times), mesh.nodeCount())
    v = pg.RMatrix(len(times), mesh.nodeCount())

    sourceID = mesh.findNearestNode(sourcePos)

    if len(uSource) != len(times):
        raise Exception("length of uSource does not fit length of times: " +
                        str(uSource) + " != " + len(times))

    A.fillStiffnessMatrix(mesh, velocities * velocities)
    M.fillMassMatrix(mesh)
#    M.fillMassMatrix(mesh, velocities)

    FV = 0
    if FV:
        A, rhs = pygimli.solver.diffusionConvectionKernel(mesh,
                                                          velocities *
                                                          velocities,
                                                          sparse=1)

        M = pygimli.solver.identity(len(rhs))

        u = pg.RMatrix(len(times), mesh.cellCount())
        v = pg.RMatrix(len(times), mesh.cellCount())
        sourceID = mesh.findCell(sourcePos).id()

    dt = times[1] - times[0]

    theta = 0.51
    #theta = 1.
    S1 = M + dt * dt * theta * theta * A
    S2 = M

    solver1 = pg.LinSolver(S1, verbose=False)
    solver2 = pg.LinSolver(S2, verbose=False)
    swatch = pg.Stopwatch(True)

#    ut = pg.RVector(mesh.nodeCount(), .0)
#    vt = pg.RVector(mesh.nodeCount(), .0)

    timeIter1 = np.zeros(len(times))
    timeIter2 = np.zeros(len(times))
    timeIter3 = np.zeros(len(times))
    timeIter4 = np.zeros(len(times))

    for n in range(1, len(times)):
        u[n - 1, sourceID] = uSource[n - 1]

        # solve for u
        tic = time.time()
        # + * dt*dt * F
        rhs = dt * M * v[n-1] + (M - dt*dt * theta * (1. - theta) * A) * u[n-1]
        timeIter1[n - 1] = time.time() - tic

        tic = time.time()
        u[n] = solver1.solve(rhs)
        timeIter2[n - 1] = time.time() - tic

        # solve for v
        tic = time.time()
        rhs = M * v[n - 1] - dt * \
            ((1. - theta) * A * u[n - 1] + theta * A * u[n])  # + dt * F
        timeIter3[n - 1] = time.time() - tic

        tic = time.time()
        v[n] = solver2.solve(rhs)
        timeIter4[n - 1] = time.time() - tic

#         same as above
#        rhs = M * v[n-1] - dt * A * u[n-1] + dt * F
#        v[n] = solver1.solve(rhs)

        t1 = swatch.duration(True)

        if verbose and (n % verbose == 0):
            print(str(n) + "/" + str(len(times)),
                  times[n],
                  dt,
                  t1,
                  min(u[n]),
                  max(u[n]))

#    plt.figure()
#    plt.plot(timeIter1, label='Ass:1')
#    plt.plot(timeIter2, label='Sol:1')
#    plt.plot(timeIter3, label='Ass:2')
#    plt.plot(timeIter4, label='Sol:2')
#    plt.legend()
#    plt.figure()
    return u

if __name__ == "__main__":
    # from pygimli.physics.seismics import drawWiggle  # alternatively ricker
    import matplotlib.pyplot as plt

    t = np.arange(0, 0.02, 1. / 5000)
    r = ricker(t, 100., 1. / 100)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    drawWiggle(ax, r, t, xoffset=0, posColor='red', negColor='blue', alpha=0.2)
    drawWiggle(ax, r, t, xoffset=1)
    drawWiggle(ax, r, t, xoffset=2, posColor='black', negColor='white',
               alpha=1.0)
    plt.show()
