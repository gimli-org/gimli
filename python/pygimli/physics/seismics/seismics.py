#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Full waveform seismics and utilities
"""

import numpy as np

def ricker(f, t, t0=0.0):
    """
    Create Ricker wavelet.
        
    Create a Ricker wavelet with a desired frequency and signal length.
        
    Parameters
    ----------
    f : float
        Frequency of the wavelet in Hz
    
    t : array [float]
        Time base definition
    
    t0 : float
        Offset time. Use 1/f to move the wavelet to start nearly from zero.
        
    Returns
    -------
    y : array_like
        Signal
        
    Examples
    --------
    Create 100 Hz Wavelet inside 1000 Hz sampled signla of length 0.1s.
    
    >>> from pygimli.physics.seismics import ricker
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> sampleFrequenz = 1000 # Hz
    >>> t = np.arange(0, 0.1, 1./sampleFrequenz)
    >>> r = ricker(100., t, 1./100)
    >>> plt.plot(t,r,'-x')
    >>> plt.show()

    """
    fact = (np.pi**2) * (f**2) * ((t-t0)**2)
    y = (1.0 - 2.0 * fact) * np.exp(-fact)
    return y

def wiggle(axes, x, t, xoffset=0.0,
           posColor='red', negColor='blue', alpha=0.5, **kwargs):
    """
    Draw signal in wiggle style into a given axes.

    Parameters
    ----------
    axes : matplotlib axes
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
    >>> from pygimli.physics.seismics import ricker, wiggle
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> t = np.arange(0, 0.02, 1./5000)
    >>> r = ricker(100., t, 1./100)
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(1,1,1)
    >>> 
    >>> wiggle(ax, r, t, xoffset=0, posColor='red', negColor='blue', alpha=0.2)
    >>> wiggle(ax, r, t, xoffset=1)
    >>> wiggle(ax, r, t, xoffset=2, posColor='black', negColor='white', alpha=1.0)
    >>> plt.show()

    """
    wiggle, = axes.plot(x + xoffset, t, color='black', **kwargs)
    
    if len(t) > 1:
        tracefill = np.array(x)
        tracefill[0] = 0.0
        tracefill[-1] = 0.0
        tracefill[np.nonzero(x > x[0])] = 0
        fill, = axes.fill(tracefill + xoffset, t, color=negColor,
                          alpha=alpha, linewidth=0)
        
        tracefill = np.array(x)
        tracefill[0] = 0.0
        tracefill[-1] = 0.0
        tracefill[np.nonzero(x < x[0])] = 0
        fill, = axes.fill(tracefill + xoffset, t, color=posColor,
                          alpha=alpha, linewidth=0)

def solvePressureWave(mesh, velocities, times, sourcePos, uSource):
    r"""
    Solve pressure wave equation.
        
    Solve pressure wave for a given source function
        
    .. math::
        \frac{\partial^2 u}{\partial t^2} & = \diverg(a\grad u) \\
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

    F = pg.RVector(grid.nodeCount(), 0.0)
    rhs = pg.RVector(grid.nodeCount(), 0.0)
    u = pg.RMatrix(len(t), grid.nodeCount())
    v = pg.RMatrix(len(t), grid.nodeCount())

    sourceID = grid.findNearestNode(sourcePos)
    
    if len(uSource) != len(times):
        raise Exception("length of uSource does not fit length of times: " +\
                        str(uSource) + " != " + len(times))
    
    u[0, sourceID] = uSource[0]

    A.fillStiffnessMatrix(grid, velocities*velocities)
    M.fillMassMatrix(grid)
    #M.fillMassMatrix(grid, velocities)
    
    S1 = M + k*k * A
    S2 = M
    solver1 = pg.LinSolver(S1, verbose=False)
    solver2 = pg.LinSolver(S2, verbose=False)
    swatch=pg.Stopwatch(True)

    dt = times[1] - times[0]

    ut = pg.RVector(grid.nodeCount(), .0 )
    vt = pg.RVector(grid.nodeCount(), .0 )
    
    for n in range(1, len(times)):
        u[n-1, source] = uSource[n-1]

        # solve for u
        rhs = M * u[n-1] + dt * M * v[n-1] + dt*dt * F
        u[n] = solver1.solve(rhs)

        # solve for v
        rhs = M * v[n-1] - dt * A * u[n] + dt * F
        v[n] = solver2.solve(rhs)
    
        t1 = swatch.duration(True)
    
        if verbose and (n%verbose == 0):
            print(n, t[n], k, len(y), len(x), t1, t2, min(u[n]), max(u[n]))
    
    return u

    
if __name__ == "__main__":
    from pygimli.physics.seismics import ricker, wiggle
    import matplotlib.pyplot as plt
    import numpy as np
    t = np.arange(0, 0.02, 1./5000)
    r = ricker(100., t, 1./100)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
   
    wiggle(ax, r, t, xoffset=0, posColor='red', negColor='blue', alpha=0.2)
    wiggle(ax, r, t, xoffset=1)
    wiggle(ax, r, t, xoffset=2, posColor='black', negColor='white', alpha=1.0)
    plt.show()

