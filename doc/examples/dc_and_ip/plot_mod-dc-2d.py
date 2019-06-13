#!/usr/bin/env python
# encoding: utf-8

r"""
Geoelectrics in 2.5D
--------------------

Geoelectrical modelling example in 2.5D. CR"""

###############################################################################
# Let us start with a mathematical formulation ...
#
# .. math::
#
#     \nabla\cdot(\sigma\nabla u)=-I\delta(\vec{r}-\vec{r}_{\text{s}}) \in R^3
#
# The source term is 3 dimensional but the distribution of the electrical
# conductivity :math:`\sigma(x,y)` should by 2 dimensional so we need a
# Fourier-Cosine-Transform from :math:`u(x,y,z) \mapsto u(x,y,k)` with the
# wave number :math:`k`. :math:`D^(a)(u(x,y,z)) \mapsto i^|a|k^a u(x,y)`
#
# .. math::
#     \nabla\cdot( \sigma \nabla u ) - \sigma k^2 u
#     &=-I\delta(\vec{r}-\vec{r}_{\text{s}}) \in R^2 \\
#     \frac{\partial }{\partial x} \left(\cdot( \sigma
#     \frac{\partial u}{\partial x}\right) + \frac{\partial }{\partial y} \left(\cdot(\sigma
#     \frac{\partial u}{\partial y}\right) - \sigma k^2 u & =
#     -I\delta(x-x_{\text{s}})\delta(y-y_{\text{s}}) \in R^2 \\
#     \frac{\partial u}{\partial \vec{n}} & = 0 \quad\mathrm{on}\quad\text{Surface} z=0
#

import numpy as np
import pygimli as pg

from pygimli.solver import solve

from pygimli.viewer import show
from pygimli.mplviewer import drawStreams


def uAnalytical(p, sourcePos, k):
    """Calculates the analytical solution for the 2.5D geoelectrical problem.
    
    Solves the 2.5D geoelectrical problem for one wave number k.
    It calculates the normalized (for injection current 1 A and sigma=1 S/m) 
    potential at position p for a current injection at position sourcePos.
    Injection at the subsurface is recognized via mirror sources along the
    surface at depth=0.
    
    Parameters
    ----------
    p : pg.Pos
        Position for the sought potential
    sourcePos : pg.Pos
        Current injection position.
    k : float
        Wave number 

    Returns
    -------
    u : float
        Solution u(p)
    """
    r1A = (p - sourcePos).abs()
    # Mirror on surface at depth=0
    r2A = (p - pg.RVector3(1.0, -1.0, 1.0) * sourcePos).abs()

    if r1A > 1e-12 and r2A > 1e-12:
        return (pg.math.besselK0(r1A * k) + pg.math.besselK0(r2A * k)) / (2.0 * np.pi)
    else:
        return 0.


def mixedBC(boundary, userData):
    """Mixed boundary conditions.

    Define the derivative of the analytical solution regarding the outer normal
    direction :math:`\vec{n}`. So we can define the value for the Neumann type
    Boundary conditions for the boundaries in the subsurface.
    """
    ### ignore surface boundaries for wildcard boundary condition
    if boundary.norm()[1] == 1.0:
        return 0

    sourcePos = userData['sourcePos']
    k = userData['k']
    r1 = boundary.center() - sourcePos
    
    # Mirror on surface at depth=0
    r2 = boundary.center() - pg.RVector3(1.0, -1.0, 1.0) * sourcePos
    r1A = r1.abs()
    r2A = r2.abs()

    n = boundary.norm()
    if r1A > 1e-12 and r2A > 1e-12:
        return k * ((r1.dot(n)) / r1A * pg.math.besselK1(r1A * k) +
                    (r2.dot(n)) / r2A * pg.math.besselK1(r2A * k)) / \
            (pg.math.besselK0(r1A * k) + pg.math.besselK0(r2A * k))
    else:
        return 0.


###############################################################################
#
#
def pointSource(mesh, source):
    """Define function for the current source term.

    :math:`\delta(x-pos), \int f(x) \delta(x-pos)=f(pos)=N(pos)`
    Right hand side entries will be shape functions(pos)
    """
    rhs = pg.Vector(mesh.nodeCount())
    
    cell = mesh.findCell(source)
    rhs.setVal(cell.N(cell.shape().rst(source)), cell.ids())
    return rhs

grid = pg.createGrid(x=np.linspace(-10.0, 10.0, 41),
                     y=np.linspace(-15.0,  0.0, 31))

grid = grid.createP2()

sourcePosA = [-5.0, -4.0]
sourcePosB = [+5.0, -4.0]

k = 1e-3
sigma = 1
u = solve(grid, a=sigma, b=-sigma * k*k, f=pointSource(grid, sourcePosA),
          bc={'Robin': [[1,2,4], mixedBC]},
          userData={'sourcePos': sourcePosA, 'k': k},
          verbose=True)

u -= solve(grid, a=sigma, b=-sigma * k*k, f=pointSource(grid, sourcePosB),
           bc={'Robin': ['*', mixedBC]},
           userData={'sourcePos': sourcePosB, 'k': k},
           verbose=True)

# uAna = pg.Vector(map(lambda p__: uAnalytical(p__, sourcePosA, k),
#                       grid.positions()))
# uAna -= pg.Vector(map(lambda p__: uAnalytical(p__, sourcePosB, k),
#                        grid.positions()))

# err = (1.0 -u/uAna) * 100.0

# print("error min max", min(err), max(err))

ax = show(grid, data=u, fillContour=True, colorBar=True, cMap="RdBu_r",
          orientation='horizontal', label='Solution u', nLevs=11,
          logScale=False, hold=True, showMesh=True)[0]

# Additional to the image of the potential we want to see the current flow too.
# The current flows along the gradient of our solution and can be plotted as
# stream lines. On default the drawStreams method draws one segment of a
# stream line per cell of the mesh. This can be a little confusing for dense
# meshes so we can give a second (coarse) mesh as a new cell basis to draw the
# streams. If the drawStreams get scalar data the gradients will be calculated.

gridCoarse = pg.createGrid(x=np.linspace(-10.0, 10.0, 20),
                           y=np.linspace(-15.0,   .0, 20))
drawStreams(ax, grid, u, coarseMesh=gridCoarse, color='Black')

pg.wait()
