#!/usr/bin/env python
# encoding: utf-8

r"""
Geoelectrics in 2.5D
--------------------

Geoelectrical modelling example in 2.5D. CR

"""
###############################################################################
# Let us start with the governing mathematical formulation:
#
# .. math::
#
#     \nabla\cdot(\sigma\nabla u)=-I\delta(\vec{r}-\vec{r}_{\text{s}}) \in R^3
#
# The source term is 3 dimensional but the distribution of the electrical
# conductivity :math:`\sigma(x,y)` should be 2D so we apply a
# Fourier-Cosine-Transform from :math:`u(x,y,z) \mapsto u(x,k,z)` with the
# wave number :math:`k`. (:math:`D^{(a)}(u(x,y,z)) \mapsto i^{|a|}k^a u(x,z)`)
#
# .. math::
#     \nabla\cdot( \sigma \nabla u ) - \sigma k^2 u
#     &=-I\delta(\vec{r}-\vec{r}_{\text{s}}) \in R^2 \\
#     \frac{\partial }{\partial x}\left(\cdot \sigma \frac{\partial u}{\partial x}\right) +
#     \frac{\partial }{\partial z}\left(\cdot\sigma \frac{\partial u}{\partial z}\right) -
#     \sigma k^2 u & =
#     -I\delta(x-x_{\text{s}})\delta(z-z_{\text{s}}) \in R^2 \\
#     \frac{\partial u}{\partial \vec{n}} & = 0 \quad\mathrm{on}\quad\text{Surface}\quad z=0 \\
#     \frac{\partial u}{\partial \vec{n}} & = a u \quad\mathrm{on}\quad\text{Subsurface}\quad z<0
#

import matplotlib
import numpy as np
import pygimli as pg

from pygimli.viewer import show
from pygimli.viewer.mpl import drawStreams

###############################################################################
# We know the exact solution:
#
def uAnalytical(p, sourcePos, k, sigma=1):
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
    r2A = (p - pg.Pos([1.0, -1.0])*sourcePos).abs()

    if r1A > 1e-12 and r2A > 1e-12:
        return  1 / (2.0 * np.pi) * 1/sigma * \
            (pg.math.besselK0(r1A * k) + pg.math.besselK0(r2A * k))
    else:
        return 0.

###############################################################################
# We assume the so called mixed boundary conditions.
#
def mixedBC(boundary, userData):
    """Mixed boundary conditions.

    Define the derivative of the analytical solution regarding the outer normal
    direction :math:`\vec{n}`. So we can define the values for Robin type
    boundary conditions :math:`\frac{\partial u}{\partial \vec{n}} = -au` 
    for the boundaries on the subsurface.

    """
    ### ignore surface boundaries for wildcard boundary condition
    if boundary.norm()[1] == 1.0:
        return 0

    sourcePos = userData['sourcePos']
    k = userData['k']
    sigma = userData['s']
    r1 = boundary.center() - sourcePos

    # Mirror on surface at depth=0
    r2 = boundary.center() - pg.Pos(1.0, -1.0) * sourcePos
    r1A = r1.abs()
    r2A = r2.abs()

    n = boundary.norm()
    if r1A > 1e-12 and r2A > 1e-12:
        return sigma * k * ((r1.dot(n)) / r1A * pg.math.besselK1(r1A * k) +
                            (r2.dot(n)) / r2A * pg.math.besselK1(r2A * k)) / \
            (pg.math.besselK0(r1A * k) + pg.math.besselK0(r2A * k))
            
    else:
        return 0.

###############################################################################
# We assemble the right hand side (rhs) for the singular current term ourself 
# since this cannot be yet done efficiently by pg.solve.
#
def rhsPointSource(mesh, source):
    """Define function for the current source term.

    :math:`\delta(x-pos), \int f(x) \delta(x-pos)=f(pos)=N(pos)`
    Right hand side entries will be shape functions(pos)
    """
    rhs = pg.Vector(mesh.nodeCount())

    cell = mesh.findCell(source)
    rhs.setVal(cell.N(cell.shape().rst(source)), cell.ids())
    return rhs

###############################################################################
# No lets create a suitabel mesh and solve the equation with pg.solve
#
mesh = pg.createGrid(x=np.linspace(-10.0, 10.0, 41),
                     y=np.linspace(-15.0,  0.0, 31))
mesh = mesh.createP2()

sourcePosA = [-5.25, -3.75]
sourcePosB = [+5.25, -3.75]

k = 1e-2
sigma = 1.0
bc={'Robin': {'1,2,4': mixedBC}}
u = pg.solve(mesh, a=sigma, b=-sigma * k*k, 
             rhs=rhsPointSource(mesh, sourcePosA),
             bc=bc, userData={'sourcePos': sourcePosA, 'k': k, 's':sigma},
             verbose=True)

u -= pg.solve(mesh, a=sigma, b=-sigma * k*k, 
              rhs=rhsPointSource(mesh, sourcePosB), 
              bc=bc, userData={'sourcePos': sourcePosB, 'k': k, 's':sigma},
              verbose=True)

ax = show(mesh, data=u, cMap="RdBu_r", cMin=-1, cMax=1,
          orientation='horizontal', label='Potential $u$', 
          nCols=16, nLevs=9, logScale=False, showMesh=True)[0]

###############################################################################
# Additional to the image of the potential we want to see the current flow too.
# The current flows along the gradient of our solution and can be plotted as
# stream lines. On default the drawStreams method draws one segment of a
# stream line per cell of the mesh. This can be a little confusing for dense
# meshes so we can give a second (coarse) mesh as a new cell basis to draw the
# streams. If the drawStreams get scalar data the gradients will be calculated.

gridCoarse = pg.createGrid(x=np.linspace(-10.0, 10.0, 20),
                           y=np.linspace(-15.0,   .0, 20))
drawStreams(ax, mesh, u, coarseMesh=gridCoarse, color='Black')
###############################################################################
# We know the exact solution so we can compare the results. 
# Unfortunately, the point source singularity does not allow a good integration
# measure for the accuracy of the resulting field so we just look for the 
# differences.
#
uAna = pg.Vector(list(map(lambda p__: uAnalytical(p__, sourcePosA, k, sigma),
                      mesh.positions())))
uAna -= pg.Vector(list(map(lambda p__: uAnalytical(p__, sourcePosB, k, sigma),
                       mesh.positions())))

ax = show(mesh, data=pg.abs(uAna-u), cMap="Reds",
          orientation='horizontal', label='|$u_{exact}$ -$u$|', 
          logScale=True, cMin=1e-7, cMax=1e-1,
          contourLines=False,
          nCols=12, 
          nLevs=7, 
          showMesh=True)[0]

# print('l2:', pg.pf(pg.solver.normL2(uAna-u)))
# print('L2:', pg.pf(pg.solver.normL2(uAna-u, mesh)))
# print('H1:', pg.pf(pg.solver.normH1(uAna-u, mesh)))
np.testing.assert_approx_equal(pg.solver.normL2(uAna-u, mesh), 
                               0.02415, significant=3)

