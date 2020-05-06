#!/usr/bin/env python
# encoding: utf-8
r"""
Geostatistical regularization
-----------------------------

In this example we illustrate the use of geostatistical constraints on
irregular meshes as presented by :cite:`jordi2018geostatistical`, compared to
classical smoothness operators of first or second kind.

In order to focus purely on the role of the constraints, we use a very simple
mapping forward operator that retrieves the values in the mesh at some given
positions. The constraints are therefore used as interpolation operators. Note
that the mapping forward operator can also be used for defining prior knowledge
if combined with another forward operator in a classical joint inversion
framework.

The elements of the covariance matrix :math:`\textbf{C}_{\text{M}}` are defined
by the distances H between the model cells i and j into the three directions

.. math::

    \textbf{C}_{\text{M},ij}=\sigma^{2}\exp{\left(
        -3\sqrt{\left(\frac{\textbf{H}_{x,ij}}{I_{x}}\right)^{2}+
                \left(\frac{\textbf{H}_{y,ij}}{I_{y}}\right)^{2}+
                \left(\frac{\textbf{H}_{z,ij}}{I_{z}}\right)^{2}}\right)}.

It defines the correlation between model cells as a function of correlation
lenghts (ranges) :math:`I_x`, :math:`I_y`, and :math:`I_z`. Of course, the
orientation of the coordinate axes is arbitrary and can be chosen by rotation.

According to inverse theory, we use the square root of the covariance matrix as
single-side regularization matrix C. It is computed by using an eigenvalue
decomposition based on the numpy linalg procedure

.. math::

    \textbf{C}_\text{M}^{-0.5} = \textbf{Q}\textbf{D}^{-0.5}\textbf{Q}^{T}

In order to avoid a matrix inverse (square root), a special matrix is derived
that does the decomposition and stores the eigenvectors and eigenvalues values.
A multiplication is done by multiplying with Q and scaling with the diagonal D.
This matrix is implemented in the class `CM05matrix` in the matrix module
:mod:`pygimli.matrix`.

It turns out that this matrix does not return a zero vector for a homogeneous
model. Therefore, an additional matrix called
`GeostatisticalConstraintsMatrix` is implemented where this spur is corrected
for. It is created by a mesh, a list of correlation lengths I, a dip angle that
distorts the x/y plane and a strike angle towards the third direction.
"""
# sphinx_gallery_thumbnail_number = 2

import pygimli as pg
import pygimli.meshtools as mt


"""
First we build a very simple forward operator, whose response returns the
values at given positions pos. In the initialization, the indices are stored
and a mapping matrix is created that projects the model vector to the forward
response. This matrix is also the Jacobian matrix for the inversion.
"""


class PriorFOP(pg.core.ModellingBase):
    """Forward operator for grabbing values."""

    def __init__(self, mesh, pos, verbose=False):
        """Init with mesh and some positions that are converted into ids."""
        super().__init__(self, verbose)
        self.setMesh(mesh)
        self.ind = [mesh.findCell(po).id() for po in pos]
        self.J = pg.core.SparseMapMatrix()
        self.J.resize(len(self.ind), mesh.cellCount())
        for i, ii in enumerate(self.ind):
            self.J.setVal(i, ii, 1.0)

        self.setJacobian(self.J)

    def response(self, model):
        """Return values at the indexed cells."""
        return model[self.ind]

    def createJacobian(self, model):
        """Do nothing (linear)."""
        pass


# %% create a simple mesh (e.g. a crosshole box geometry)
rect = mt.createRectangle(start=[0, -10], end=[10, 0])
mesh = mt.createMesh(rect, quality=34.5, area=0.1)
# choose some positions and create forward operator
pos = [[2, -2], [8, -2], [5, -5], [2, -8], [8, -8]]
fop = PriorFOP(mesh, pos)
# %% define some plotting options and a figure
kw = dict(
    colorBar=True,
    cMin=30,
    cMax=300,
    orientation='vertical',
    cMap='Spectral_r',
    logScale=True)
fig, ax = pg.plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
#
startModel = pg.Vector(mesh.cellCount(), 30)
fop.regionManager().setConstraintType(2)
tLog = pg.core.TransLog()
vals = [30, 50, 300, 100, 200]
# vals = [10, 20, 50, 30, 40]
inv = pg.core.Inversion(vals, fop, tLog, tLog)
inv.setRelativeError(0.05)  # 5 % error
inv.setModel(startModel)
inv.setLambda(200)
res = inv.run()
print(('{:.1f} ' * 5).format(*fop(res)), inv.chi2())
pg.show(mesh, res, ax=ax[0, 1], **kw)
#
fop.regionManager().setConstraintType(1)
res = inv.run()
print(('{:.1f} ' * 5).format(*fop(res)), inv.chi2())
pg.show(mesh, res, ax=ax[0, 0], **kw)
# geostatistic isotropic operator with 5m correlation length
C = pg.matrix.GeostatisticConstraintsMatrix(mesh=mesh, I=5)
fop.setConstraints(C)
inv.setModel(startModel)
inv.setLambda(30)
res = inv.run()
print(('{:.1f} ' * 5).format(*fop(res)), inv.chi2())
pg.show(mesh, res, ax=ax[1, 0], **kw)
ax[0, 0].set_title("1st order")
ax[0, 1].set_title("2nd order")
ax[1, 0].set_title("I=5")
# anisotropic operator with three times higher range and axis dipping 25deg
C = pg.matrix.GeostatisticConstraintsMatrix(mesh=mesh, I=[10, 3], dip=-25)
fop.setConstraints(C)
inv.setLambda(20)
inv.setModel(startModel)
res = inv.run()
print(('{:.1f} ' * 5).format(*fop(res)), inv.chi2())
pg.show(mesh, res, ax=ax[1, 1], **kw)
ax[1, 1].set_title("I=[10/3], dip=25")
# plot the position of the priors
for ai in ax.flat:
    for po in pos:
        ai.plot(*po, marker='o', markersize=10, color='k', fillstyle='none')

pg.plt.show()
# %%  We want to have a closer look at the constraint matrix
vec = pg.Vector(mesh.cellCount())
vec[fop.ind[2]] = 1.0
col = pg.log10(pg.abs(C * vec))
pg.show(mesh, col, cMin=-6, cMax=0, cMap="magma_r")
# it has a rather small footprint
# %% On the other hand the correlation/covariance matrix is
CM = pg.matrix.covarianceMatrix(mesh, I=[10, 2], dip=-25)
col = pg.log10(CM.dot(vec))
pg.show(mesh, col, cMap="magma_r")
# demonstrates the ellipsoidal footprint covering the whole mesh.
# %%
# For generating geostatistical media, one can use the function
# generateGeostatisticalModel. It computes a correlation matrix and multiplies
# it with a pseudo-random (randn) series.
model = pg.utils.generateGeostatisticalModel(mesh, I=[20, 4])
pg.show(mesh, model)
