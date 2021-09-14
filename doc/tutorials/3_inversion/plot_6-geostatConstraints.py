#!/usr/bin/env python
# encoding: utf-8
r"""
Geostatistical regularization
-----------------------------

In this example we illustrate the use of geostatistical constraints on
irregular meshes as presented by :cite:`jordi2018geostatistical`, compared to
classical smoothness operators of first or second kind.

The elements of the covariance matrix :math:`\textbf{C}_{\text{M}}` are defined
by the distances H between the model cells i and j into the three directions

.. math::

    \textbf{C}_{\text{M},ij}=\sigma^{2}\exp{\left(
        -3\sqrt{\left(\frac{\textbf{H}^x_{ij}}{I_{x}}\right)^{2}+
                \left(\frac{\textbf{H}^y_{ij}}{I_{y}}\right)^{2}+
                \left(\frac{\textbf{H}^z_{ij}}{I_{z}}\right)^{2}}\right)}.

It defines the correlation between model cells as a function of correlation
lenghts (ranges) :math:`I_x`, :math:`I_y`, and :math:`I_z`. Of course, the
orientation of the coordinate axes is arbitrary and can be chosen by rotation.
Let us illustrate this by a simple mesh:
"""

# %%
# Computing covariance and constraint matrices
# --------------------------------------------
# We create a simple mesh using a box geometry
import matplotlib.pyplot as plt
import pygimli as pg
import pygimli.meshtools as mt

# We create a rectangular domain and mesh it with small triangles
rect = mt.createRectangle(start=[0, -10], end=[10, 0])
mesh = mt.createMesh(rect, quality=34.5, area=0.1)

# %%
# We compute such a covariance matrix by calling
CM = pg.utils.covarianceMatrix(mesh, I=5)  # I taken for both x and y
# We search for the cell where the midpoint (5, -5) is located in
ind = mesh.findCell([5, -5]).id()
# and plot the according column using index access (numpy)
ax, cb = pg.show(mesh, CM[:, ind], cMap="magma_r")

# %%
# According to inverse theory, we use the square root of the covariance matrix
# as single-side regularization matrix C. It is computed by using an eigenvalue
# decomposition
#
# .. math::
#
#     \textbf{C}_\text{M} = \textbf{Q}\textbf{D}\textbf{Q}^{T}
#
# based on LAPACK (numpy.linalg). The inverse square root is defined by
#
# .. math::
#
#     \textbf{C}_\text{M}^{-0.5} = \textbf{Q}\textbf{D}^{-0.5}\textbf{Q}^{T}
#
# In order to avoid a matrix inverse (square root), a special matrix is derived
# doing the decomposition and storing the eigenvectors and eigenvalues values.
# A multiplication is done by multiplying with Q and scaling with the diagonal.
# This matrix is implemented in the :mod:`pygimli.matrix` module
# by the class :py:mod:`pg.matrix.Cm05Matrix`

Cm05 = pg.matrix.Cm05Matrix(CM)
# %%
# However, this matrix does not return a zero vector for a constant vector
out = Cm05 * pg.Vector(mesh.cellCount(), 1.0)
print("min/max value ", min(out), max(out))

# %%
# as desired for a roughness operator. Therefore, an additional matrix called
# :py:mod:`pg.matrix.GeostatisticalConstraintsMatrix`
# was implemented where this spur is corrected for.
# It is, like the correlation matrix, created by a mesh, a list of correlation
# lengths I, a dip angle# that distorts the x/y plane and a strike angle
# towards the third direction.
#
C = pg.matrix.GeostatisticConstraintsMatrix(mesh=mesh, I=5)

# %%
# In order to extract a certain column, we generate a vector with a single 1
vec = pg.Vector(mesh.cellCount())
vec[ind] = 1.0
ax, cb = pg.show(mesh, pg.log10(pg.abs(C*vec)),
                 cMin=-6, cMax=0, cMap="magma_r")

# %%
# The constraints have a rather small footprint compared to the correlation
# (note the logarithmic scale) but still to the whole mesh unlike the classical
# constraint matrices that only include relations to neighboring cells.

# %%
# Such a matrix can also be defined for different ranges and dip angles, e.g.
Cdip = pg.matrix.GeostatisticConstraintsMatrix(mesh=mesh, I=[9, 2], dip=-25)
ax, cb = pg.show(mesh, pg.log10(pg.abs(Cdip*vec)),
                 cMin=-6, cMax=0, cMap="magma_r")

# %%
# In order to illustrate the role of the constraints, we use a very simple
# mapping forward operator that retrieves the values in the mesh at some given
# positions. The constraints are therefore used as interpolation operators.
# Note that the mapping forward operator can also be used for defining prior
# knowledge if combined with another forward operator in a classical joint
# inversion framework.
# In the initialization, the indices are stored and a mapping matrix is created
# that projects the model vector to the forward response.
# This matrix is also the Jacobian matrix for the inversion.


class PriorFOP(pg.Modelling):
    """Forward operator for grabbing values."""

    def __init__(self, mesh, pos, **kwargs):
        """Init with mesh and some positions that are converted into ids."""
        super().__init__(**kwargs)
        self.setMesh(mesh)
        self.ind = [mesh.findCell(po).id() for po in pos]
        self.J = pg.SparseMapMatrix()
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


# %%
# Inversion with geostatistical constraints
# -----------------------------------------
# We choose some positions and initialize the forward operator
pos = [[2, -2], [8, -2], [5, -5], [2, -8], [8, -8]]
fop = PriorFOP(mesh, pos)
# For plotting the results, we create a figure and define some plotting options
fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
kw = dict(
    colorBar=True,
    cMin=30,
    cMax=300,
    orientation='vertical',
    cMap='Spectral_r',
    logScale=True)
# We want to use a homogenenous starting model
tLog = pg.trans.TransLog()
vals = [30, 50, 300, 100, 200]
# We assume a 5% relative accuracy of the values
error = pg.Vector(len(vals), 0.05)
# set up data and model transformation log-scaled
inv = pg.Inversion(fop=fop)
inv.transData = tLog
inv.transModel = tLog
inv.lam = 40
startModel = pg.Vector(mesh.cellCount(), 30)
inv.startModel = startModel
# Initially, we use the first-order constraints (default)
res = inv.run(vals, error, cType=1, lam=35)
print(('Ctype=1: ' + '{:.1f} ' * 6).format((*fop(res)), inv.chi2()))
pg.show(mesh, res, ax=ax[0, 0], **kw)
ax[0, 0].set_title("1st order")
# Next, we use the second order (curvature) constraint type
res = inv.run(vals, error, cType=2, lam=1000)
print(('Ctype=2: ' + '{:.1f} ' * 6).format((*fop(res)), inv.chi2()))
pg.show(mesh, res, ax=ax[0, 1], **kw)
ax[0, 1].set_title("2nd order")
# Now we set the geostatistic isotropic operator with 5m correlation length
fop.setConstraints(C)
res = inv.run(vals, error, lam=25)
print(('Cg-5/5m: ' + '{:.1f} ' * 6).format((*fop(res)), inv.chi2()))
pg.show(mesh, res, ax=ax[1, 0], **kw)
ax[1, 0].set_title("I=5")
# and finally we use the dipping constraint matrix
fop.setConstraints(Cdip)
res = inv.run(vals, error, lam=35)
print(('Cg-9/2m: ' + '{:.1f} ' * 6).format((*fop(res)), inv.chi2()))
pg.show(mesh, res, ax=ax[1, 1], **kw)
ax[1, 1].set_title("I=[10/2], dip=25")
# plot the position of the priors
for ai in ax.flat:
    for po in pos:
        ai.plot(*po, marker='o', markersize=10, color='k', fillstyle='none')
# %%
# Note that all four regularization operators fit the data equivalently but
# the images (i.e. how the gaps between the data points are filled) are quite
# different. This is something we should have in mind using regularization.

# %%
# Generating geostatistical media
# -------------------------------
# For generating geostatistical media, one can use the function
# generateGeostatisticalModel. It computes a correlation matrix and multiplies
# it with a pseudo-random (randn) series. The arguments are the same as for the
# correlation or constraint matrices.
model = pg.utils.generateGeostatisticalModel(mesh, I=[20, 4])
ax, cb = pg.show(mesh, model)
