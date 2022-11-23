#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
Regularization - concepts explained
-----------------------------------

In geophysical inversion, we minimize the data objective functional as
the L2 norm of the misfit between data $d$ and the forward
response $f$ of the model $m$, weighted by the data error
$\epsilon$:

.. math:: \Phi_d = \sum\limits_i^N \left(\frac{d_i-f_i(m)}{\epsilon_i}\right)^2=\|W_d(d-f(m))\|^2

As this minimization problem is non-unique and ill-posed, we introduce a
regularization term $\Phi$, weighted by a regularization parameter
$\lambda$:

.. math:: \Phi = \Phi_d + \lambda \Phi_m

The regularization strength $\lambda$ should be chosen so that the
data are fitted within noise, i.e. $\chi^2=\Phi_d/N=1$.

In the term $\Phi-m$ we put our expectations to the model, e.g. to
be close to any prior model. In many cases we do not have much
information and aim for the smoothest model that is able to fit our
data. We decribe it by the operator $W_m$:

.. math:: \Phi_m=\|W_m (m-m_{ref})\|^2

The regularization operator is defined by some constraint operator
$C$ weighted by some weighting function $w$ so that
$W_m=\mbox{diag}(w) C$. The operator $C$ can be a discrete
smoothness operator, or the identity to keep the model close to the
reference model $m_{ref}$.

"""

# %%%
# We start with importing the numpy, matplotlib and pygimli libraries
#

import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.core.math import symlog

# %%%
# Regularisation drives the model where the data are too weak to constrain
# the model. In order to explain different kinds of regularization (also
# called constraints), we use a very simple mapping forward operator: The
# values at certain positions are picked.
#

from pygimli.frameworks import PriorModelling

# %%%
# Implementation 1. determine the indices where the cells are
#
# ::
#
#    ind = [mesh.findCell(po).id() for po in pos]
#
# 2. forward response: take the model at indices
#
# ::
#
#    response = model[ind]
#
# 3. Jacobian matrix
#
# ::
#
#    J = pg.SparseMapMatrix()
#    J.resize(len(ind), mesh.cellCount())
#    for i, n in enumerate(self.ind):
#        self.J.setVal(i, n, 1.0)
#

# %%%
# We exemplify this on behalf of a simple triangular mesh in a rectangular
# domain.
#

rect = mt.createRectangle(start=[0, -10], end=[10, 0])
mesh = mt.createMesh(rect, quality=34.5, area=0.3)
print(mesh)
pg.show(mesh);

# %%%
# We define two positions where we associate two arbitrary values.
#

pos = [[3, -3], [7, -7]]
vals = np.array([20., 15.])
fop = PriorModelling(mesh, pos)

# %%%
# We set up an inversion instance with the forward operator and prepare
# the keywords for running the inversion always the same way: - the data
# vector - the error vector (as relative error) - a starting model value
# (could also be vector)
#

inv = pg.Inversion(fop=fop, verbose=False)
invkw = dict(dataVals=vals, errorVals=np.ones_like(vals)*0.03, startModel=10)

# %%%
# Classical smoothness constraints
# --------------------------------
#

inv.setRegularization(cType=1)  # the default
result = inv.run(**invkw)
pg.show(mesh, result);

# %%%
# We will have a closer look at the regularization matrix C.
#

C = fop.constraints()
print(C.rows(), C.cols(), mesh)
ax, _ = pg.show(fop.constraints(), markersize=1);

row = C.row(111)
nz = np.nonzero(row)[0]
print(nz, row[nz])

inv.setRegularization(cType=1, zWeight=0.2)  # the default
result = inv.run(**invkw)
pg.show(mesh, result);

RM = fop.regionManager()
cw = RM.constraintWeights()
print(min(cw), max(cw))

# %%%
# Now we try some other regularization options.
#

inv.setRegularization(cType=0)  # damping difference to starting model
result = inv.run(**invkw)
pg.show(mesh, result);

# %%%
# Obviously, the damping keeps the model small ($\log 1=0$) as the
# starting model is NOT a reference model by default. We will enable this
# by specifying the ``isReference`` switch.
#

invkw["isReference"] = True
result = inv.run(**invkw)
pg.show(mesh, result);

inv.setRegularization(cType=10)  # mix of 1st order smoothing and damping
result = inv.run(**invkw)
pg.show(mesh, result);

C = fop.constraints()
print(C.rows(), C.cols())
print(mesh)
pg.show(fop.constraints(), markersize=1);

# %%%
# We see that we have the first order smoothness and the identity matrix
# below each other. We can also use a second-order (-1 2 -1) smoothness
# operator by ``cType=2``.
#

inv.setRegularization(cType=2)  # 2nd order smoothing
result = inv.run(**invkw)
pg.show(mesh, result);

# %%%
# We have a closer look at the constraints matrix
#

C = fop.constraints()
print(C.rows(), C.cols(), mesh)
pg.show(C, markersize=1);

# %%%
# It looks like a Laplace operator and seems to have a wider range
# compared to first-order smoothness.
#

# %%%
# Geostatistic regularization
# ---------------------------
#

# %%%
# The idea is that not only neighbors are correlated to each other but to
# have a wider correlation by using an operator
#
# .. math::
#
#        \textbf{C}_{\text{M},ij}=\sigma^{2}\exp{\left(
#            -\sqrt{
#            \left(\frac{\textbf{H}^x_{ij}}{I_{x}}\right)^{2}+
#            \left(\frac{\textbf{H}^y_{ij}}{I_{y}}\right)^{2}
#    }\right)}.
#
# More details can be found in
# https://www.pygimli.org/_tutorials_auto/3_inversion/plot_6-geostatConstraints.html
#

ind = mesh.findCell([5, -5]).id()
print(ind)
CM = pg.utils.covarianceMatrix(mesh, I=[8, 4], dip=-20)
pg.show(mesh, CM[:, ind], cMap="magma_r");

# %%%
# We have a look at the matrix and multiply it with a zero vector of just one 1
#

C = pg.matrix.Cm05Matrix(CM)
print(C)

vec = pg.Vector(mesh.cellCount())
vec[ind] = 1.0
pg.show(mesh, C*vec, cMin=-1, cMax=1, cMap="bwr");

# %%%
# For displaying the wide range of magnitudes we use the symlog function
#

pg.show(mesh, symlog(C*vec, 1e-2), cMin=-2, cMax=2, cMap="bwr");

# %%%
# For comparison, we use a much finer mesh
#

fineMesh = mt.createMesh(rect, area=0.03)
pg.show(fineMesh);

# %%%
# and do the same
#

Cfine = pg.matrix.GeostatisticConstraintsMatrix(
    mesh=fineMesh, I=[8, 4], dip=-20)
indFine = fineMesh.findCell([5, -5]).id()
vec = pg.Vector(fineMesh.cellCount())
vec[indFine] = 1.0
pg.show(fineMesh, symlog(Cfine*vec, 1e-2), cMin=-1, cMax=1, cMap="bwr");

# %%%
# Application
# ~~~~~~~~~~~
# We now pass the correlation length directly to the inversion instance
#

inv.setRegularization(correlationLengths=[2, 2, 2])
result = inv.run(**invkw)
pg.show(mesh, result);

# %%%
# This look structurally similar to the second-order smoothness, but can
# drive values outside the expected range in regions of no data coverage.
#

inv.setRegularization(correlationLengths=[2, 0.5, 2])
result = inv.run(**invkw)
pg.show(mesh, result);

# %%%
# We change the dip to be inclining
#

inv.setRegularization(correlationLengths=[2, 0.5, 2], dip=-20)
result = inv.run(**invkw)
pg.show(mesh, result);

# %%%
# We now add many more points.
#

N = 30
x = np.random.rand(N) * 10
y = -np.random.rand(N) * 10
v = np.random.rand(N)*10 + 10
plt.plot(x, y, "*")

# %%%
# and repeat the above computations
#

fop = PriorModelling(mesh, zip(x, y))
inv = pg.Inversion(fop=fop, verbose=True)
inv.setRegularization(correlationLengths=[4, 4])
result = inv.run(v, np.ones_like(v)*0.03, startModel=10)
pg.show(mesh, result);

# %%%
# Comparing the data with the model response is always a good idea.
#

plt.plot(v, inv.response, "*");

# %%%
# Individual regularization operators
# -----------------------------------
#
# Say you want to combine geostatistic operators with a damping, you can
# create a block matrix pasting the matric vertically.
#

C = pg.matrix.BlockMatrix()
G = pg.matrix.GeostatisticConstraintsMatrix(mesh=mesh, I=[2, 0.5], dip=-20)
I = pg.matrix.IdentityMatrix(mesh.cellCount(), val=0.1)
C.addMatrix(G, 0, 0)
C.addMatrix(I, mesh.cellCount(), 0)
pg.show(C);

# %%%
# Note that in `pg.matrix` you find a lot of matrices and matrix generators.
#
# We set this matrix directly and do the inversion.
#

fop.setConstraints(C)
result = inv.run(v, np.ones_like(v)*0.03, startModel=10, isReference=True)
pg.show(mesh, result);

# %%%
# If you are using a method manager, you access the inversion instance by
# `mgr.inv` and the forward operator by `mgr.fop`.
#

# %%%
# Take-away messages
# ==================
#
# -  regularization drives the model where data are weak
# -  think and play with your assumptions to the model
# -  there are several predefined options
# -  geostatistic regularization can be superior
#
#    -  it is mesh-independent
#    -  it better fills the data gaps (e.g. 3D inversion of 2D profiles)
#
