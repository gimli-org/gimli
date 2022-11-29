#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
Incorporating prior data into ERT inversion
-------------------------------------------
Prior data can often help to overcome ambiguity in the inversion process.
Here we demonstrate the use of direct push (DP) data in an ERT inversion of
data collected at the surface.
"""

import matplotlib.pyplot as plt
import numpy as np
import pygimli as pg
from pygimli.physics import ert
from pygimli.frameworks import PriorModelling, JointModelling

# %%%
# The prior data
# ~~~~~~~~~~~~~~
#
# This field data is from a site with layered sands and clays over a
# resistive bedrock. We load it from the example repository.
#
# As a position of x=155m (center of the profile) we have a
# borehole/direct push with known in-situ-data. We load the three-column
# file using numpy.
#

x, z, r = pg.getExampleData("ert/bedrock.txt").T
fig, ax = plt.subplots()
ax.semilogx(r, z, "*-")
ax.set_xlabel(r"$\rho$ ($\Omega$m)")
ax.set_ylabel("depth (m)")
ax.grid(True)

# %%%
# We mainly see four layers: 1. a conductive (clayey) overburden of about
# 17m thickness, 2. a medium resistivity interbedding of silt and sand,
# about 7m thick 3. again clay with 8m thickness 4. the resistive bedrock
# with a few hundred :math:`\Omega`m
#

# %%%
# The ERT data
# ~~~~~~~~~~~~
# We load the ERT data from the example repository and plot the pseudosection.
#

data = pg.getExampleData("ert/bedrock.dat")
print(data)
ax, cb = ert.show(data)

# %%%
# The apparent resistivities show increasing values with larger spacings
# with no observable noise. We first compute geometric factors and
# estimate an error model using rather low values for both error parts.
#

data["k"] = ert.geometricFactors(data)
data["err"] = ert.estimateError(data, relativeError=0.02, absoluteUError=50e-6)
# ax, cb = ert.show(data, data["err"]*100, label="error (%)")
mgr = ert.ERTManager(data, verbose=True)
mgr.invert(paraDepth=80, quality=34.6, paraMaxCellSize=100)
# mgr.showFit()

# %%%
# For reasons of comparability, we define a unique colormap and store all
# options in a dictionary to be used in subsequent show commands.
#
# We plot the result with these and plot the DP points onto the mesh.
#

kw = dict(cMin=20, cMax=500, logScale=True, cMap="Spectral_r",
          xlabel="x (m)", ylabel="y (m)")
ax, cb = mgr.showResult(**kw)
ax.plot(x, z, ".", color="black", linestyle="None", markersize=2)
ax.grid(True)

# %%%
# We want to extract the resistivity from the mesh at the positions where
# the prior data are available. To this end, we create a list of positions
# (``pg.Pos`` class) and use a forward operator that picks the values from a
# model vector according to the cell where the position is in. See the
# regularization tutorial for details about that.
#

posVec = [pg.Pos(pos) for pos in zip(x, z)]
para = pg.Mesh(mgr.paraDomain)  # make a copy
para.setCellMarkers(pg.IVector(para.cellCount()))
fopDP = PriorModelling(para, posVec)

# %%%
# We can now use it to retrieve the model values, store it and plot it along
# with the measured values.
#

fig, ax = plt.subplots()
ax.semilogx(r, z, label="borehole")
res1 = fopDP(mgr.model)
ax.semilogx(res1, z, label="ERT")
ax.set_xlabel(r"$\rho$ ($\Omega$m)")
ax.set_ylabel("depth (m)")
ax.grid(True)
ax.legend()

# %%%
# Even though the tomogram looks interesting, the resistivity seems to
# follow a simple gradient. This is apparently a lack of resolution. Our
# assumption of an overall smooth image is wrong. Therefore we use an
# anisotropic smoothness using the vertical weighting factor ``zWeight``.
#

mgr.inv.setRegularization(zWeight=0.2)
mgr.invert()
ax, cb = mgr.showResult(**kw)
res2 = fopDP(mgr.model)

# %%%
# We observe a much more structured result with stronger vertical
# gradients that are, however, not continuous. In the middle of the
# profile we can see a short layer of increased resistivity.
#
# As alternative, we can use a geostatistic model. The vertical range can
# be well estimated from the DP data using a variogram analysis, we guess
# 5m. For the horizontal one, we can only guess a 10m higher value.
#

mgr.inv.setRegularization(2, correlationLengths=[50, 5])
mgr.invert()
ax, cb = mgr.showResult(**kw)
res3 = fopDP(mgr.model)

# %%%
# Let's compare the three resistivity soundings with the ground truth.
#

fig, ax = plt.subplots()
ax.semilogx(r, z, label="borehole")
ax.semilogx(res1, z, label="ERT smooth")
ax.semilogx(res2, z, label="ERT aniso")
ax.semilogx(res3, z, label="ERT geostat")
ax.set_xlabel(r"$\rho$ ($\Omega$m)")
ax.set_ylabel("depth (m)")
ax.grid()
ax.legend()

# %%%
# The anisotropic regularization starts to see the good conductor, but only
# the geostatistical regularization operator is able to retrieve values that
# are close to the direct push. Both show the conductor too deep.
#
# One alternative could be to use the interfaces as structural constraints in
# the neighborhood of the borehole. See ERT with structural constraints example
#
# As the DP data is not only good for comparison, we want to use its values as
# data in inversion. This is easily accomplished by taking the mapping operator
# that we already use for interpolation as a forward operator.
#
# We set up an inversion with this mesh, logarithmic transformations and
# invert the model.
#

inv = pg.Inversion(fop=fopDP, verbose=True)
inv.mesh = para
tLog = pg.trans.TransLog()
inv.modelTrans = tLog
inv.dataTrans = tLog
# inv.setRegularization(zWeight=0.2)
inv.setRegularization(correlationLengths=[50, 5])
rError = np.ones_like(r)*0.1
model = inv.run(r, rError)
ax, cb = pg.show(para, model, **kw)


# %%%
# Apparently, the geostatistical operator can be used to extrapolate
# values with given assumptions.
#

# %%%
# Joint inversion of ERT and DP data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# We now use the framework ``JointModelling`` to combine the ERT and the
# DP forward operators. So we set up a new ERT modelling operator and join
# it with ``fopDP``.
#

# fopERT = ert.ERTModelling()
# fopERT.setMesh(mesh)
# fopERT.setData(data) # not necessary as done by JointModelling
# fopJoint = JointModelling([fopERT, fopDP])
fopJoint = JointModelling([mgr.fop, fopDP])
# fopJoint.setMesh(para)
fopJoint.setData([data, pg.Vector(r)])  # needs to have .size() attribute!

# %%%
# We first test the joint forward operator. We create a modelling vector
# of constant resistivity and distribute the model response into the two
# parts that can be looked at individually.
#

model = pg.Vector(para.cellCount(), 100)
response = fopJoint(model)
respERT = response[:data.size()]
respDP = response[data.size():]
print(respDP)
ax, cb = ert.show(data, respERT)

# %%%
# The jacobian can be looked up
#

# test Jacobian
fopJoint.createJacobian(model)  # works
J = fopJoint.jacobian()
print(type(J))  # wrong type

ax, cb = pg.show(J)
print(J.mat(0))
ax, cb = pg.show(J.mat(1), markersize=4)

# %%%
# For the joint inversion, concatenate the data and error vectors, create a new
# inversion instance, set logarithmic transformations and run the inversion.
#

dataVec = np.concatenate((data["rhoa"], r))
errorVec = np.concatenate((data["err"], rError))
inv = pg.Inversion(fop=fopJoint, verbose=True)
transLog = pg.trans.TransLog()
inv.modelTrans = transLog
inv.dataTrans = transLog
inv.run(dataVec, errorVec, startModel=model)

ax, cb = pg.show(para, inv.model, **kw)

# %%%
# We have a local improvement of the model in the neighborhood of the
# borehole. Now we want to use geostatistics to get them further into the
# model.
#

inv.setRegularization(2, correlationLengths=[50, 5])
model = inv.run(dataVec, errorVec, startModel=model)
ax, cb = pg.show(para, model, **kw)

# %%%
# This model much better resembles the subsurface from all data and our
# expectations to it.
#
# We split the model response in the ERT part and the DP part. The first
# is shown as misfit.
#

respERT = inv.response[:data.size()]
ert.show(data, respERT)
misfit = - respERT / data["rhoa"] * 100 + 100
# ax, cb = ert.show(data, misfit, cMap="bwr", cMin=-5, cMax=5)

# %%%
# The second is shown as depth profile.
#

respDP = inv.response[data.size():]
fig, ax = plt.subplots()
ax.semilogx(r, z, label="borehole")
# resMesh = pg.interpolate(srcMesh=para, inVec=inv.model, destPos=posVec)
# ax.semilogx(resMesh, z, label="ERT+DP")
ax.semilogx(respDP, z, label="response")
ax.grid(True)
ax.legend()

# %%%
# The model response can much better resemble the given data compared to
# pure interpolation.
#

# %%%
# Take-away messages
# ~~~~~~~~~~~~~~~~~~
#
# -  (ERT) data inversion is highly ambiguous, particularly for hidden
#    layers
# -  prior data can help to improve regularization
# -  structural data can be of great help, but only if extended
# -  point data improve images, but only locally with smoothness
#    constraints
# -  geostatistical regularization can extrapolate point data
# -  incorporation of prior data with geostatistic regularization is best
#    and simple
#
