#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DC1dBlock
---------

This tutorial shows how an built-in forward operator is used for inversion.
A DC 1D (VES) modelling is used to generate data, noisify and invert them."""

###############################################################################
# We import numpy, matplotlib and the 1D plotting function
import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.mplviewer import drawModel1D

###############################################################################
# some definitions before
nlay = 4  # number of layers
lam = 200.  # (initial) regularization parameter
errPerc = 3.  # relative error of 3 percent
###############################################################################
# generate current and potential length vector
ab2 = np.logspace(-1, 2, 50)
mn2 = ab2 / 3.
###############################################################################
# initialize the forward modelling operator
f = pg.DC1dModelling(nlay, ab2, mn2)
###############################################################################
# other ways are by specifying a Data Container or am/an/bm/bn distances
synres = [100., 500., 20., 800.]  # synthetic resistivity
synthk = [0.5, 3.5, 6.]  # synthetic thickness (lay layer is infinite)
###############################################################################
# the forward operator can be called by f.response(model) or simply f(model)
rhoa = f(synthk+synres)
rhoa = rhoa * (pg.randn(len(rhoa)) * errPerc / 100. + 1.)
###############################################################################
# create some transformations used for inversion
transThk = pg.RTransLog()  # log-transform ensures thk>0
transRho = pg.RTransLogLU(1, 1000)  # lower and upper bound
transRhoa = pg.RTransLog()  # log transformation also for data
###############################################################################
# set model transformation for thickness and resistivity
f.region(0).setTransModel(transThk)  # 0=thickness
f.region(1).setTransModel(transRho)  # 1=resistivity
###############################################################################
# generate start model from median app. resistivity & current bipole length
paraDepth = max(ab2) / 3.  # rule-of-thumb for Wenner/Schlumberger
f.region(0).setStartValue(paraDepth / nlay / 2)
f.region(1).setStartValue(np.median(rhoa))
###############################################################################
# set up inversion
inv = pg.RInversion(rhoa, f, transRhoa, True)  # data vector, fop, verbose
###############################################################################
# could also be set by inv.setTransData(transRhoa)
inv.setRelativeError(errPerc / 100.0)
inv.setLambda(lam)  # (initial) regularization parameter
inv.setMarquardtScheme(0.9)  # decrease lambda by factor 0.9
model = f.createStartVector()  # creates from region start value
###############################################################################
# optionally change default model by changing a layer resistivity
model[nlay] *= 1.5
inv.setModel(model)  #
###############################################################################
# run actual inversion
model = inv.run()  # result is a pg.RVector, but compatible to numpy array
res, thk = model[nlay-1:nlay*2-1], model[0:nlay-1]
###############################################################################
# show everything
fig, ax = plt.subplots(ncols=2, figsize=(8, 6))  # two-column figure
# plot model (inverted and synthetic)
drawModel1D(ax[0], thk, res, color='r')  # r'\rho in \Omega m')
drawModel1D(ax[0], synthk, synres, color='b')


print(thk)

ax[0].grid(True, which='both')
# plot sounding curve data and model response
ax[1].loglog(rhoa, ab2, 'rx-', label='measured')
ax[1].loglog(inv.response(), ab2, 'b-', label='fitted')
ax[1].set_ylim((max(ab2), min(ab2)))
ax[1].grid(True, which='both')
ax[1].set_xlabel(r'$\rho_a$ [$\Omega$m]')
ax[1].set_ylabel('AB/2 [m]')
ax[1].legend(loc='best')

plt.show()
