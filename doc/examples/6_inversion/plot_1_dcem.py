#!/usr/bin/env python
# encoding: utf-8
"""
DC-EM Joint inversion
---------------------

This is an old script from early pyGIMLi jointly inverting direct current (DC)
and electromagnetic (EM) soundings on the modelling abstraction level.
Note that this is not recommended as a basis for programming, because there
is a dedicated framework for classical joint inversion. However, it explains
what happens under the hood in the much simpler script that follows."""

###############################################################################
# The case has been documented by :cite:`Guenther2013NSG`.
import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.physics.em import HEMmodelling
from pygimli.physics.ves import VESModelling
from pygimli.frameworks import JointModelling
from pygimli.frameworks import MarquardtInversion
from pygimli.viewer.mpl import drawModel1D


###############################################################################
# The actual script starts here. There are some options to play with
errorEMabs = 1.  # absolute (ppm of primary signal)
errorDCrel = 3.  # in per cent

###############################################################################
# First we create a synthetic model.

nlay = 4  # number of layers
synThk = [5, 15, 15]
synRes = [1000, 100, 500, 20]

# %%%
# We first set up EM forward operator and generate synthetic data with noise
#

# MaxMin-10/Promys instrument, 1m above the ground
nf = 10
freq = 2**np.arange(nf) * 110.
fEM = HEMmodelling(nlay=nlay, height=1, f=freq, r=100, scaling="%")
# alternative: RESOLVE airborne system
# freq = [338]
# fEM = HEMmodelling(nlay=nlay, height=50, f=freq, r=10)

dataEM = fEM(synThk + synRes)
errorEM = np.ones_like(dataEM) * errorEMabs
dataEM += pg.randn(len(dataEM), seed=1234) * errorEM
print(dataEM)

# %%%
# We define model transformations: logarithms and log with upper+lower bounds
#

# transRhoa = pg.trans.TransLog()
# transEM = pg.trans.Trans()

###############################################################################
# We set up the independent EM inversion and run the model.
startModel = [10]*(nlay-1) + [100]*nlay

invEM = MarquardtInversion(fop=fEM, verbose=False)
# invEM.setRegularization(0, limits=[1, 30], cType=0)
# invEM.setRegularization(1, limits=[10, 1000])
modelEM = invEM.run(dataEM, np.abs(errorEM/dataEM), startModel=startModel)

# %%%
# Next we set up the DC forward operator and generate synthetic data with noise
#

ab2 = pg.Vector(20, 3.)
na = len(ab2)
mn2 = pg.Vector(na, 1.0)
for i in range(na-1):
    ab2[i+1] = ab2[i] * 1.3
fDC = VESModelling(ab2=ab2, mn2=mn2)
dataDC = fDC(synThk+synRes)
errorDC = np.ones_like(dataDC) * errorDCrel / 100.
dataDC *= 1. + pg.randn(len(dataDC), seed=1234) * errorDC

# fDC.region(0).setTransModel(transThk)
# fDC.region(1).setTransModel(transRes)

# We set up the independent DC inversion and let it run.
invDC = MarquardtInversion(fop=fDC, verbose=False)
# invDC.dataTrans = pg.trans.TransLog()
# invDC.setRegularization(0, limits=[1, 100], cType=0)
# invDC.setRegularization(1, limits=[1, 1000], cType=0)
modelDC = invDC.run(dataDC, errorDC, startModel=startModel)

# %%%
# We create a the joint forward operator using the Joint inversion framework.
#

fDCEM = JointModelling([fDC, fEM])
fDCEM.setData([dataDC, dataEM])  # just for sizes!
# transDCEM = pg.trans.TransCumulative()
# transDCEM.add(transRhoa, na)
# transDCEM.add(transEM, nf*2)
jointData = pg.cat(dataDC, dataEM)
jointError = pg.cat(errorDC, np.abs(errorEM/dataEM))
invDCEM = MarquardtInversion(fop=fDCEM, verbose=False)
# invDCEM.setRegularization(0, limits=[1, 100], cType=0)
# invDCEM.setRegularization(1, limits=[1, 1000], cType=0)
modelDCEM = invDCEM.run(jointData, jointError,
                        startModel=startModel)

# %%%
# The results of the inversion are plotted for comparison.
#

for inv in [invEM, invDC, invDCEM]:
    inv.echoStatus()

print([invEM.chi2(), invDC.chi2(), invDCEM.chi2()])  # chi-square values

# %%%
# We finally plot the results along with data and model responses.
#

fig, (ax1, ax2, ax3) = plt.subplots(figsize=(10, 5), ncols=3)
drawModel1D(ax1, synThk, synRes, plot='semilogx', color='C0', label="synth")
drawModel1D(ax1, model=modelEM, color='C1', label="DC")
drawModel1D(ax1, model=modelDC, color='C2', label="EM")
drawModel1D(ax1, model=modelDCEM, color='C3', label="DC-EM")
ax1.legend()
ax1.set_xlim((10., 1000.))
ax1.set_ylim((40., 0.))
ax1.grid(which='both')
ax2.semilogy(dataEM[0:nf], freq, 'x', color="C0", label='syn IP')
ax2.semilogy(dataEM[nf:nf*2], freq, 'o', color="C0", label='syn OP')
ax2.semilogy(invEM.response[0:nf], freq, '--', color="C2", label='EM')
ax2.semilogy(invEM.response[nf:nf*2], freq, '--', color="C2")
ax2.semilogy(invDCEM.response[na:na+nf], freq, ':', color="C2", label='DCEM')
ax2.semilogy(invDCEM.response[na+nf:na+nf*2], freq, '2:', color="C2")
ax2.set_ylim((min(freq), max(freq)))
ax2.set_xlabel("IP/OP in %")
ax2.set_ylabel("$f$ in Hz")
ax2.yaxis.set_label_position("right")
ax2.grid(which='both')
ax2.legend(loc="best")
ax3.loglog(dataDC, ab2, 'bx-', label='syn')
ax3.loglog(invDC.response, ab2, '-', label='DC', color="C1")
ax3.loglog(invDCEM.response[0:na], ab2, '-', label='DCEM', color="C3")
ax3.set_ylim((max(ab2), min(ab2)))
ax3.grid(which='both')
ax3.set_xlabel(r"$\rho_a$ in $\Omega$m")
ax3.set_ylabel("AB/2 in m")
ax3.yaxis.set_ticks_position("right")
ax3.yaxis.set_label_position("right")
ax3.legend()

###############################################################################

# Günther, T. (2013): On Inversion of Frequency Domain Electromagnetic Data in
# Salt Water Problems - Sensitivity and Resolution. Ext. Abstr., 19th European
# Meeting of Environmental and Engineering Geophysics, Bochum, Germany.
