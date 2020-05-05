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
np.random.seed(1337)

import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.viewer.mpl import drawModel1D


###############################################################################
# First, we define a modelling class that calls two other classes and pastes
# their results to one vector.
class DCEM1dModelling(pg.core.ModellingBase):
    """Modelling jointing DC and EM 1Dforward operators."""

    def __init__(self, nlay, ab2, mn2, freq, coilspacing, verbose=False):
        """Init number of layers, AB/2, MN/2, frequencies & coil spacing."""
        pg.core.ModellingBase.__init__(self, verbose)
        self.nlay_ = nlay
        self.fDC_ = pg.core.DC1dModelling(nlay, ab2, mn2, verbose)
        self.fEM_ = pg.core.FDEM1dModelling(nlay, freq, coilspacing, verbose)
        self.mesh_ = pg.meshtools.createMesh1DBlock(nlay)
        self.setMesh(self.mesh_)

    def response(self, model):
        """Return concatenated response of DC and EM FOPs."""
        return pg.cat(self.fDC_(model), self.fEM_(model))


###############################################################################
# The actual script starts here. There are some options to play with
noiseEM = 1.  # absolute (per cent of primary signal)
noiseDC = 3.  # in per cent
lamEM, lamDC, lamDCEM = 300., 500., 500.  # regularization strength
verbose = False

###############################################################################
# First we create a synthetic model.

nlay = 3  # number of layers
thk = pg.Vector(nlay - 1, 15.0)  # 15m thickness each
res = pg.Vector(nlay, 200.0)  # 200 Ohmm
res[1] = 10.
res[2] = 50.
model = pg.cat(thk, res)  # paste together to one model

###############################################################################
# We first set up EM forward operator and generate synthetic data with noise

coilspacing = 50.
nf = 10
freq = pg.Vector(nf, 110.)
for i in range(nf-1):
    freq[i+1] = freq[i] * 2.

fEM = pg.core.FDEM1dModelling(nlay, freq, coilspacing)
dataEM = fEM(model)
for i in range(len(dataEM)):
    dataEM[i] += pg.randn(1)[0] * noiseEM

###############################################################################
# We define model transformations: logarithms and log with upper+lower bounds

transRhoa = pg.trans.TransLog()
transThk = pg.trans.TransLog()
transRes = pg.trans.TransLogLU(1., 1000.)
transEM = pg.trans.Trans()
fEM.region(0).setTransModel(transThk)
fEM.region(1).setTransModel(transRes)

###############################################################################
# We set up the independent EM inversion and run the model.

invEM = pg.core.Inversion(dataEM, fEM, transEM, True, True)
modelEM = pg.Vector(nlay * 2 - 1, 50.)
invEM.setModel(modelEM)
invEM.setAbsoluteError(noiseEM)
invEM.setLambda(lamEM)
invEM.setMarquardtScheme(0.9)
modelEM = invEM.run()
respEM = invEM.response()

###############################################################################
# Next we set up the DC forward operator and generate synthetic data with noise

ab2 = pg.Vector(20, 3.)
na = len(ab2)
mn2 = pg.Vector(na, 1.0)
for i in range(na-1):
    ab2[i+1] = ab2[i] * 1.3
fDC = pg.core.DC1dModelling(nlay, ab2, mn2)
dataDC = fDC(model)
for i in range(len(dataDC)):
    dataDC[i] *= 1. + pg.randn(1)[0] * noiseDC / 100.

fDC.region(0).setTransModel(transThk)
fDC.region(1).setTransModel(transRes)

# We set up the independent DC inversion and let it run.
invDC = pg.core.Inversion(dataDC, fDC, transRhoa, verbose)
modelDC = pg.Vector(nlay*2-1, 20.)
invDC.setModel(modelDC)
invDC.setRelativeError(noiseDC/100.)
invDC.setLambda(lamDC)
invDC.setMarquardtScheme(0.9)
modelDC = invDC.run()
respDC = invDC.response()

###############################################################################
# Next we create a the joint forward operator (see class above).

fDCEM = DCEM1dModelling(nlay, ab2, mn2, freq, coilspacing)
fDCEM.region(0).setTransModel(transThk)
fDCEM.region(1).setTransModel(transRes)

###############################################################################
# We setup the joint inversion combining, transformations, data and errors.

transData = pg.trans.TransCumulative()
transData.add(transRhoa, na)
transData.add(transEM, nf*2)
invDCEM = pg.core.Inversion(pg.cat(dataDC, dataEM), fDCEM, transData, verbose)
modelDCEM = pg.Vector(nlay * 2 - 1, 20.)
invDCEM.setModel(modelDCEM)
err = pg.cat(dataDC * noiseDC / 100., pg.Vector(len(dataEM), noiseEM))
invDCEM.setAbsoluteError(err)
invDCEM.setLambda(lamDCEM)
invDCEM.setMarquardtScheme(0.9)
modelDCEM = invDCEM.run()
respDCEM = invDCEM.response()

###############################################################################
# The results of the inversion are plotted for comparison.

for inv in [invEM, invDC, invDCEM]:
    inv.echoStatus()
print([invEM.chi2(), invDC.chi2(), invDCEM.chi2()])  # chi-square values

###############################################################################
# %% We finally plot the results
fig = plt.figure(1, figsize=(10, 5))
ax1 = fig.add_subplot(131)
drawModel1D(ax1, thk, res, plot='semilogx', color='blue')
drawModel1D(ax1, modelEM(0, nlay-1), modelEM(nlay-1, nlay*2-1), color='green')
drawModel1D(ax1, modelDC(0, nlay-1), modelDC(nlay-1, nlay*2-1), color='cyan')
drawModel1D(ax1, modelDCEM(0, nlay-1), modelDCEM(nlay-1, nlay*2-1),
            color='red')
ax1.legend(('syn', 'EM', 'DC', 'JI'))
ax1.set_xlim((10., 1000.))
ax1.set_ylim((40., 0.))
ax1.grid(which='both')
ax2 = fig.add_subplot(132)
ax2.semilogy(dataEM(0, nf), freq, 'bx', label='syn IP')
ax2.semilogy(dataEM(nf, nf*2), freq, 'bo', label='syn OP')
ax2.semilogy(respEM(0, nf), freq, 'g--', label='EM')
ax2.semilogy(respEM(nf, nf*2), freq, 'g--')
ax2.semilogy(respDCEM(na, na+nf), freq, 'r:', label='DCEM')
ax2.semilogy(respDCEM(na+nf, na+nf*2), freq, 'r:')
ax2.set_ylim((min(freq), max(freq)))
ax2.set_xlabel("IP/OP in %")
ax2.set_ylabel("$f$ in Hz")
ax2.yaxis.set_label_position("right")
ax2.grid(which='both')
ax2.legend(loc="best")
ax3 = fig.add_subplot(133)
ax3.loglog(dataDC, ab2, 'bx-', label='syn')
ax3.loglog(respDC, ab2, 'c-', label='DC')
ax3.loglog(respDCEM(0, na), ab2, 'r:', label='DCEM')
# ax3.axis('tight')
ax3.set_ylim((max(ab2), min(ab2)))
ax3.grid(which='both')
ax3.set_xlabel(r"$\rho_a$ in $\Omega$m")
ax3.set_ylabel("AB/2 in m")
ax3.yaxis.set_ticks_position("right")
ax3.yaxis.set_label_position("right")
ax3.legend(loc="best")
pg.wait()

###############################################################################

# Günther, T. (2013): On Inversion of Frequency Domain Electromagnetic Data in
# Salt Water Problems - Sensitivity and Resolution. Ext. Abstr., 19th European
# Meeting of Environmental and Engineering Geophysics, Bochum, Germany.
