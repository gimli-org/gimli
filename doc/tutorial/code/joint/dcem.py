#!/usr/bin/env python
# -*- coding: utf-8 -*-
# this is a redundant file but it will be recycled for the examples
import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.viewer.mpl import drawModel1D


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


# options to play with
noiseEM = 1.  # absolute (per cent of primary signal)
noiseDC = 3.  # in per cent
lamEM, lamDC, lamDCEM = 300., 500., 500.
verbose = False

# synthetic model
nlay = 3
model = pg.Vector(nlay * 2 - 1)
for i in range(nlay - 1):
    model[i] = 15.

for i in range(nlay):
    model[i+nlay-1] = 200.

model[nlay] = 10.
model[nlay + 1] = 50.
thk = model(0, nlay-1)
res = model(nlay-1, nlay*2-1)

# EM forward operator and synthetic data
coilspacing = 50.
nf = 10
freq = pg.Vector(nf, 110.)
for i in range(nf-1):
    freq[i+1] = freq[i] * 2.

fEM = pg.core.FDEM1dModelling(nlay, freq, coilspacing)
dataEM = fEM(model)
for i in range(len(dataEM)):
    dataEM[i] += np.random.randn(1)[0] * noiseEM

# model transformations
transRhoa = pg.trans.TransLog()
transThk = pg.trans.TransLog()
transRes = pg.trans.TransLogLU(1., 1000.)
transEM = pg.trans.Trans()
fEM.region(0).setTransModel(transThk)
fEM.region(1).setTransModel(transRes)

# independent EM inversion
invEM = pg.Inversion(dataEM, fEM, transEM, verbose)
modelEM = pg.Vector(nlay * 2 - 1, 50.)
invEM.setModel(modelEM)
invEM.setAbsoluteError(noiseEM)
invEM.setLambda(lamEM)
invEM.setMarquardtScheme(0.9)
modelEM = invEM.run()
respEM = invEM.response()

# DC forward operator and synthetic data
ab2 = pg.Vector(20, 3.)
na = len(ab2)
mn2 = pg.Vector(na, 1.0)
for i in range(na-1):
    ab2[i+1] = ab2[i] * 1.3
fDC = pg.core.DC1dModelling(nlay, ab2, mn2)
dataDC = fDC(model)
for i in range(len(dataDC)):
    dataDC[i] *= 1. + np.random.randn(1)[0] * noiseDC / 100.

fDC.region(0).setTransModel(transThk)
fDC.region(1).setTransModel(transRes)

# independent DC inversion
invDC = pg.Inversion(dataDC, fDC, transRhoa, verbose)
modelDC = pg.Vector(nlay*2-1, 20.)
invDC.setModel(modelDC)
invDC.setRelativeError(noiseDC/100.)
invDC.setLambda(lamDC)
invDC.setMarquardtScheme(0.9)
modelDC = invDC.run()
respDC = invDC.response()

# joint forward operator
fDCEM = DCEM1dModelling(nlay, ab2, mn2, freq, coilspacing)
fDCEM.region(0).setTransModel(transThk)
fDCEM.region(1).setTransModel(transRes)

# joint inversion
transData = pg.trans.TransCumulative()
transData.add(transRhoa, na)
transData.add(transEM, nf*2)
invDCEM = pg.Inversion(pg.cat(dataDC, dataEM), fDCEM, transData, verbose)
modelDCEM = pg.Vector(nlay * 2 - 1, 20.)
invDCEM.setModel(modelDCEM)
err = pg.cat(dataDC * noiseDC / 100., pg.Vector(len(dataEM), noiseEM))
invDCEM.setAbsoluteError(err)
invDCEM.setLambda(lamDCEM)
invDCEM.setMarquardtScheme(0.9)
modelDCEM = invDCEM.run()
respDCEM = invDCEM.response()

# comparison
invDC.echoStatus()
invEM.echoStatus()
invDCEM.echoStatus()
[invDC.chi2(), invEM.chi2(), invDCEM.chi2()]

# plot results
fig = plt.figure(1)
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
ax2.semilogy(dataEM(0, nf), freq, 'bx', dataEM(nf, nf*2), freq, 'bo')
ax2.semilogy(respEM(0, nf), freq, 'g--', respEM(nf, nf*2), freq, 'g--')
ax2.semilogy(respDCEM(na, na+nf), freq, 'r:',
             respDCEM(na+nf, na+nf*2), freq, 'r:')
ax2.set_ylim((min(freq), max(freq)))
ax2.grid(which='both')
ax2.legend(("syn", "", "EM", "", "DCEM", ""), loc="lower left")
ax3 = fig.add_subplot(133)
ax3.loglog(dataDC, ab2, 'bx-', label='syn')
ax3.loglog(respDC, ab2, 'c-', label='DC')
ax3.loglog(respDCEM(0, na), ab2, 'r:', label='EM')
# ax3.axis('tight')
ax3.set_ylim((max(ab2), min(ab2)))
ax3.grid(which='both')
ax3.set_xlabel(r"\rho_a in \Omegam")
ax3.set_ylabel("AB/2 in m")
ax3.legend(loc="best")
ax3.legend(("syn", "DC", "EMDC"), loc="lower left")
plt.show()
