#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DCEM1dJointInversion
--------------------

This tutorial demonstrates how to do a classical joint inversion, i.e.
inverting two different types of data for the same physical parameter.
We use 1D versions of Direct Current (DC) and electromagnetic (EM) data that
are both sensitive to resistivity (other electromagnetic parameters neglected).
First we create a joint modelling operator, then we generate synthetic data.
We set up individual and joint inversion and compare results and uncertainties.

A ready framework for this and other joint inversions will be in pg.inversion.
"""
###############################################################################
# We import numpy numerics, mpl plotting, pygimli and the 1D plotting function
import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.mplviewer import drawModel1D, drawModel1DErr
from pygimli.inversion import modCovar, iterateBounds


# For a combined use of both methods we create a joint modelling operator that
# must be derived from the base modelling operator pg.ModellingBase.
# The individual operators are created as object members for later use.
# Another way is be setting up individual operators and passing these to init.
# The joint response is obtained by pasting (pg.cat) the individual responses.
class DCEM1dModelling(pg.ModellingBase):
    def __init__(self, nlay, ab2, mn2, freq, coilspacing, height, verbose=0):
        pg.ModellingBase.__init__(self, verbose)
        self.nlay_ = nlay
        self.fDC_ = pg.DC1dModelling(nlay, ab2, mn2, verbose)
        self.fEM_ = pg.FDEM1dModelling(nlay, freq, coilspacing, height,
                                       verbose)
        self.mesh_ = pg.createMesh1DBlock(nlay)
        self.setMesh(self.mesh_)

    def response(self, model):
        """ return combined forward response by pasting both responses."""
        return pg.cat(self.fDC_(model), self.fEM_(model))

zmax = 60
resSyn = np.array([100., 8., 30.])
thkSyn = np.array([30., 10.])
nlay = len(resSyn)
modelSyn = pg.cat(thkSyn, resSyn)  # paste thickness and resistivity
# %%
freq = [110.*2**i for i in range(10)]  # Maxmin/Promis frequencies
coilspacing = 100.0
height = 1.0
# error levels for DC (relative) and EM (absolute)
errDC = 0.02  # relative error of 2%
errEM = 2.0  # absolute error of 1 (Vs/Vp)
# EM fwd OP and synthetic model
fEM = pg.FDEM1dModelling(nlay, freq, coilspacing, height)
dataEM = fEM(modelSyn)
noiseEM = pg.RVector(len(dataEM))
pg.randn(noiseEM)
dataEM += noiseEM*errEM
# DC fwd op and synthetic
ab2 = np.logspace(0, 3, 8*3+1)
mn2 = ab2 / 3
fDC = pg.DC1dModelling(nlay, ab2, mn2)
dataDC = fDC(modelSyn)
noiseDC = pg.RVector(len(dataDC))
pg.randn(noiseDC)
dataDC *= (noiseDC * errDC + 1.)
# combined fwd op
fCO = DCEM1dModelling(nlay, ab2, mn2, freq, coilspacing, height)
dataCO = pg.cat(dataDC, dataEM)
errCO = pg.cat(dataDC*errDC, pg.RVector(len(dataEM), errEM))  # absolute error
# transformations
transThk, transRes = pg.RTransLog(), pg.RTransLog()
transEM, transDC = pg.RTrans(), pg.RTransLog()
transCO = pg.RTransCumulative()
transCO.add(transDC, len(dataDC))
transCO.add(transEM, len(dataEM))
for f in [fEM, fDC, fCO]:
    f.region(0).setTransModel(transThk)
    f.region(1).setTransModel(transRes)
# EM inversion
iEM = pg.RInversion(dataEM, fEM, transEM, False)
iEM.setAbsoluteError(errEM)
# DC inversion
iDC = pg.RInversion(dataDC, fDC, transDC, False)
iDC.setRelativeError(errDC)
# DC+EM inversion
iCO = pg.RInversion(dataCO, fCO, transCO, False)
iCO.setAbsoluteError(errCO)
# create starting model
startModel = pg.RVector(nlay*2-1, 1.0) * 40.
for i in range(nlay-1):  # thickness values
    startModel[i] = 20.
# set-up inversion schemes
for INV in [iEM, iDC, iCO]:
    INV.setMarquardtScheme(0.8)
    INV.setLambda(1000.)
    INV.setModel(startModel)
# run inversions and collect models, model bounds and uncertainties
MOD, MODL, MODU, VAR = [], [], [], []
names = ['EM', 'DC', 'DCEM']
for i, INV in enumerate([iEM, iDC, iCO]):
    print(names[i]+" inversion:")
    model = INV.run()
    INV.echoStatus()
    modelL, modelU = iterateBounds(INV, change=1.02, dchi2=1.5)
    MOD.append(model)
    MODL.append(modelL)
    MODU.append(modelU)
    var, MCMs = modCovar(INV)
    VAR.append(var)
# plot all results
xt = [1., 3., 10., 30., 100., 300.]
xtl = [str(int(xti)) for xti in xt]
fig = plt.figure(figsize=(15, 5))
for i in range(len(MOD)):
    ax = fig.add_subplot(1, len(MOD), i+1)
    ax.set_title(names[i])
    drawModel1D(ax, thkSyn, resSyn, color='k', label='syn',
                plotfunction='semilogx')
    drawModel1DErr(ax, MOD[i], MODL[i], MODU[i], label='obs')
    ax.set_xlim(1., 300.)
    ax.set_ylim(zmax, 0.)
    ax.set_xticks(xt)
    ax.set_xticklabels(xtl)
    ax.grid(True)
    ax.legend(numpoints=1, loc='upper center')

    thk = MOD[i][:nlay-1]
    zl = np.cumsum(thk)
    zvec = np.hstack((zl, zl-thk/2, zl[-1]+thk[-1]/2))
    for j in range(nlay*2-1):
        v = np.log(1.+VAR[i][j])
        if j < nlay-1:
            plt.text(1., zvec[j], '$\delta$='+str(np.round_(v*thk[j], 1))+'m')
        else:
            plt.text(1., zvec[j], '$\delta$='+str(np.round_(v*100., 1))+'$\%$')

plt.show()
# We clearly see that DC resolves the near-surface layer better, but EM sees
# the conductor better. Both methods overestimate the thickness of the second
# layer. The joint inversion combines the advantages and is closest to reality.
# Günther, T. (2014): On Inversion of Frequency Domain Electromagnetic Data in
# Salt Water Problems - Sensitivity and Resolution. Ext. Abstr. EAGE Near
# Surface 2013, doi:10.3997/2214-4609.20131387.
