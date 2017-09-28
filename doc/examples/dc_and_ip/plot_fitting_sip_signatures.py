#!/usr/bin/env python
# encoding: utf-8

r"""
Fitting SIP signatures
----------------------

This example highlights some of the capabilities of pyGimli to analyze spectral
induced polarization (SIP) signatures."""

###############################################################################
# Given is a SIP signature for N frequencies in the form of resistivity
# magnitude and phase values

from pygimli.physics import SIPSpectrum
from pygimli.physics.SIP import ColeColeRho
import numpy as np
import pygimli as pg

f = np.logspace(-2, 5, 100)
m = np.linspace(0.1, 0.9, 5)
tau = 0.01

# generate synthetic data to fit
Z = ColeColeRho(f, rho=100, m=m[0], tau=tau, c=0.5)
amplitude = np.abs(Z)
phase = np.angle(Z)

sip = SIPSpectrum(
    f=f,
    amp=np.abs(Z),
    # note the minus sign: we need to provide -phase[rad]
    phi=-np.angle(Z),
)
sip.showData()  # znorm=True)
sip.showDataKK()  # check Kramers-Kronig relations
if True:  # determine constant epsilon and remove i omega term
    sip.removeEpsilonEffect()
    sip.showDataKK()  # check Kramers-Kronig again after removing epsilon
    sip.fitColeCole(useCond=False)  # works for both rho and sigma models
else:  # classical Pelton approach
    sip.fitCCEM()  # fit an SIP Cole-Cole term and an EM term (also Cole-Cole)

# %% create titles and plot data, fit and model
sip.showAll()  # save=True)
# %%
sip.fitDebyeModel(new=True, showFit=True)
pg.wait()
