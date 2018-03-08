#!/usr/bin/env python
# encoding: utf-8

r"""
Fitting SIP signatures
----------------------

This example highlights some of the capabilities of pyGimli to analyze spectral
induced polarization (SIP) signatures."""

###############################################################################
# Generate synthetic data and initialize an SIPSpectrum object
from pygimli.physics.SIP import ColeColeRho
from pygimli.physics import SIPSpectrum
import numpy as np
import pygimli as pg

f = np.logspace(-2, 5, 100)

# term1
Z1 = ColeColeRho(f, rho=1, m=0.1, tau=0.5, c=0.5)
# term2
Z2 = ColeColeRho(f, rho=1, m=0.25, tau=1e-6, c=1.0)
# create sum
rho0 = 100 # (Ohm m)
Z = rho0 * (Z1 + Z2)
amplitude = np.abs(Z)
phase = np.angle(Z)

sip = SIPSpectrum(f=f,
                  amp=np.abs(Z),
                  # note the minus sign: we need to provide -phase[rad]
                  phi=-np.angle(Z))
sip.showData()
sip.showDataKK()  # check Kramers-Kronig relations

###############################################################################
# Fit a Cole-Cole model
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

sip = SIPSpectrum(f=f,
                  amp=np.abs(Z),
                  # note the minus sign: we need to provide -phase[rad]
                  phi=-np.angle(Z),
)
sip.fitColeCole(useCond=False)  # works for both rho and sigma models
sip.showAll()  # save=True)
pg.wait()

###############################################################################
# Fit a double Cole-Cole model
from pygimli.physics.SIP import ColeColeRho
import numpy as np
import pygimli as pg

f = np.logspace(-2, 5, 100)

# term1
Z1 = ColeColeRho(f, rho=1, m=0.1, tau=0.5, c=0.5)
# term2
Z2 = ColeColeRho(f, rho=1, m=0.25, tau=1e-6, c=1.0)
# create sum
rho0 = 100 #(Ohm m)
Z = rho0 * (Z1 + Z2)
amplitude = np.abs(Z)
phase = np.angle(Z)

sip = SIPSpectrum(
    f=f,
    amp=np.abs(Z),
    # note the minus sign: we need to provide -phase[rad]
    phi=-np.angle(Z),
)
sip.fitCCEM()  # fit an SIP Cole-Cole term and an EM term (also Cole-Cole)

# %% create titles and plot data, fit and model
sip.showAll()  # save=True)
# %%
pg.wait()

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
# sip.showData()  # znorm=True)
# sip.showDataKK()  # check Kramers-Kronig relations
# if True:  # determine constant epsilon and remove i omega term
#     sip.removeEpsilonEffect()
#     sip.showDataKK()  # check Kramers-Kronig again after removing epsilon
#   # sip.fitColeCole(useCond=False)  # works for both rho and sigma models

# %% create titles and plot data, fit and model
sip.showAll()  # save=True)
# %%
sip.fitDebyeModel(new=True, showFit=True)
pg.wait()
