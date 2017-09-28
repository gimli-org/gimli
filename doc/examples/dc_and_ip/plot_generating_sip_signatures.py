#!/usr/bin/env python
# encoding: utf-8

r"""
Generating SIP signatures
-------------------------

This example highlights some of the capabilities of pyGimli to generate
spectral induced polarization (SIP) signatures."""

###############################################################################
# Generate a Cole-Cole signature

from pygimli.physics.SIP import ColeColeRho
import numpy as np
import pygimli as pg

f = np.logspace(-2, 5, 100)
m = np.linspace(0.1, 0.9, 5)
tau = 0.01

Z = ColeColeRho(f, rho=1, m=m[0], tau=tau, c=0.5)

fig, axes = pg.plt.subplots(2, 2, figsize=(15 / 2.54, 10 / 2.54))
ax = axes[0, 0]
ax.semilogx(f, np.abs(Z), '.-')
ax = axes[0, 1]
ax.semilogx(f, -np.angle(Z) * 1e3, '.-')
ax = axes[1, 0]
Y = 1 / Z
ax.semilogx(f, np.real(Y), '.-')
ax = axes[1, 1]
ax.semilogx(f, np.imag(Y), '.-')
fig.tight_layout()
pg.plt.show()
