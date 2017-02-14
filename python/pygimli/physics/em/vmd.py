#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Electromagnetics for Vertical Magnetic Dipole (VMD) functions and class."""

from math import sqrt, pi

import numpy as np

import pygimli as pg

from pygimli.physics import constants


class VMDModelling(pg.ModellingBase):
    """Modelling operator for a Vertical Magnetic Dipole (VMD).

    Modelling operator for a Vertical Magnetic Dipole (VMD) to calculate the
    electromagnetic response in cylindric coordinates
    ::math::`H_z, H_r, E_{\phi} \in I\!C`
    for a layered halfspace (1D) using Hankel transformation.

    VMD is at the origin ::math::`r_s = (0.0)` and ::math::`z_s=0`.

    """
    def __init__(self, f=None, **kwargs):
        """Initialize forward operator.

        Initialize forward operator with optional halfspace geometry
        and measurement spezifications.

        Parameters
        ----------
        layers : iterable [float]
            Array of layer thinnumber of layers

        f : iterable [float]
            Frequency vector

        **kwargs : dict()
            Forward to ModellingBase
        """
        pg.ModellingBase.__init__(self, **kwargs)

        self.nlay = 0     # Number of halfspace layers
        self.wem = 0      # ??
        self.iwm = 0      # ??

        self.setFrequencies(f)

    def setFrequencies(self, f):
        """Set measurement frequencies.

        Parameters
        ----------
        f : iterable [float]
            Frequency vector, Default is [1000] i.e. 1kHz
        """
        if f is None:
            self.f = np.asarray([1000])
        else:
            self.f = np.asarray(f)

        self.wem = (2.0 * pi * self.f) ** 2 * constants.e0 * constants.mu0
        self.iwm = 2.0 * pi * self.f * constants.mu0 * 1j


    def response(self, par):
        """Compute response vector for a set of model parameter.

        Parameters
        ----------
        par : iterabale
            DOCUMENTME
        """
        THROW_TO_IMPL

    @staticmethod
    def fdResponse(f, rho, d, ze, zs, rmin, nr, tm):
        """Calculate frequency domain response for a given model.

        Parameters
        ----------

        f :

        depths : iterable [float]
            Array of layer depths in positive z coordinates.

        rho : iterable [float]
            Array of electrical resistivities for the given layers.
            Must be length of depths plus one for the lower background

        Returns
        -------
        hz : np.array [complex]
            Vertical component of the magnetic field ::math::`H_z(f)` in [A/m]
        hr : np.array [complex]
            Radial component of the magnetic field ::math::`H_r(f)` in [A/m]
        ephi : np.array [complex]
            Tangential component ::math::`E_{\phi}(f)` of the electric field in [V/m]

        """
        epr = 1.  # epsilonR
        mur = 1.  # muR
        quasistatic = False
        h = 0.01
        rho = np.array(rho)
        #thickness = np.array(pg.abs(pg.utils.diff(np.array([0.] +  depths))))
        self.r = np.array([1.])

        # relative electrical permeabilty for each layer
        if isinstance(epr, float):
            epr = np.ones(len(rho)) * epr
        # relative magnetic permeabilty for each layer
        if isinstance(mur, float):
            mur = np.ones(len(rho)) * mur

        fc0, nc0 = pg.utils.hankelFC(3)  # J0
        fc1, nc0 = pg.utils.hankelFC(4)  # J1
        nc = len(fc0)

        # arrays anlegen
        nf = len(f)
        lam = np.zeros((1, nc, nf))
        alpha0 = np.zeros((1, nc, nf), np.complex)
        delta0 = np.zeros((1, nc, nf), np.complex)
        delta1 = np.zeros((1, nc, nf), np.complex)
        delta2 = np.zeros((1, nc, nf), np.complex)
        delta3 = np.zeros((1, nc, nf), np.complex)
        aux0 = np.zeros((1, nf), np.complex)
        aux1 = np.zeros((1, nf), np.complex)
        aux2 = np.zeros((1, nf), np.complex)
        aux3 = np.zeros((1, nf), np.complex)
        Z = np.zeros(nf, np.complex)

        r0 = np.copy(self.r)

        # determine optimum r0 (shift nodes) for f > 1e4 and h > 100
        if quasistatic:
            index = np.zeros(self.f.shape, np.bool)
        else:
            index = np.logical_and(self.f >= 1e4, h >= 100.0)

        if np.any(index):
            opt = np.floor(10.0 * np.log10(self.r[index] * 2.0 * pi *
                           self.f[index] / constants.c) + nc0)
            r0[index] = constants.c / (2.0 * pi * self.f[index]) * \
                        10.0 ** ((opt + 0.5 - nc0) / 10.0)

        # prepare wavenumbers
        n = np.arange(nc0 - nc, nc0, 1, np.float)
        q = 0.1 * np.log(10)
        lam = np.reshape(np.exp(-n[np.newaxis, :, np.newaxis] * q) /
                         r0[np.newaxis, np.newaxis, :], (-1, nc, r0.size))

        # (1, 100, nfreq)
        # wavenumbers in air, quasistatic approximation
        alpha0 = lam * np.complex(1, 0)  # (1, 100, nfreq)

        # wavenumbers in air, complete solution (f > 1e4 Hz)
        if quasistatic:
            index = np.zeros(self.f.shape, np.bool)
        else:
            index = self.f >= 1e4

        if np.any(index):
            alpha0[:, :, index] = np.sqrt(lam[:, :, index]**2 - \
                np.tile(self.wem[index], (1, nc, 1)) +
                np.tile(self.iwm[index], (1, nc, 1)) / 1e9)

        # (1, 100 , nfreq)
        # Admittanzen on the surface of the layered Halfspace
        b1, aa, aap = self.downward(rho, thickness, 0.0, epr, mur, lam)

        # Kernel functions
        e = np.exp(-2.0 * h * alpha0)  # (1, 100, nfreq)
        delta0 = (b1 - alpha0 * mur[0]) / (b1 + alpha0 * mur[0]) * e
        delta1 = (2. * mur[0]) / (b1 + alpha0 * mur[0]) * e  # (1, 100, nfreq)

        delta2 = 1. / h * e  # (1, 100, nfreq)
        delta3 = 1. / (2. * h) * e  # (1, 100, nfreq)

        # Convolution
        # quasistatic approximation
        aux0 = np.sum(delta0 * lam ** 3 / alpha0 *
                      np.tile(fc0[::-1].T[:, :, np.newaxis],
                              (1, 1, self.f.size)), 1, np.complex) / r0

        Z = self.r ** 3 * aux0 * self.scaling

        # full solution, partial integration
        if np.any(index):
            aux1 = np.sum(delta1 * lam ** 3 *
                          np.tile(fc0[::-1].T[:, :, np.newaxis],
                                  (1, 1, self.f.size)), 1, np.complex) / r0

            aux2 = np.sum(delta2 * lam *
                          np.tile(fc0[::-1].T[:, :, np.newaxis],
                                  (1, 1, self.f.size)), 1, np.complex) / r0

            aux3 = np.sum(delta3 * lam ** 2 *
                          np.tile(fc1[::-1].T[:, :, np.newaxis],
                                  (1, 1, self.f.size)), 1, np.complex) / r0

            Z[:, index] = (-self.r[index]**3 * aux1[:, index] +
                            self.r[index]**3 * aux2[:, index] -
                            self.r[index]**4 * aux3[:, index]) * self.scaling

        return np.real(Z[0]), np.imag(Z[0])


    def downward(self, rho, d, z, epr, mur, lam):
        """Downward continuation of fields.

        For internal use .. DOCUMENTME

        Parameters
        ----------

        Returns
        -------

        AUTHORS: RUB, AU, TG ?

        """
        # Anzahl der Schichten
        nl = len(rho)

        # arrays anlegen
        alpha = np.zeros((nl, lam.shape[1], self.f.size), np.complex)
        b = np.zeros((nl, lam.shape[1], self.f.size), np.complex)
        aa = np.zeros((nl, lam.shape[1], self.f.size), np.complex)
        aap = np.zeros((nl, lam.shape[1], self.f.size), np.complex)

        rho = rho[:, np.newaxis, np.newaxis] * np.ones(
            (rho.size, lam.shape[1], self.f.size), np.float)

        d = d[:, np.newaxis, np.newaxis] * np.ones(
            (d.size, lam.shape[1], self.f.size), np.float)
        h = np.insert(np.cumsum(d[:, 0, 0]), 0, 0)
        epr = epr[:, np.newaxis, np.newaxis] * np.ones(
            (epr.size, lam.shape[1], self.f.size), np.float)
        mur = mur[:, np.newaxis, np.newaxis] * np.ones(
            (mur.size, lam.shape[1], self.f.size), np.float)
        lam = np.tile(lam, (nl, 1, 1))
        # Ausbreitungskonstante
        alpha = np.sqrt(lam ** 2 - np.tile(self.wem, (nl, lam.shape[1], 1)) *
                        epr * mur + np.tile(self.iwm, (nl, lam.shape[1], 1)) *
                        mur / rho)
        if nl == 1:  # homogenous halfspace
            b1 = alpha.copy()  # (1, 100, nfreq)
            a = np.exp(-alpha * z)
            ap = a.copy()
            return b1, a, ap
        elif nl > 1:  # multi-layer case
            # tanh num instabil tanh(x)=(exp(x)-exp(-x))/(exp(x)+exp(-x))
            ealphad = np.exp(-2.0 * alpha[0:-1, :, :] * d)
            talphad = (1.0 - ealphad) / (1.0 + ealphad)
            b[-1, :, :] = np.copy(alpha[-1, :, :])
            # rekursive Berechnung der Admittanzen an der Oberkante der Schicht
            # von unten nach oben, für nl-1 Schichten
            for n in range(nl-2, -1, -1):
                b[n, :, :] = alpha[n, :, :] * (b[n+1, :, :] + alpha[n, :, :] *
                                               talphad[n, :, :]) / (
                    alpha[n, :, :] + b[n+1, :, :] * talphad[n, :, :])
            # Impedance
            c = 1.0 / b
            # b1 == 1. row in b (nl, 100, nfreq)
            b1 = np.copy(b[0, :, :][np.newaxis, :, :])  # (1, 100, nfreq)
            # Variation from one layer boundary to the other
            for n in range(0, nl-1):
                aa[n, :, :] = (b[n, :, :] + alpha[n, :, :]) / (
                    b[n+1, :, :] + alpha[n, :, :]) * \
                    np.exp(-alpha[n, :, :] * d[n, :, :])
                aap[n, :, :] = (1.0 + alpha[n, :, :] * c[n, :, :]) / (
                    1.0 + alpha[n, :, :] * c[n+1, :, :]) * \
                    np.exp(-alpha[n, :, :] * d[n, :, :])
            # Determin layer Index where z is

            for n in range(0, nl-1):
                if np.logical_and(z >= h[n], z < h[n+1]):
                    ind = n
            try:
                ind
            except NameError:
                ind = nl - 1

            print("check ind", ind, nl-1)


            if (ind + 1) < nl:
                a = np.prod(aa[:ind, :, :], 0) * 0.5 * (
                    1.0 + b[ind, :, :] / alpha[ind, :, :]) * \
                    (np.exp(-alpha[ind, :, :] * (z - h[ind])) -
                        (b[ind+1, :, :] - alpha[ind, :, :]) /
                        (b[ind+1, :, :] + alpha[ind, :, :]) *
                        np.exp(-alpha[ind, :, :] *
                               (d[ind, :, :] + h[ind+1] - z)))
                ap = np.prod(aap[:ind, :, :], 0) * 0.5 * \
                    (1.0 + alpha[ind, :, :] * c[ind, :, :]) * \
                    (np.exp(-alpha[ind, :, :] * (z - h[ind])) +
                        (1.0 - alpha[ind, :, :] * c[ind+1, :, :]) /
                        (1.0 + alpha[ind, :, :] * c[ind+1, :, :]) *
                        np.exp(-alpha[ind, :, :] *
                               (d[ind, :, :] + h[ind+1] - z)))
            else:
                a = np.prod(aa, 0) * np.exp(-alpha[ind, :, :] * (z - h[ind]))
                ap = np.prod(aap, 0) * np.exp(-alpha[ind, :, :] * (z - h[ind]))
            a = a[np.newaxis, :, :]  # (1, 100, nfreq)
            ap = ap[np.newaxis, :, :]  # (1, 100, nfreq)
            return b1, a, ap

