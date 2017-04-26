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
    def __init__(self, **kwargs):
        """Initialize forward operator.

        Initialize forward operator with optional halfspace geometry
        and measurement spezifications.

        Parameters
        ----------

        **kwargs : dict()
            Forward to ModellingBase
        """
        pg.ModellingBase.__init__(self, **kwargs)


    def response_mt(self, par, i=0):
        """Compute response vector for a set of model parameter.

        Parameters
        ----------
        par : iterabale
            DOCUMENTME
        """
        THROW_TO_IMPL

    def response(self, par, i=0):
        return self.response_mt(par)

    def calcEPhiF(self, f, rho, d, rmin=1, nr=41, ze=0, zs=0, tm=1):
        """VMD E-phi field."""
        #ze    z-Koordinate des Empfaengers in Meter (default: ZE=0)
        #zs    z-Koordinate des Senders in Meter (default: ZS=0)

        # Filterkoeffizienten

        q = 10**0.1
        rr = rmin * q **(np.array(nr) - 1)

        he = ze
        hs = zs
        zm = hs - he
        zp = hs + he

        rm = np.sqrt(rr**2 + zm**2)
        rp = np.sqrt(rr**2 + zp**2)

        fcJ1, nc0 = pg.utils.hankelFC(4)
        nc = len(fcJ1)

        ncnr = nc + nr - 1

        # Vakuumwerte, nicht normiert
        if ze <= 0:
            ePhi = rr / rp**3

        # Create Wavenumbers
        nu = np.arange(1, ncnr + 1)

        n = nc0 - nc + nu
        q = np.log(10) * 0.1
        k = np.exp(-(n-1) * q) / rmin

        if ze <= 0:
            ePhi = rr / rp **3.0

        # Admittanzen for halfspace borders for each Wavenumbers k
        bt = np.zeros(ncnr) * 1j
        aa = np.zeros(ncnr) * 1j
        for i in range(ncnr):
            if ze <= 0:
                bt[i] = self.btp(k[i], f, rho, d, type=1)
            else:
                raise Exception('NeedTests')
                aa[i], aap[i], bt[i] = downward(k[i], f, rho, d, ze)

        # Kernel functions
        if ze <= 0:
            e = np.exp(k * ze) * np.exp(k * zs)
            delta = e * (bt - k) / (bt + k)
        else:
            raise Exception('NeedTests')
            e = np.exp(k * zs)
            delta = 2. * k**2. / (bt + k)*e

        # convolution
        aux3 = np.zeros(nr) * 1j

        for n in range(nr):
            for nn in range(nc):
                nu = nn + n;
                mn = nc - nn;

                nnn = nc0 - nc + nu;
                k = np.exp( - (nnn) * q) / rmin

                del0 = delta[nu]

                if ze <= 0:
                    del1 = del0 * k
                    aux3[n] = aux3[n] + del1 * fcJ1[mn-1]
                else:
                    raise Exception('NeedTests')
                    del3=del0 * aa(nu);
                    aux3[n] = aux3[n] + del3 * fcJ1[mn]

            if nr > 1:
                aux3[n] = aux3[n] / rr[n]
                raise Exception("more then one r .. check code here")
            else:
                aux3[n] = aux3[n] / rr

        if ze <= 0:
            # Air
            ePhi = ePhi - aux3
        else:
            raise Exception("more then one r .. check code here")
            # Halfspace
            ePhi = aux3

        # Normalization
        ePhi *= -tm * 1j * (2*pi*f) * pg.physics.constants.mu0 / (4*pi)
        return ePhi

    def btp(self, k, f, rho, d, type=1):
        """Admittance of a layered halfspace
        for TE (TYP=1) and TM (TYP=2) mode.
        k: wave number (1/m)
        f: frequency (1/s)
        rho: layer resitivities
        d: layer thicknesses
        """
        nl = len(rho)
        mu0 = 4e-7 * pi
        c = 1j * pg.physics.constants.mu0 * 2 * pi * f
        b = np.sqrt(k**2 + c/rho[nl-1])

        if nl > 1:
            beta = 1
            for nn in range(nl-2, -1, -1):
                alpha = np.sqrt(k**2 + c / rho[nn])
                if type == 2:
                    beta = rho[nn] / rho[nn + 1]

                cth = np.exp(-2. * d[nn] * alpha)
                cth = (1-cth) / (1+cth)
                b = (b + alpha * beta * cth) / (beta + cth * b / alpha)

        return b


class VMDTimeDomainModelling(VMDModelling):
    """
    """
    def __init__(self, t, txarea, rxarea=None, **kwargs):
        """
        """
        VMDModelling.__init__(self, **kwargs)

        self.t = t
        self.txarea = txarea

        if rxarea is None:
            self.rxarea = txarea
        else:
            self.rxarea = rxarea

    def response_mt(self, par, i=0):
        """
        par = [thicknesses, res]
        """
        nLay = (len(par)-1)/2
        thk = par[0:nLay]
        res = par[nLay:]

        return self.calcRhoa(thk, res)

    def response(self, par):
        """par = [thicknesses(nLay), res(nlay + 1)]"""
        return self.response_mt(par, 0)

    def calcRhoa(self, thk, res):
        """
        """
        a = sqrt(self.txarea / pi) # TX coil radius

        ePhiTD, tD = self.calcEphiT(tMin=min(self.t), tMax=max(self.t),
                                    rho=res, d=thk, rMin=a, rMax=a, z=0,
                                    dipm=self.rxarea)

        ePhi = np.exp(np.interp(np.log(self.t),
                                np.log(tD), np.log(ePhiTD[:,0])))

        tmp = a**(4./3) * self.rxarea**(2./3) * \
            pg.physics.constants.mu0**(5./3)/ (20**(2./3) * pi**(1./3))
        rhoa = tmp / (self.t**(5./3) * (ePhi * 2 * pi * a)**(2./3))

        return rhoa


    def calcEphiT(self, tMin, tMax, rho, d, rMin, rMax, z, dipm):
        """
        """
        nl = len(rho)

        r = pg.utils.niceLogspace(rMin, rMax, nDec=10)
        r = [rMin]

        t = pg.utils.niceLogspace(tMin, tMax, nDec=10)

        fcS, nc0 = pg.utils.hankelFC(1)

        nr = len(r)
        nc = len(fcS)
        nt = len(t)
        ncnt = nc + nt

        fef = np.zeros((ncnt,nr), np.complex)
        ePhi = np.zeros((nt,nr))

        omega = 10. ** (0.1 * (1 - (-nc + nc0 + np.arange(1, ncnt)))) / t[0]

        for nn, o in enumerate(omega):
            f = o / (2. * pi)
            ePhiF = self.calcEPhiF(f, rho, d, ze=z, zs=0.,
                                   rmin=r[0], nr=len(r), tm=1.0)

            fef[nn,:] = ePhiF / sqrt(o)

            ita = max(0, nn-nc)
            ite = min(nt, nn+1)
            for it in range(ita, ite):
                itn = nc-nn+it-1
                ePhi[it] = ePhi[it] - np.real(fef[nn]) * fcS[itn]

        for it in range(nt):
            ePhi[it] = ePhi[it] * dipm * sqrt(2/pi / t[it])

        return ePhi, t
