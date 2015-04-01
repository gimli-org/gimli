#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Spectral induced polarisation (SIP) relaxations models
"""

from math import pi
import numpy as np
from scipy.integrate import simps
import pygimli as pg


# basic functions used by the modelling operators
def relaxationTerm(f, tau, c=1., a=1.):
    """ auxiliary function for Debye type relaxation term """
    return 1. / ((f * 2. * pi * tau * 1j)**c + 1)**a


def DebyeRelaxation(f, tau, m):
    """ complex-valued single Debye relaxation term """
    return 1. - (1.-relaxationTerm(f, tau))*m


def ColeColeRho(f, R, m, tau, c, a=1):
    """ Complex valued Cole-Cole model """
    return (1. - m * (1. - relaxationTerm(f, tau, c, a))) * R


def ColeColeSigma(f, R, m, tau, c, a=1):
    """ Complex valued Cole-Cole model """
    return (1. + m / (1-m) * (1. - relaxationTerm(f, tau, c, a))) * R


def ColeCole(f, R, m, tau, c, a=1):
    """ for backward compatibility """
    return ColeColeRho(f, R, m, tau, c, a)


# modelling operators for use with pygimli inversion
class ColeColePhi(pg.ModellingBase):
    """" Cole-Cole model with EM term after Pelton et al. (1978)"""
    def __init__(self, f, verbose=False):  # initialize class
        pg.ModellingBase.__init__(self, verbose)  # call default constructor
        self.f_ = f                               # save frequencies
        self.setMesh(pg.createMesh1D(1, 3))       # 4 single parameters

    def response(self, par):
        """ phase angle of the model """
        spec = ColeCole(self.f_, 1.0, par[0], par[1], par[2])
        return -np.angle(spec)


class ColeColeComplex(pg.ModellingBase):
    """" Cole-Cole model with EM term after Pelton et al. (1978)"""
    def __init__(self, f, verbose=False):  # initialize class
        pg.ModellingBase.__init__(self, verbose)  # call default constructor
        self.f_ = f                               # save frequencies
        self.setMesh(pg.createMesh1D(1, 4))       # 4 single parameters

    def response(self, par):
        """ phase angle of the model """
        spec = ColeColeRho(self.f_, *par)
        return pg.cat(np.abs(spec), -np.angle(spec))


class ColeColeComplexSigma(pg.ModellingBase):
    """" Cole-Cole model with EM term after Pelton et al. (1978)"""
    def __init__(self, f, verbose=False):  # initialize class
        pg.ModellingBase.__init__(self, verbose)  # call default constructor
        self.f_ = f                               # save frequencies
        self.setMesh(pg.createMesh1D(1, 4))       # 4 single parameters

    def response(self, par):
        """ phase angle of the model """
        spec = ColeColeSigma(self.f_, *par)
        return pg.cat(np.real(spec), np.imag(spec))


class PeltonPhiEM(pg.ModellingBase):
    """" Cole-Cole model with EM term after Pelton et al. (1978)"""
    def __init__(self, f, verbose=False):  # initialize class
        pg.ModellingBase.__init__(self, verbose)  # call default constructor
        self.f_ = f                               # save frequencies
        self.setMesh(pg.createMesh1D(1, 4))       # 4 single parameters

    def response(self, par):
        """ phase angle of the model """
        spec = ColeCole(self.f_, 1.0, par[0], par[1], par[2]) * \
            relaxationTerm(self.f_, par[3])  # pure EM has c=1
        return -np.angle(spec)


class DebyePhi(pg.ModellingBase):
    """ Debye decomposition (smooth Debye relaxations) phase only """
    def __init__(self, fvec, tvec, verbose=False):  # save reference in class
        """ constructor with frequecy and tau vector """
        self.f_ = fvec
        self.nf_ = len(fvec)
        self.t_ = tvec
        mesh = pg.createMesh1D(len(tvec))  # standard 1d discretization
        pg.ModellingBase.__init__(self, mesh, verbose)

    def response(self, par):
        """ amplitude/phase spectra as function of spectral chargeabilities """
        y = np.ones(self.nf_, dtype=np.complex)  # 1 -
        for (tau, mk) in zip(self.t_, par):
            y -= (1. - relaxationTerm(self.f_, tau)) * mk

        return -np.angle(y)


class DebyeComplex(pg.ModellingBase):
    """ Debye decomposition (smooth Debye relaxations) of complex data """
    def __init__(self, fvec, tvec, verbose=False):  # save reference in class
        """ constructor with frequecy and tau vector """
        self.f = fvec
        self.nf = len(fvec)
        self.t = tvec
        self.nt = len(tvec)
        mesh = pg.createMesh1D(len(tvec))  # standard 1d discretization
        pg.ModellingBase.__init__(self, mesh, verbose)
        T, W = np.meshgrid(tvec, fvec * 2. * pi)
        WT = W*T
        self.A = WT**2 / (WT**2 + 1)
        self.B = WT / (WT**2+1)
        self.J = pg.RMatrix()
        self.J.resize(len(fvec)*2, len(tvec))
        for i in range(self.nf):
            wt = fvec[i] * 2.0 * pi * tvec
            wt2 = wt**2
            self.J[i] = wt2 / (wt2 + 1.0)
            self.J[i+self.nf] = wt / (wt2 + 1.0)

        self.setJacobian(self.J)

    def response(self, par):
        """ amplitude/phase spectra as function of spectral chargeabilities """
        return self.J * par

    def createJacobian(self, par):
        """ linear jacobian after Nordsiek&Weller (2008) """
        pass


if __name__ == "__main__":
    pass
