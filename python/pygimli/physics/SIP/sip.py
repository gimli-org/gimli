import numpy as np
import pygimli as pg
import matplotlib.pyplot as plt
from math import pi  # , log10, exp


def showAmplitudeSpectrum(ax, freq, amp, ylabel=r'$\rho_a$ in $\Omega$m',
                          grid=True, **kwargs):
    """ show amplitude spectrum """
    ax.loglog(freq, amp, **kwargs)
    ax.set_xlabel('f in Hz')
    ax.set_ylabel(ylabel)
    ax.set_grid(grid)


def showPhaseSpectrum(ax, freq, phi, ylabel=r'$\phi_a$ in mrad',
                      grid=True, **kwargs):
    """ show amplitude spectrum """
    ax.loglog(freq, phi, **kwargs)
    ax.set_xlabel('f in Hz')
    ax.set_ylabel(ylabel)
    ax.set_grid(grid)


def showSpectrum(freq, amp, phi):
    fig, ax = plt.subplots(nrows=2, sharex=True)
    showAmplitudeSpectrum(ax[0], freq, amp)
    showPhaseSpectrum(ax[1], freq, phi)


def relaxationTerm(f, tau, c=1.):
    ''' relaxation term as used by the Cole-Cole model (and the EM term) '''
    return 1. / ((f * 2. * pi * tau * 1j)**c + 1)


def ColeCole(f, R, m, tau, c):
    ''' Complex valued Cole-Cole model without R '''
    return (1. - m * (1. - relaxationTerm(f, tau, c))) * R


class PeltonPhiEM(pg.ModellingBase):
    ''' modelling base that returns amplitude and phase '''
    def __init__(self, f, verbose=False):  # initialize class
        pg.ModellingBase.__init__(self, verbose)  # call default constructor
        self.f_ = f                               # save frequencies
        self.setMesh(pg.createMesh1D(1, 4))       # 4 single parameters

    def response(self, par):
        ''' Cole-Cole model with EM term after Pelton et al. (1978) '''
        spec = ColeCole(self.f_, 1.0, par[0], par[1], par[2]) * \
            relaxationTerm(self.f_, par[3])  # pure EM has c=1
        return -np.angle(spec)


class DebyePhi(pg.ModellingBase):
    """forward operator for Debye decomposition."""
    def __init__(self, fvec, tvec, verbose=False):  # save reference in class
        self.f_ = fvec
        self.nf_ = len(fvec)
        self.t_ = tvec
        mesh = pg.createMesh1D(len(tvec))  # standard 1d discretization
        pg.ModellingBase.__init__(self, mesh, verbose)

    def response(self, par):
        """amplitude/phase spectrum as function of spectral chargeabilities."""
        y = np.ones(self.nf_, dtype=np.complex)  # 1 -
        for (tau, mk) in zip(self.t_, par):
            y -= (1. - relaxationTerm(self.f_, tau)) * mk

        return -np.angle(y)


if __name__ == "__main__":
    pass
