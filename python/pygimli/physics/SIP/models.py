#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Spectral induced polarisation (SIP) relaxations models."""

from math import pi
import numpy as np
import pygimli as pg


def ColeColeRho(f, R, m, tau, c, a=1):
    r"""Frequency domain Cole-Cole impedance model after Pelton et al. (1978) :cite:`PeltonWarHal1978`

    .. math::

        Z(\omega) & = \rho_0\left[1 - m \left(1 - \frac{1}{1+(\text{i}\omega\tau)^c}\right)\right] \\
        \quad \text{with}\quad m & = \frac{1}{1+\frac{\rho_0}{\rho_1}} \quad \text{and}\quad \omega =2\pi f

    * :math:`Z(\omega)` - Complex impedance per 1A current injection
    * :math:`f` - Frequency
    * :math:`\rho_0` -- Background resistivity states the unblocked pore path
    * :math:`\rho_1` -- Resistance of the solution in the blocked pore passages
    * :math:`m` -- Chargeability proposed by Seigel (1959) :cite:`Seigel1959` as
      being the ratio of voltage immediately after, to the voltage immediately
      before cessation of an infinitely long charging current.
    * :math:`\tau` -- 'Time constant' [s] determines the length of time required
      for the decay in the time domain
    * :math:`c` - Frequency dependence, typical [0.1 .. 0.6]

    Examples
    --------
    >>> # no need to import matplotlib. pygimli's show does
    >>> import numpy as np
    >>> import pygimli as pg
    >>> from pygimli.physics.SIP import ColeColeRho
    >>> f = np.logspace(-2, 5, 100)
    >>> m = np.linspace(0.1, 0.9, 5)
    >>> fig, ax1 = pg.plt.subplots()
    >>> ax2 = ax1.twinx()
    >>> for i in range(len(m)):
    ...     Z = ColeColeRho(f, R=1, m=m[i], tau=0.01, c=0.5)
    ...     ax1.loglog(f, np.abs(Z), color='b')
    ...     print(np.angle(Z))
    ...     ax2.loglog(f, -np.angle(Z)*1000, color='g')
    >>> ax1.set_ylim(1e-3, 1)
    >>> ax1.set_ylabel('Amplitude $|Z(f)|$', color='b')
    >>> ax2.set_ylim(1, 1e3)
    >>> ax2.set_ylabel(r'- Phase $\varphi$ [mrad]', color='g')
    >>> ax1.set_xlabel('Frequency $f$ [Hz]')
    >>> pg.plt.show()
    """
    return (1. - m * (1. - relaxationTerm(f, tau, c, a))) * R


def relaxationTerm(f, tau, c=1., a=1.):
    """Auxiliary function for Debye type relaxation term."""
    return 1. / ((f * 2. * pi * tau * 1j)**c + 1.)**a


def DebyeRelaxation(f, tau, m):
    """Complex-valued single Debye relaxation term with chargeability."""
    return 1. - (1. - relaxationTerm(f, tau)) * m


def WarbugRelaxation(f, tau, m):
    """Complex-valued single Debye relaxation term with chargeability."""
    return 1. - (1. - relaxationTerm(f, tau, c=0.5)) * m


def ColeColeEpsilon(f, e0, eInf, tau, alpha):
    """Original complex-valued permittivity formulation (Cole&Cole, 1941)."""
    return (e0 - eInf) * relaxationTerm(f, tau, c=1./alpha) + eInf


def ColeColeSigma(f, R, m, tau, c, a=1):
    """Complex-valued conductivity Cole-Cole model"""
    return (1. + m / (1-m) * (1. - relaxationTerm(f, tau, c, a))) * R


def ColeCole(f, R, m, tau, c, a=1):
    """For backward compatibility.

    DEPRECATED?
    """
    return ColeColeRho(f, R, m, tau, c, a)


def ColeDavidson(f, R, m, tau, a=1):
    """For backward compatibility.

    DEPRECATED?
    """
    return ColeCole(f, R, m, tau, c=1, a=1)


class ColeColePhi(pg.ModellingBase):
    r"""Cole-Cole model with EM term after Pelton et al. (1978) :cite:`PeltonWarHal1978`

    Modelling operator for the Frequency Domain Cole-Cole impedance model
    :py:mod:`pygimli.physics.SIP.ColeColeRho`

    * :math:`\textbf{m} =\{ m, \tau, c\}`

        Modelling parameter for the Cole-Cole model with :math:`\rho_0 = 1`

    * :math:`\textbf{d} =\{\varphi_i(f_i)\}`

        Modeling Response for all given frequencies as negative the phase angles
        with :math:`\varphi(f) = -tan^{-1}\frac{\text{Im}\,Z(f)}{\text{Re}\,Z(f)}`
        and :math:`Z(f, \rho_0, m, \tau, c) =` ColeCole impedance.

    """
    def __init__(self, f, verbose=False):
        pg.ModellingBase.__init__(self, verbose)
        self.f_ = f
        self.setMesh(pg.createMesh1D(1, 3))

    def response(self, par):
        """Phase angle of the model."""
        spec = ColeCole(self.f_, 1.0, par[0], par[1], par[2])
        return -np.angle(spec)


class DoubleColeColePhi(pg.ModellingBase):
    r"""Double Cole-Cole model with EM term after Pelton et al. (1978) :cite:`PeltonWarHal1978`

    Modelling operator for the Frequency Domain Cole-Cole impedance model
    :py:mod:`Cole-Cole <pygimli.physics.SIP.ColeColeRho>`

    * :math:`\textbf{m} =\{ m_1, \tau_1, c_1, m_2, \tau_2, c_2\}`

        Modelling parameter for the Cole-Cole model with :math:`\rho_0 = 1`

    * :math:`\textbf{d} =\{\varphi_i(f_i)\}`

        Modeling Response for all given frequencies as negative the phase angles
        with :math:`\varphi(f) = \varphi_1(Z_1(f))+\varphi_2(Z_2(f)) = -tan^{-1}\frac{\text{Im}\,(Z_1(f)Z_2(f))}{\text{Re}\,(Z_1(f)Z_2(f))}`
        and :math:`Z_1(f, 1, m_1, \tau_1, c_1)` and
        :math:`Z_2(f, 1, m_2, \tau_2, c_2)` ColeCole impedances.

    """
    def __init__(self, f, verbose=False):  # initialize class
        pg.ModellingBase.__init__(self, verbose)  # call default constructor
        self.f_ = f                               # save frequencies
        self.setMesh(pg.createMesh1D(1, 6))       # 4 single parameters

    def response(self, par):
        """phase angle of the model"""
        spec1 = ColeCole(self.f_, 1.0, par[0], par[1], par[2])
        spec2 = ColeCole(self.f_, 1.0, par[3], par[4], par[5])
        return -np.angle(spec1 * spec2)
#        return -np.angle(spec1) - np.angle(spec2)


class ColeColeAbs(pg.ModellingBase):
    """Cole-Cole model with EM term after Pelton et al. (1978)"""
    def __init__(self, f, verbose=False):  # initialize class
        pg.ModellingBase.__init__(self, verbose)  # call default constructor
        self.f_ = f                               # save frequencies
        self.setMesh(pg.createMesh1D(1, 4))       # 3 single parameters

    def response(self, par):
        """phase angle of the model"""
        spec = ColeCole(self.f_, par[0], par[1], par[2], par[3])
        return np.abs(spec)


class ColeColeComplex(pg.ModellingBase):
    """Cole-Cole model with EM term after Pelton et al. (1978)"""
    def __init__(self, f, verbose=False):  # initialize class
        pg.ModellingBase.__init__(self, verbose)  # call default constructor
        self.f_ = f                               # save frequencies
        self.setMesh(pg.createMesh1D(1, 4))       # 4 single parameters

    def response(self, par):
        """phase angle of the model"""
        spec = ColeColeRho(self.f_, *par)
        return pg.cat(np.abs(spec), -np.angle(spec))


class ColeColeComplexSigma(pg.ModellingBase):
    """Cole-Cole model with EM term after Pelton et al. (1978)"""
    def __init__(self, f, verbose=False):  # initialize class
        pg.ModellingBase.__init__(self, verbose)  # call default constructor
        self.f_ = f                               # save frequencies
        self.setMesh(pg.createMesh1D(1, 4))       # 4 single parameters

    def response(self, par):
        """phase angle of the model"""
        spec = ColeColeSigma(self.f_, *par)
        return pg.cat(np.real(spec), np.imag(spec))


class PeltonPhiEM(pg.ModellingBase):
    """Cole-Cole model with EM term after Pelton et al. (1978)"""
    def __init__(self, f, verbose=False):  # initialize class
        pg.ModellingBase.__init__(self, verbose)  # call default constructor
        self.f_ = f                               # save frequencies
        self.setMesh(pg.createMesh1D(1, 4))       # 4 single parameters

    def response(self, par):
        """phase angle of the model"""
        spec = ColeCole(self.f_, 1.0, par[0], par[1], par[2]) * \
            relaxationTerm(self.f_, par[3])  # pure EM has c=1
        return -np.angle(spec)


class DebyePhi(pg.ModellingBase):
    """Debye decomposition (smooth Debye relaxations) phase only"""
    def __init__(self, fvec, tvec, verbose=False):  # save reference in class
        """constructor with frequecy and tau vector"""
        pg.ModellingBase.__init__(self, verbose)
        self.f_ = fvec
        self.nf_ = len(fvec)
        self.t_ = tvec
        self.mesh = pg.createMesh1D(len(tvec))  # standard 1d discretization
        self.setMesh(self.mesh)

    def response(self, par):
        """amplitude/phase spectra as function of spectral chargeabilities"""
        y = np.ones(self.nf_, dtype=np.complex)  # 1 -
        for (tau, mk) in zip(self.t_, par):
            y -= (1. - relaxationTerm(self.f_, tau)) * mk

        return -np.angle(y)


class DebyeComplex(pg.ModellingBase):
    """Debye decomposition (smooth Debye relaxations) of complex data"""
    def __init__(self, fvec, tvec, verbose=False):  # save reference in class
        """constructor with frequecy and tau vector"""
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
        """amplitude/phase spectra as function of spectral chargeabilities"""
        return self.J * par

    def createJacobian(self, par):
        """linear jacobian after Nordsiek&Weller (2008)"""
        pass


if __name__ == "__main__":
    pass
