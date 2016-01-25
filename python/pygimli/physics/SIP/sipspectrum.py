#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Spectral induced polarisation (SIP) module
"""

from math import log10, exp, pi
import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
from .importexport import readTXTSpectrum
from .plotting import showAmplitudeSpectrum, showSpectrum
from .models import DebyePhi, DebyeComplex, relaxationTerm
from .tools import KramersKronig, fitCCEMPhi, fitCCC, fitCCCC, fitCCPhi


class SIPSpectrum():
    """SIP spectrum data analysis"""
    def __init__(self, filename=None, unify=False, onlydown=True,
                 f=None, amp=None, phi=None, k=1, basename='new'):
        """init SIP class with either filename to read or data vectors

        Examples
        --------
        >>> #sip = SIPSpectrum('sipexample.txt', unify=True) # unique f values
        >>> #sip = SIPSpectrum(f=f, amp=R, phi=phase, basename='new')
        """
        self.basename = basename
        if filename is not None:
            if filename.rfind('.txt') > 0:
                self.basename = filename[:-4]
                self.f, self.amp, self.phi = readTXTSpectrum(filename)
                self.amp *= k
        if f is not None:
            self.f = np.asarray(f)
        if amp is not None:
            self.amp = np.asarray(amp)
        if phi is not None:
            self.phi = np.asarray(phi)
        if unify:
            self.unifyData(onlydown)

    def __repr__(self):
        """String representation of the class."""
        return self.__str__()

    def __str__(self):
        """Human readable string representation of the class."""
        out = self.__class__.__name__ + " object"
        if hasattr(self, 'f'):
            if hasattr(self.f, '__iter__'):
                out += "\nnf=" + str(len(self.f)) + " min/max="
                out += str(min(self.f)) + "/" + str(max(self.f))
        return out

    def unifyData(self, onlydown=False):
        """ unify data (only one value per frequency) by mean or selection"""
        fu = np.unique(self.f)
        if len(fu) < len(self.f) or onlydown:
            if onlydown:
                wende = min(np.nonzero(np.diff(self.f) > 0)[0])
                if wende > 0:
                    self.f = self.f[wende::-1]
                    self.amp = self.amp[wende::-1]
                    self.phi = self.phi[wende::-1]
            else:
                amp = np.zeros(fu.shape)
                phi = np.zeros(fu.shape)
                for i in range(len(fu)):
                    ind = np.nonzero(self.f == fu[i])[0]
                    amp[i] = np.mean(self.amp[ind])
                    phi[i] = np.mean(self.phi[ind])

                self.f = fu
                self.amp = amp
                self.phi = phi

    def sortData(self):
        """sort data along increasing frequency (e.g. useful for KK)"""
        ind = np.argsort(self.f)
        self.amp = self.amp[ind]
        self.phi = self.phi[ind]
        self.f = self.f[ind]

    def cutF(self, fcut=1e99):
        """ cut frequencies above a certain value """
        self.amp = self.amp[self.f <= fcut]
        self.phi = self.phi[self.f <= fcut]
        if hasattr(self, 'phiOrg'):
            self.phiOrg = self.phiOrg[self.f <= fcut]
        # finally cut f
        self.f = self.f[self.f <= fcut]

    def realimag(self, cond=False):
        """real and imaginary part"""
        if cond:
            amp = 1. / self.amp
        else:
            amp = self.amp
        return amp * np.cos(self.phi), amp * np.sin(self.phi)

    def zNorm(self):
        """normalized real (difference) and imag. z (Nordsiek&Weller, 2008)"""
        re, im = self.realimag()
        R0 = max(self.amp)
        zNormRe = 1. - re / R0
        zNormIm = im / R0
        return zNormRe, zNormIm

    def showData(self, reim=False, znorm=False, nrows=2):
        """Show amplitude and phase spectrum in two subplots

        Parameters
        ----------
        reim : bool
            show real/imaginary part instead of amplitude/phase

        znorm : bool (true forces reim)
            use normalized real/imag parts after Nordsiek&Weller (2008)

        nrows - use nrows subplots (default=2)

        Returns
        -------
        fig, ax : mpl.figure, mpl.axes array
        """
        if reim or znorm:
            addstr = ''
            if znorm:
                re, im = self.zNorm()
                addstr = ' (norm)'
            else:
                re, im = self.realimag()

            fig, ax = showSpectrum(self.f, re, im, ylog=False, nrows=nrows)
            ax[0].set_ylabel('real part'+addstr)
            ax[1].set_ylabel('imaginary part'+addstr)
        else:
            fig, ax = showSpectrum(self.f, self.amp, self.phi)

        plt.show(block=False)
        return fig, ax

    def getKK(self, use0=False):
        """retrieve Kramers-Kronig values (re->im and im->re)"""
        re, im = self.realimag()
        if False:
            ind = np.argsort(self.f)
            reKK, imKK = KramersKronig(self.f[ind], re[ind], im[ind],
                                       usezero=use0)
        else:
            fsort, ind = np.unique(self.f, return_index=True)

        reKK, imKK = KramersKronig(fsort, re[ind], im[ind], usezero=use0)
        re[ind] = reKK  # sort back
        im[ind] = imKK
        return re, im

    def getPhiKK(self, use0=False):
        """ get phase from Kramers-Kronig relations """
        re, im = self.realimag()
        reKK, imKK = self.getKK(use0)
        return np.arctan2(imKK, re)

    def showDataKK(self, use0=False):
        """show data as real/imag subplots along with Kramers-Kronig curves"""
        fig, ax = self.showData(reim=True)
        reKK, imKK = self.getKK(use0)
        ax[0].plot(self.f, reKK, label='KK')
        ax[1].plot(self.f, imKK, label='KK')
        for i in (0, 1):
            ax[i].set_yscale('linear')
            ax[i].legend()

    def checkCRKK(self, useEps=False, use0=False, ax=None):
        """check coupling removal (CR) by Kramers-Kronig (KK) relation"""
        if ax is None:
            fig, ax = plt.subplots()
        ax.semilogx(self.f, self.phi*1000, label='org')
        ax.semilogx(self.f, self.getPhiKK(use0)*1000, label='orgKK')
        if useEps:
            self.removeEpsilonEffect()
        else:
            self.fitCCEM()
        ax.semilogx(self.f, self.phi*1000, label='corr')
        ax.semilogx(self.f, self.getPhiKK(use0)*1000, label='corrKK')
        ax.grid(True)
        ax.legend(loc='best')

    def showPolarPlot(self):
        """show data in a polar plot (real against imaginary parts)"""
        re, im = self.realimag()
        fig, ax = plt.subplots()
        ax.plot(re, im, 'b.')
        ax.set_aspect(1)
        ax.grid(True)
        for i in range(0, len(re), 5):
            fi = self.f[i]
            mul = 10**np.floor(np.log10(fi) - 2)  # 3 counting digits
            ax.text(re[i], im[i], str(int(fi / mul) * mul))

        return fig, ax

    def epsilonR(self):
        """calculate relative permittivity from imaginary conductivity"""
        ECr, ECi = self.realimag(cond=True)
        eps0 = 8.854e-12
        we0 = self.f*2*pi*eps0  # Omega epsilon_0
        return ECi/we0

    def removeEpsilonEffect(self, er=None, mode=0):
        """remove effect of (constant high-frequency) epsilon from sigma

        Parameters
        ----------
        er : float
            relative epsilon to correct for (else automatically determined)

        Returns
        -------
        er : float
            determined er
        """
        ECr, ECi = self.realimag(cond=True)
        we0 = self.f*2*pi*8.854e-12  # Omega epsilon_0
        if er is None:  #
            epsr = ECi / we0
            nmax = np.argmax(self.f)
            if mode == 0:
                f1 = self.f * 1
                f1[nmax] = 0
                nmax1 = np.argmax(f1)
                er = 2*epsr[nmax] - epsr[nmax1]
            else:
                er = epsr[nmax]
            print("detected epsilon of ", er)

        ECi -= er * we0
        self.phiOrg = self.phi
        self.phi = np.arctan(ECi/ECr)
        self.ampOrg = self.amp
        self.amp = 1. / np.sqrt(ECr**2 + ECi**2)
        return er

    def fitCCPhi(self, ePhi=0.001, lam=1000., mpar=(0.2, 0, 1),
                 taupar=(1e-2, 1e-5, 100), cpar=(0.3, 0, 1)):
        """fit a Cole-Cole term to phase only

        Parameters
        ----------
        ePhi : float
            absolute error of phase angle
        lam : float
            regularization parameter
        mpar, taupar, cpar : list[3]
            inversion parameters (starting value, lower bound, upper bound)
            for Cole-Cole parameters (m, tau, c) and EM relaxation time (em)

        """
        self.mCC, self.phiCC = fitCCPhi(self.f, self.phi, ePhi, lam, mpar=mpar,
                                        taupar=taupar, cpar=cpar)

    def fitCCEM(self, ePhi=0.001, lam=1000., remove=True,
                mpar=(0.2, 0, 1), taupar=(1e-2, 1e-5, 100),
                cpar=(0.25, 0, 1), empar=(1e-7, 1e-9, 1e-5)):
        """fit a Cole-Cole term with additional EM term to phase

        Parameters
        ----------
        ePhi : float
            absolute error of phase angle
        lam : float
            regularization parameter
        remove: bool
            remove EM term from data
        mpar, taupar, cpar, empar : list[3]
            inversion parameters (starting value, lower bound, upper bound)
            for Cole-Cole parameters (m, tau, c) and EM relaxation time (em)
        """
        self.mCC, self.phiCC = fitCCEMPhi(self.f, self.phi, ePhi, lam, mpar,
                                          taupar, cpar, empar)
        # correct EM term from data
        if remove:
            self.phiOrg = self.phi
            self.phi = self.phi + \
                np.angle(relaxationTerm(self.f, self.mCC[3]))

    def fitColeCole(self, useCond=False, **kwargs):
        """fit a Cole-Cole model to the data

        Parameters
        ----------
        useCond : bool
            use conductivity form of Cole-Cole model instead of resistivity
        """
        if useCond:  # use conductivity formulation instead of resistivity
            self.mCC, self.ampCC, self.phiCC = fitCCCC(self.f, self.amp,
                                                       self.phi, **kwargs)
            self.mCC[0] = 1. / self.mCC[0]
        else:
            self.mCC, self.ampCC, self.phiCC = fitCCC(self.f, self.amp,
                                                      self.phi, **kwargs)

    def fitDebyeModel(self, ePhi=0.001, lam=1e3, lamFactor=0.8,
                      mint=None, maxt=None, nt=None, new=False,
                      showFit=False, cType=1):
        """fit a (smooth) continuous Debye model (Debye decomposition)

        Parameters
        ----------
        ePhi : float
            absolute error of phase angle
        lam : float
            regularization parameter
        lamFactor : float
            regularization factor for subsequent iterations
        mint/maxt : float
            minimum/maximum tau values to use (else automatically from f)
        nt : int
            number of tau values (default number of frequencies * 2)
        new : bool
            new implmentation (experimental)
        showFit : bool
            show fit
        cType : int
            constraint type (1/2=smoothness 1st/2nd order, 0=minimum norm)
        """
        nf = len(self.f)
        if mint is None:
            mint = .1 / max(self.f)
        if maxt is None:
            maxt = .5 / min(self.f)
        if nt is None:
            nt = nf*2
        # discretize tau, setup DD and perform DD inversion
        self.tau = np.logspace(log10(mint), log10(maxt), nt)
        phi = self.phi
        tLin, tLog, tM = pg.RTrans(), pg.RTransLog(), pg.RTransLog()
        # pg.RTransLogLU(0., 1.)
        if new:
            reNorm, imNorm = self.zNorm()
            fDD = DebyeComplex(self.f, self.tau)
            Znorm = pg.cat(reNorm, imNorm)
            IDD = pg.RInversion(Znorm, fDD, tLog, tM, False)
            IDD.setAbsoluteError(max(Znorm)*0.003+0.01)
        else:
            fDD = DebyePhi(self.f, self.tau)
            IDD = pg.RInversion(phi, fDD, tLin, tM, True)
            IDD.setAbsoluteError(ePhi)  # 1 mrad

        fDD.regionManager().setConstraintType(cType)
        IDD.stopAtChi1(False)
        startModel = pg.RVector(nt, 0.01)
        IDD.setModel(startModel)
        IDD.setLambda(lam)
        IDD.setLambdaFactor(lamFactor)
        self.mDD = IDD.run()
        IDD.echoStatus()
        if new:
            resp = np.array(IDD.response())
            respRe = resp[:nf]
            respIm = resp[nf:]
            respC = ((1 - respRe) + respIm * 1j) * max(self.amp)
            self.phiDD = np.angle(respC)
            self.ampDD = np.abs(respC)
            if showFit:
                fig, ax = self.showData(znorm=True, nrows=3)
                ax[0].plot(self.f, respRe, 'r-')
                ax[1].plot(self.f, respIm, 'r-')
                ax[2].semilogx(self.tau, self.mDD, 'r-')
                ax[2].set_xlim(max(self.tau), min(self.tau))
                ax[2].set_ylim(0., max(self.mDD))
                ax[2].grid(True)
                ax[2].set_xlabel(r'$\tau$ [s]')
                ax[2].set_xlabel('$m$ [-]')
        else:
            self.phiDD = IDD.response()
            if showFit:
                fig, ax = self.showData(nrows=3)
                ax[2].semilogx(self.tau, self.mDD, 'r-')

    def totalChargeability(self):
        """total chargeability from as Debye curve as curve integral"""
        return sum(self.mDD)

    def logMeanTau(self):
        """mean logarithmic relaxation time as 50% cumulative log curve"""
        return exp(np.sum(np.log(self.tau) * self.mDD) / sum(self.mDD))

    def showAll(self, save=False):
        """plot spectrum, Cole-Cole fit and Debye distribution"""
        # generate title strings
        if hasattr(self, 'mCC'):
            tstr = r'CC: m={:.3f} $\tau$={:.1e}s c={:.2f}'
            mCC = self.mCC
            if len(mCC) > 3:
                tstr += ' $\tau_2$={:.1e}s'
            if mCC[0] > 1:
                tstr = r'CC: $\rho$={:.1f} m={:.3f} $\tau$={:.1e}s c={:.2f}'
            tCC = tstr.format(*mCC)
        fig, ax = plt.subplots(nrows=2+hasattr(self, 'mDD'), figsize=(12, 12))
        fig.subplots_adjust(hspace=0.25)
        # amplitude
        showAmplitudeSpectrum(ax[0], self.f, self.amp, label='data', ylog=0)
        if hasattr(self, 'ampDD'):
            ax[0].plot(self.f, self.ampDD, 'm-', label='DD response')
        if hasattr(self, 'ampCC'):
            ax[0].semilogx(self.f, self.ampCC, 'r-', label='CC+E model')
        ax[0].legend(loc='best')
        # phase
        if hasattr(self, 'phiOrg'):
            ax[1].semilogx(self.f, self.phiOrg * 1e3, 'c+-', label='org. data')
        ax[1].semilogx(self.f, self.phi * 1e3, 'b+-', label='data')
        if hasattr(self, 'phiCC'):
            ax[1].semilogx(self.f, self.phiCC * 1e3, 'r-', label='CC+E model')
        if hasattr(self, 'phiDD'):
            ax[1].semilogx(self.f, self.phiDD * 1e3, 'm-', label='DD model')
        ax[1].grid(True)
        ax[1].legend(loc='best')
        ax[1].set_xlabel('f [Hz]')
        ax[1].set_ylabel('-phi [mrad]')
        if hasattr(self, 'mCC'):
            ax[1].set_title(tCC, loc='left')
        if hasattr(self, 'mDD'):  # relaxation time
            mtot = self.totalChargeability()
            lmtau = self.logMeanTau()
            tDD = r'DD: m={:.3f} $\tau$={:.1e}s'.format(mtot, lmtau)
            ax[2].semilogx(self.tau, self.mDD * 1e3)
            ax[2].set_xlim(ax[2].get_xlim()[::-1])
            ax[2].grid(True)
            ax[2].set_xlabel(r'$\tau$ [ms]')
            ax[2].set_ylabel('m [mV/V]')
            ax[2].set_title(tDD, loc='left')
        if save:
            if isinstance(save, str):
                savename = save
            else:
                savename = self.basename + '.pdf'
            fig.savefig(savename, bbox_inches='tight')
        plt.show(block=False)
        return fig, ax

if __name__ == "__main__":
    filename = 'sipexample.txt'
    sip = SIPSpectrum(filename)
#    sip.showData(znorm=True)
    if True:  # Pelton
        sip.fitCCEM()
    else:
        sip.removeEpsilonEffect()
        sip.fitColeCole(useCond=False)
#    sip.showDataKK()
    # %%
    sip.fitDebyeModel()  # , showFit=True)
    # %% create titles and plot data, fit and model
    sip.showAll(save=True)
    plt.show()
