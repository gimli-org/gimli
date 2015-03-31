from math import pi, log10, exp
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt
import pygimli as pg


def readTXTSpectrum(filename):
    """ read spectrum from ZEL device output (txt) data file """
    fid = open(filename)
    lines = fid.readlines()
    fid.close()
    f, amp, phi = [], [], []
    for line in lines[1:]:
        snums = line.replace(';', ' ').split()
        if len(snums) > 3:
            f.append(float(snums[0]))
            amp.append(float(snums[1]))
            phi.append(-float(snums[3]))
        else:
            break

    return np.asarray(f), np.asarray(amp), np.asarray(phi)


def showAmplitudeSpectrum(ax, freq, amp, ylabel=r'$\rho_a$ in $\Omega$m',
                          grid=True, marker='+', ylog=True, **kwargs):
    """ show amplitude spectrum """
    ax.semilogx(freq, amp, marker=marker, label='obs', **kwargs)
    if ylog:
        ax.set_yscale('log')
    ax.set_ylim(min(amp) * .99, max(amp * 1.01))
    ax.set_xlabel('f in Hz')
    ax.set_ylabel(ylabel)
    ax.grid(grid)


def showPhaseSpectrum(ax, freq, phi, ylabel=r'$\phi_a$ in mrad',
                      grid=True, marker='+', ylog=False, **kwargs):
    """ show phase spectrum """
    ax.semilogx(freq, phi, marker=marker, label='obs', **kwargs)
    if ylog:
        ax.set_yscale('log')
    ax.set_xlabel('f in Hz')
    ax.set_ylabel(ylabel)
    ax.grid(grid)


def plotSpectrum(ax, freq, vals, ylabel=r'$\phi_a$ in mrad',
                      grid=True, marker='+', ylog=True, **kwargs):
    """ show phase spectrum """
    ax.loglog(freq, phi, marker=marker, label='obs', **kwargs)
    if ylog:
        ax.set_yscale('log')
    ax.set_xlabel('f in Hz')
    ax.set_ylabel(ylabel)
    ax.grid(grid)


def showSpectrum(freq, amp, phi, nrows=2, ylog=None):
    """ show amplitude and phase spectra in two subplots """
    fig, ax = plt.subplots(nrows=nrows, sharex=(nrows==2))
    showAmplitudeSpectrum(ax[0], freq, amp, ylog=ylog)
    showPhaseSpectrum(ax[1], freq, phi, ylog=ylog)
    return fig, ax


def relaxationTerm(f, tau, c=1., a=1.):
    """ auxiliary function for Debye type relaxation term """
    return 1. / ((f * 2. * pi * tau * 1j)**c + 1)**a


def DebyeRelaxation(f, tau, m):
    """ complex-valued single Debye relaxation term """
    return 1. - (1.-relaxationTerm(f, tau))*m


def ColeCole(f, R, m, tau, c, a=1):
    """ Complex valued Cole-Cole model """
    return (1. - m * (1. - relaxationTerm(f, tau, c, a))) * R


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


def readSIP256file(resfile, verbose=False):
    """
    read SIP256 file (RES format) - mostly used for 2d SIP by pybert.sip

    Parameters
    ----------
        filename - *.RES file (SIP256 raw output file)
        verbose - do some output [False]

    Returns
    -------
        header - dictionary of measuring setup
        DATA - data AB-list of MN-list of matrices with f, Z, phi, dZ, dphi
        AB - list of current injection
        RU - list of remote units

    Examples
    --------
        header, DATA, AB, RU = readSIP256file('myfile.res', True)
    """
    activeBlock = ''
    header = {}
    LINE = []
    dataAct = False
    with open(resfile, 'r') as f:
        for line in f:
            if dataAct:
                LINE.append(line)
            elif len(line):
                if line[0] == '[':
                    token = line[1:line.rfind(']')].replace(' ', '_')
                    if token == 'Messdaten_SIP256':
                        dataAct = True
                    elif token[:3] == 'End':
                        header[activeBlock] = np.array(header[activeBlock])
                        activeBlock = ''
                    elif token[:5] == 'Begin':
                        activeBlock = token[6:]
                        header[activeBlock] = []
                    else:
                        value = line[line.rfind(']') + 1:]
                        try:  # direct line information
                            if '.' in value:
                                num = float(value)
                            else:
                                num = int(value)
                            header[token] = num
                        except Exception:  # maybe beginning or end of a block
                            pass
                else:
                    if activeBlock:
                        nums = np.array(line.split(), dtype=float)
                        header[activeBlock].append(nums)

    DATA, Data, data, AB, RU, ru = [], [], [], [], [], []
    for line in LINE:
        sline = line.split()
        if line.find('Reading') == 0:
            rdno = int(sline[1])
            if rdno:
                AB.append((int(sline[4]), int(sline[6])))
            if ru:
                RU.append(ru)
                ru = []
            if rdno > 1 and Data:
                Data.append(np.array(data))
                DATA.append(Data)
                Data, data = [], []
                if verbose:
                    print('Reading ' + str(rdno - 1) + ':' + str(len(Data)) +
                          ' RUs')
        elif line.find('Remote Unit') == 0:
            ru.append(int(sline[2]))
            if data:
                Data.append(np.array(data))
                data = []
        elif line.find('Freq') >= 0:
            pass
        elif len(sline) > 1 and rdno > 0:  # some data present
            for c in range(6):
                if len(sline[c]) > 11:  # too long line / missing space
                    if c == 0:
                        part1 = sline[c][:-15]
                        part2 = sline[c][10:]
                    else:
                        part1 = sline[c][:-10]
                        part2 = sline[c][9:]
                    sline = sline[:c] + [part1] + [part2] + sline[c + 1:]
            data.append(np.array(sline[:5], dtype=float))

    Data.append(np.array(data))
    DATA.append(Data)
    if verbose:
        print('Reading ' + str(rdno) + ':' + str(len(Data)) + ' RUs')

    return header, DATA, AB, RU


def KramersKronig(f, re, im, usezero=False):
    """ return real/imaginary parts retrieved by Kramers-Kronig relations

        formulas including singularity removal according to Boukamp (1993)
    """
    x = f * 2. * pi
    im2 = np.zeros(im.shape)
    re2 = np.zeros(im.shape)
    re3 = np.zeros(im.shape)
    drdx = np.diff(re) / np.diff(x)
    dredx = np.hstack((drdx[0], (drdx[:-1] + drdx[1:]) / 2, drdx[-1]))
    didx = np.diff(im) / np.diff(x)
    dimdx = np.hstack((didx[0], (didx[:-1] + didx[1:]) / 2, didx[-1]))
    for num in range(len(x)):
        w = x[num]
        x2w2 = x**2 - w**2
        x2w2[num] = 1e-12
        fun1 = (re - re[num]) / x2w2
        fun1[num] = dredx[num] / 2 / w
        im2[num] = -simps(fun1, x) * 2. * w / pi
        fun2 = (im * w / x - im[num]) / x2w2
        re2[num] = simps(fun2, x) * 2. * w / pi + re[0]
        fun3 = (im * x - im[num] * w) / x2w2
        fun3[num] = (im[num] / w + dimdx[num]) / 2
        re3[num] = simps(fun3, x) * 2. / pi + re[-1]

    if usezero:
        return re2, im2
    else:
        return re3, im2


class SIPSpectrum():
    """ SIP spectrum data analysis """
    def __init__(self, filename=None, unify=False, onlydown=True,
                 f=None, amp=None, phi=None, basename='new'):
        """init SIP class with either filename to read or data vectors

        Examples
        --------
        >>> sip = SIPSpectrum('myfilename.txt', unify=True) # unique f values
        >>> sip = SIPSpectrum(f=f, amp=R, phi=phase, basename='new')
        """
        self.basename = basename
        if filename is not None:
            if filename.rfind('.txt') > 0:
                self.basename = filename[:-4]
                self.f, self.amp, self.phi = readTXTSpectrum(filename)
        if f is not None:
            self.f = f
        if amp is not None:
            self.amp = amp
        if phi is not None:
            self.phi = phi
        if unify:
            self.unifyData()

    def unifyData(self, onlydown=False):
        fu = np.unique(self.f)
        if len(fu) < len(self.f):
            if onlydown:
                wende = min(np.nonzero(np.diff(self.f) > 0)[0])
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

    def realimag(self):
        """real and imaginary part"""
        return self.amp * np.cos(self.phi), self.amp * np.sin(self.phi)

    def zNorm(self):
        """ normalized real (difference) and imag. z (Nordsiek&Weller, 2008)"""
        re, im = self.realimag()
        R0 = max(self.amp)
        zNormRe = 1. - re / R0
        zNormIm = im / R0
        return zNormRe, zNormIm

    def showData(self, reim=False, znorm=False, nrows=2):
        """ show amplitude and phase spectrum in two subplots

        Parameters
        ----------
        reim - show real/imaginary part instead of amplitude/phase
        znorm - use normalized real/imag parts after Nordsiek&Weller (2008)
        nrows - use nrows subplots (default=2)
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

    def showDataKK(self):
        """show data as real/imag subplots along with Kramers-Kronig curves"""
        fig, ax = self.showData(reim=True)
        re, im = self.realimag()
        reKK, imKK = KramersKronig(self.f, re, im)

        ax[0].plot(self.f, reKK, label='KK')
        ax[1].plot(self.f, imKK, label='KK')
        for i in (0, 1):
            ax[i].set_yscale('linear')
            ax[i].legend()

    def showPolarPlot(self):
        """ show data in a polar plot (real against imaginary parts) """
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

    def fitCCEM(self, ePhi=0.001, lam=1000., remove=True,
                mpar=(0.2, 0, 1), taupar=(1e-2, 1e-5, 100),
                cpar=(0.25, 0, 1), empar=(1e-7, 1e-9, 1e-5)):
        """ fit a Cole-Cole term with additional EM term to phase """
        # %% Cole-Cole forward model
        fCCEM = PeltonPhiEM(self.f)
        fCCEM.region(0).setParameters(*mpar)    # m (start,lower,upper)
        fCCEM.region(1).setParameters(*taupar)  # tau
        fCCEM.region(2).setParameters(*cpar)   # c
        fCCEM.region(3).setParameters(*empar)   # tau-EM
        ICC = pg.RInversion(self.phi, fCCEM, False)  # set up inversion class
        ICC.setAbsoluteError(ePhi)  # 1 mrad
        ICC.setLambda(1000.)  # start with large damping and cool later
        ICC.setMarquardtScheme(0.8)  # lower lambda by 20%/it., no stop chi=1
        self.mCC = ICC.run()  # run inversion
        ICC.echoStatus()
        self.rCC = np.asarray(ICC.response())  # get model response for display
        # %% correct EM term from data
        if remove:
            self.phiOrg = self.phi
            self.phi = self.phi + \
                np.angle(relaxationTerm(self.f, self.mCC[3]))

    def fitDebyeModel(self, ePhi=0.001, lam=1e3, lamFactor=0.8,
                      mint=None, maxt=None, nt=None, new=False,
                      showFit=False, cType=1):
        """ fit a (smooth) continuous Debye model (Debye decomposition) """
        nf = len(self.f)
        if mint is None:
            mint = .1 / max(self.f)
        if maxt is None:
            maxt = .5 / min(self.f) * 30
        if nt is None:
            nt = nf*2
        # %% discretize tau, setup DD and perform DD inversion
        self.tau = np.logspace(log10(mint), log10(maxt), nt)
        phi = self.phi
        tLin, tLog, tM = pg.RTrans(), pg.RTransLog(), pg.RTransLogLU(0., 1.)
        if new:
            reNorm, imNorm = self.zNorm()
            fDD = DebyeComplex(self.f, self.tau)
            Znorm = pg.cat(reNorm, imNorm)
            IDD = pg.RInversion(Znorm, fDD, tLog, tM, False)
            IDD.setAbsoluteError(max(Znorm)*0.003+0.01)
        else:
            fDD = DebyePhi(self.f, self.tau)
            IDD = pg.RInversion(phi, fDD, tLin, tM, False)
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
        """ total chargeability from as Debye curve as curve integral """
        return sum(self.mDD)

    def logMeanTau(self):
        """ mean logarithmic relaxation time as 50% cumulative log curve """
        return exp(np.sum(np.log(self.tau) * self.mDD) / sum(self.mDD))

    def showAll(self):
        """ plot spectrum, Cole-Cole fit and Debye distribution """
        mtot = self.totalChargeability()
        lmtau = self.logMeanTau()
        if hasattr(self, 'mCC'):
            tstr = r'CC: m={:.3f} $\tau$={:.1e}s c={:.2f} $\tau_2$={:.1e}s'
            mCC = self.mCC
            tCC = tstr.format(mCC[0], mCC[1], mCC[2], mCC[3])
        tDD = r'DD: m={:.3f} $\tau$={:.1e}s'.format(mtot, lmtau)
        fig, ax = plt.subplots(nrows=3, figsize=(12, 12))
        fig.subplots_adjust(hspace=0.25)
        showAmplitudeSpectrum(ax[0], self.f, self.amp)
        if hasattr(self, 'ampDD'):
            ax[0].plot(self.f, self.ampDD, 'r-')
        if hasattr(self, 'phiOrg'):
            ax[1].semilogx(self.f, self.phiOrg * 1e3, 'c+-', label='org. data')
        ax[1].semilogx(self.f, self.phi * 1e3, 'b+-', label='data')
        if hasattr(self, 'rCC'):
            ax[1].semilogx(self.f, self.rCC * 1e3, 'r-', label='CC+EM model')
        ax[1].semilogx(self.f, self.phiDD * 1e3, 'm-', label='DD model')
        ax[1].grid(True)
        ax[1].legend(loc='best')
        ax[1].set_xlabel('f [Hz]')
        ax[1].set_ylabel('-phi [mrad]')
        if hasattr(self, 'mCC'):
            ax[1].set_title(tCC)
        ax[2].semilogx(self.tau, self.mDD * 1e3)
        ax[2].set_xlim(ax[2].get_xlim()[::-1])
        ax[2].grid(True)
        ax[2].set_xlabel(r'$\tau$ [ms]')
        ax[2].set_ylabel('m [mV/V]')
        ax[2].set_title(tDD)
        fig.savefig(self.basename + '.pdf', bbox_inches='tight')
        plt.show(block=False)
        return fig, ax

if __name__ == "__main__":
    # %%
    filename = 'sipexample.txt'
    sip = SIPSpectrum(filename)
#    sip.showData(znorm=True)
    sip.fitCCEM()
    # %%
    sip.fitDebyeModel(new=True)  # , showFit=True)
    # %% create titles and plot data, fit and model
    sip.showAll()
    plt.show()
