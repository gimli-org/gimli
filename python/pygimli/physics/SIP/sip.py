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
                          grid=True, marker='+', **kwargs):
    """ show amplitude spectrum """
    ax.loglog(freq, amp, marker=marker, label='obs', **kwargs)
    ax.set_ylim(min(amp) * .99, max(amp * 1.01))
    ax.set_xlabel('f in Hz')
    ax.set_ylabel(ylabel)
    ax.grid(grid)


def showPhaseSpectrum(ax, freq, phi, ylabel=r'$\phi_a$ in mrad',
                      grid=True, marker='+', **kwargs):
    """ show phase spectrum """
    ax.loglog(freq, phi, marker=marker, label='obs', **kwargs)
    ax.set_xlabel('f in Hz')
    ax.set_ylabel(ylabel)
    ax.grid(grid)


def showSpectrum(freq, amp, phi):
    """ show amplitude and phase spectra in two subplots """
    fig, ax = plt.subplots(nrows=2, sharex=True)
    showAmplitudeSpectrum(ax[0], freq, amp)
    showPhaseSpectrum(ax[1], freq, phi)
    return fig, ax


def relaxationTerm(f, tau, c=1.):
    """ auxiliary function for Debye type relaxation term """
    return 1. / ((f * 2. * pi * tau * 1j)**c + 1)


def ColeCole(f, R, m, tau, c):
    """ Complex valued Cole-Cole model """
    return (1. - m * (1. - relaxationTerm(f, tau, c))) * R


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

    """ Debye decomposition (smooth Debye relaxation model) """

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


def readSIP256file(resfile, verbose=False):
    """
    read SIP256 file (RES format) - can be used for 2d SIP by pybert/sip

    Parameters:
    -----------
        filename - *.RES file (SIP256 raw output file)
        verbose - do some output [False]

    Returns:
    --------
        header - dictionary of measuring setup
        DATA - data matrix
        AB - list of current injection
        RU - list of remote units

    Examples:
    --------
        header, DATA, AB, RU = readSIP256file(filename, TrueS)
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

    """ SIP spectrum data and analyses """

    def __init__(self, filename=None, unify=False, onlydown=True,
                 f=None, amp=None, phi=None, basename='new'):
        """init SIP class with either filename to read or data vectors

        Examples:
        ---------
        >>> sip = SIPSpectrum('myfilename.txt', unify=True) # unique f values
        >>> sip = SIPSpectrum(f=f, amp=R, phi=phase, basename='new')
        """
        self.basename = basename
        if filename is not None:
            if filename.rfind('.txt') > 0:
                self.basename = filename[:-4]
                self.f, self.amp, self.phi = readTXTSpectrum(filename)
                if unify:
                    self.unifyData()
                if onlydown:
                    wende = min(np.nonzero(np.diff(self.f) > 0)[0])
                    self.f = self.f[wende::-1]
                    self.amp = self.amp[wende::-1]
                    self.phi = self.phi[wende::-1]
        if f is not None:
            self.f = f
        if amp is not None:
            self.amp = amp
        if phi is not None:
            self.phi = phi

    def unifyData(self):
        fu = np.unique(self.f)
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
        return self.amp * np.cos(self.phi), self.amp * np.sin(self.phi)

    def showData(self, reim=False):
        """ show amplitude and phase spectrum in two subplots """
        if reim:
            re, im = self.realimag()
            fig, ax = showSpectrum(self.f, re, im)
            ax[0].set_ylabel('real part')
            ax[1].set_ylabel('imaginary part')
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
            self.phiC = self.phi + \
                np.angle(relaxationTerm(self.f, self.mCC[3]))

    def fitDebyeModel(self, ePhi=0.001, lam=1e5, lamFactor=0.8,
                      mint=None, maxt=None, nt=None):
        """ fit a (smooth) continuous Debye model (Debye decomposition) """
        if mint is None:
            mint = .1 / max(self.f)
        if maxt is None:
            maxt = .5 / min(self.f)
        if nt is None:
            nt = len(self.f)
        # %% discretize tau, setup DD and perform DD inversion
        self.tau = np.logspace(log10(mint), log10(maxt), nt)
        fDD = DebyePhi(self.f, self.tau)
        fDD.region(0).setConstraintType(1)
        tD, tM = pg.RTrans(), pg.RTransLogLU(0., 1.)
        phi = self.phi
        if hasattr(self, 'phiC'):  # corrected data present
            phi = self.phiC
        IDD = pg.RInversion(phi, fDD, tD, tM, False)
        IDD.setAbsoluteError(ePhi)  # 1 mrad
        IDD.stopAtChi1(False)
        startModel = pg.RVector(nt, 0.01)
        IDD.setModel(startModel)
        IDD.setLambda(lam)
        IDD.setLambdaFactor(lamFactor)
        self.mDD = IDD.run()
        IDD.echoStatus()
        self.rDD = IDD.response()

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
        fig, ax = plt.subplots(2, 1, figsize=(12, 10))
        fig.subplots_adjust(hspace=0.25)
        ax[0].semilogx(self.f, self.phi * 1e3, 'b+-', label='data')
        ax[0].semilogx(self.f, self.rCC * 1e3, 'r-', label='CC+EM model')
        if hasattr(self, 'phiC'):
            ax[0].semilogx(self.f, self.phiC * 1e3, 'c+-', label='corr. data')
        ax[0].semilogx(self.f, self.rDD * 1e3, 'm-', label='DD model')
        ax[0].grid(True)
        ax[0].legend(loc='best')
        ax[0].set_xlabel('f [Hz]')
        ax[0].set_ylabel('-phi [mrad]')
        if hasattr(self, 'mCC'):
            ax[0].set_title(tCC)
        ax[1].semilogx(self.tau, self.mDD * 1e3)
        ax[1].set_xlim(ax[1].get_xlim()[::-1])
        ax[1].grid(True)
        ax[1].set_xlabel(r'$\tau$ [ms]')
        ax[1].set_ylabel('m [mV/V]')
        ax[1].set_title(tDD)
        fig.savefig(self.basename + '.pdf', bbox_inches='tight')
        plt.show(block=False)


if __name__ == "__main__":
    # %%
    filename = 'sipexample.txt'
    sip = SIPSpectrum(filename)
    sip.showData()
    print('Cole-Cole inversion')
    sip.fitCCEM()
    print('Debye inversion')
    sip.fitDebyeModel()
    # %% create titles and plot data, fit and model
    sip.showAll()
