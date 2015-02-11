import numpy as np
import pygimli as pg
import matplotlib.pyplot as plt
from math import pi, log10, exp


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


def readSIP256file(resfile, verbose=False):
    """ read SIP256 file (RES format) """
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
                        value = line[line.rfind(']')+1:]
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
                    sline = sline[:c] + [part1] + [part2] + sline[c+1:]
            data.append(np.array(sline[:5], dtype=float))

    Data.append(np.array(data))
    DATA.append(Data)
    if verbose:
        print('Reading ' + str(rdno) + ':' + str(len(Data)) + ' RUs')

    return header, DATA, AB, RU


class SIPSpectrum():
    """ SIP spectrum data and analyses """
    def __init__(self, filename=None):  # initialize class
        if filename.rfind('.txt') > 0:
            self.basename = filename[:-4]
            self.f, self.amp, self.phi = readTXTSpectrum(filename)

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
            mint = .1/max(self.f)
        if maxt is None:
            maxt = .5/min(self.f)
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
        return sum(self.mDD)

    def logMeanTau(self):
        return exp(np.sum(np.log(self.tau) * self.mDD) / sum(self.mDD))

    def showAll(self):
        """ plot the spectrum with its fitting and the Debye distribution """
        mtot = self.totalChargeability()
        lmtau = self.logMeanTau()
        if hasattr(self, 'mCC'):
            tstr = r'CC: m={:.3f} $\tau$={:.1e}s c={:.2f} $\tau_2$={:.1e}s'
            mCC = self.mCC
            tCC = tstr.format(mCC[0], mCC[1], mCC[2], mCC[3])
        tDD = r'DD: m={:.3f} $\tau$={:.1e}s'.format(mtot, lmtau)
        fig, ax = plt.subplots(2, 1, figsize=(12, 10))
        fig.subplots_adjust(hspace=0.25)
        ax[0].semilogx(self.f, self.phi*1e3, 'b+-', label='data')
        ax[0].semilogx(self.f, self.rCC*1e3, 'r-', label='CC+EM model')
        if hasattr(self, 'phiC'):
            ax[0].semilogx(self.f, self.phiC*1e3, 'c+-', label='corr. data')
        ax[0].semilogx(self.f, self.rDD*1e3, 'm-', label='DD model')
        ax[0].grid(True)
        ax[0].legend(loc='best')
        ax[0].set_xlabel('f [Hz]')
        ax[0].set_ylabel('-phi [mrad]')
        if hasattr(self, 'mCC'):
            ax[0].set_title(tCC)
        ax[1].semilogx(self.tau, self.mDD*1e3)
        ax[1].set_xlim(ax[1].get_xlim()[::-1])
        ax[1].grid(True)
        ax[1].set_xlabel(r'$\tau$ [ms]')
        ax[1].set_ylabel('m [mV/V]')
        ax[1].set_title(tDD)
        fig.savefig(self.basename+'.pdf', bbox_inches='tight')
        plt.show(block=False)


if __name__ == "__main__":
    # %%
    filename = 'sipexample.txt'
    sip = SIPSpectrum(filename)
    print('Cole-Cole inversion')
    sip.fitCCEM()
    print('Debye inversion')
    sip.fitDebyeModel()
    # %% create titles and plot data, fit and model
    sip.showAll()
