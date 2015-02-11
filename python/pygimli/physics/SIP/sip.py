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


if __name__ == "__main__":
    pass
