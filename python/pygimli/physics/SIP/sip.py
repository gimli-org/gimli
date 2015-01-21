import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import pygimli as pg
from pygimli.mplviewer import createColorbar

import pybert as pb
from pybert.dataview import drawDataAsPatches


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

class SIP():
    '''
    class for managing spectral induced polarisation data
    '''
    def __init__(self, filename=None, verbose=True, **kwargs):
        ''' init function with optional data load '''
        if filename is not None:
            self.load(filename)

    def __repr__(self):  # for print function
        out = ['SIP data: nf=' + str(len(self.freq)) + ' nc=' +
               str(self.nc) + ' nv=' + str(self.nv)]
        for key in self.header:
            val = self.header[key]
            if isinstance(val, int) or isinstance(val, float):
                out.append(key+' = '+str(val))
            else:
                out.append(key+' = array('+str(val.shape)+')')
        return "\n".join(out)

    def load(self, filename):
        ''' load data from file '''
        self.header, self.DATA, self.AB, self.RU = readSIP256file(filename)
        self.basename = filename.replace('.res', '').replace('.RES', '')
        self.nc = self.header['Number_of_Readings']
        self.nv = self.header['Number_of_Remote_Units']
        self.organiseSIP256data()

    def organiseSIP256data(self):
        ''' builds up empty data container with the quadrupoles '''
        # retrieve frequencies
        self.freq = []
        for line in self.header['FrequencyParameter']:
            if line[-1] == 1.:
                self.freq.append(line[0].round(3))
        self.freq = np.array(self.freq)
        # assemble measurement logics
        aa, bb, mm, nn, ii, iu = [], [], [], [], [], []
        for ir in range(len(self.DATA)):
            readings = self.header['Readings'][ir]
            leftout = readings[3:]
            iA, iB = self.AB[ir]
            ru = self.RU[ir]
            for iru in range(len(ru)):
                iM = ru[iru]
                iN = iM + 1
                while iN in leftout:
                    iN += 1
                if iM > iB and iN-iM == iB-iA:
                    aa.append(iA)
                    bb.append(iB)
                    mm.append(iM)
                    nn.append(iN)
                    ii.append(ir)
                    iu.append(iru)

        self.ABMN = np.column_stack((aa, bb, mm, nn))
        # create data container
        self.data = pb.DataContainerERT()
#        self.data = pg.DataContainer()
#        self.data.setSensorIndexOnFileFromOne(True)
#        for c in 'abmn':
#            self.data.registerSensorIndex(c)
        for line in self.header['Layout']:
            self.data.createSensor(pg.RVector3(line[1], 0., 0.))
        self.data.resize(len(aa))
        self.data.set('a', np.array(aa)-1)
        self.data.set('b', np.array(bb)-1)
        self.data.set('m', np.array(mm)-1)
        self.data.set('n', np.array(nn)-1)
        self.data.markValid(self.data('a') > -1)
        self.data.set('k', pb.geometricFactor(self.data))
        # assemble data matrices
        self.RHOA = np.ones((self.data.size(), len(self.freq)))*np.nan
        self.PHIA = self.RHOA * 1.
        for i in range(len(ii)):
            A = self.DATA[ii[i]][iu[i]]
            for ifr in range(len(self.freq)):
                line = A[A[:, 0].round(3) == self.freq[ifr]]
                if len(line):
                    self.RHOA[i, ifr] = line[0][1]
                    self.PHIA[i, ifr] = -line[0][2]

    def writeAllData(self):
        ''' output the data as matrices '''
        fmt = '%d\t%d\t%d\t%d'
        for i in range(len(self.freq)):
            fmt += '\t%.2f'
        np.savetxt(self.basename + '_rhoa.all',
                   np.column_stack((self.ABMN, self.RHOA)), fmt=fmt)
        np.savetxt(self.basename + '_phia.all',
                   np.column_stack((self.ABMN, self.PHIA)), fmt=fmt)

    def singleFrequencyData(self, ifr, kmax=None):
        ''' return filled data container for one frequency '''
        if isinstance(ifr, float):  # choose closest frequency
            ifr = np.argmin(np.abs(self.freq - ifr))
        data1 = pb.DataContainerERT(self.data)
        data1.set('rhoa', self.RHOA[:, ifr])
        data1.set('ip', self.PHIA[:, ifr])
        if kmax:
            data1.markInvalid(pg.abs(data1('k')) > kmax)
        data1.removeInvalid()

        return data1

    def writeSingleFrequencyData(self, kmax=None):
        ''' write single frequency data in unified data format '''
        for ifr, fri in enumerate(self.freq):
            data1 = self.singleFrequencyData(ifr, kmax=kmax)
            data1.checkDataValidity()
            if fri > 1.:
                fname = '{:02d}-{:d}Hz.dat'.format(ifr, int(fri))
            else:
                fname = '{:02d}-{:d}mHz.dat'.format(ifr, int(fri*1e3))
            data1.save(self.basename+'_'+fname, 'a b m n rhoa ip')

    def generateSpectraPDF(self, useall=False, maxphi=20., maxdist=999):
        ''' make pdf file containing all spectra '''
        pdf = PdfPages(self.basename + '-spectra.pdf')
        fig, ax = plt.subplots(figsize=(8.5, 11), nrows=2, sharex=True)
        colors = 'bgrcmyk'
        markers = ('x', 'o', 'v', '^', 's', 'p', '>', '<', '+', 'd')
        for i, Data in enumerate(self.DATA):
            act = 0
            for j, data in enumerate(Data):
                freq = data[:, 0]
                rhoa, drhoa = data[:, 1], data[:, 3]
                phi, dphi = -data[:, 2], data[:, 4]
                dist = j - self.AB[i][1]
                useit = (dist > 0) and (dist < maxdist) or useall
                if useit and np.isfinite(phi[0]) and np.isfinite(rhoa[0]):
                    col = colors[j % 7]
                    ma = markers[j / 7]
                    ax[0].errorbar(freq, np.abs(rhoa), yerr=drhoa*rhoa/100.,
                                   color=col, label=str(j), marker=ma)
                    ax[0].set_xscale('log')
                    ax[0].set_yscale('log')
                    ax[1].errorbar(freq, phi, yerr=dphi,
                                   color=col, label=str(j), marker=ma)
                    ax[1].set_xscale('log')
                    act += 1

            if act:
                if i == 0:
                    ax[0].legend(loc='upper right', numpoints=1, ncol=3)
                ax[0].set_xlabel('f in Hz')
                ax[1].set_xlabel('f in Hz')
                ax[0].set_ylabel(r'$\rho_a$ in $\Omega$m')
                ax[1].set_ylabel(r'-$\phi_a$ in grad')
                ax[0].set_title(str(self.AB[i][0])+'-'+str(self.AB[i][1]))
                ax[1].set_ylim(0., maxphi)
                fig.savefig(pdf, format='pdf')
                ax[0].cla()
                ax[1].cla()

        pdf.close()

    def generateDataPDF(self, kmax=None, ipmax=15, rmin=None, rmax=None):
        ''' generate a pdf file with appar. resistivity&phase pseudosections'''
        xl = self.header['Layout'][[0, -1], 1]
        pdf = PdfPages(self.basename + '-data.pdf')
        fig, ax = plt.subplots(nrows=2, figsize=(8, 11.5), sharex=True)
        for fri in self.freq[::-1]:
            fstr = '(f={:d}Hz)'.format(int(fri))
            if fri < 1.:
                fstr = '(f={:d}mHz)'.format(int(fri*1e3))
            data = self.singleFrequencyData(fri, kmax=kmax)
            if rmin is None:
                rmin = min(data('rhoa'))
            if rmax is None:
                rmax = max(data('rhoa'))
            gci = drawDataAsPatches(ax[0], data, data('rhoa'), schemetype=5)
            createColorbar(gci, nLevs=5, cMin=rmin, cMax=rmax)
            ax[0].set_title('apparent resistivity ' + fstr)
            ax[0].set_xlim(xl)
            gci = drawDataAsPatches(ax[1], data, data('ip'), schemetype=5,
                                    logScale=False)
            createColorbar(gci, nLevs=5, cMin=0, cMax=ipmax)
            ax[1].set_title('apparent phase ' + fstr)
            fig.savefig(pdf, format='pdf')

        pdf.close()


if __name__ == "__main__":
    sip = SIP(sys.argv[1])
    print(sip)
    sip.generatePDF()
