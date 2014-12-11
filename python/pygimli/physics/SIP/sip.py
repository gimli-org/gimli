import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def readSIP256file(resfile, verbose=True):
    """ read SIP256 file (RES format) """
    activeBlock = ''
    data = {}
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
                        data[activeBlock] = np.array(data[activeBlock])
                        activeBlock = ''
                    elif token[:5] == 'Begin':
                        activeBlock = token[6:]
                        data[activeBlock] = []
                    else:
                        value = line[line.rfind(']')+1:]
                        try:  # direct line information
                            if '.' in value:
                                num = float(value)
                            else:
                                num = int(value)
                            data[token] = num
                        except:  # maybe beginning or end of a block
                            pass
                else:
                    if activeBlock:
                        nums = np.array(line.split(), dtype=float)
                        data[activeBlock].append(nums)

    RHOA, PHASE, AB = [], [], []
    Rhoa, Phase, rhoa, phase = [], [], [], []
    for line in LINE:
        sline = line.split()
        if line.find('Reading') == 0:
            rdno = int(sline[1])
            if rdno:
                AB.append((int(sline[4]), int(sline[6])))
            if rdno > 1:
                RHOA.append(Rhoa)
                PHASE.append(Phase)
                if verbose:
                    print('Reading ' + str(rdno - 1) + ':' + str(len(Rhoa)) +
                          ' RUs')
            Rhoa, Phase = [], []
        elif line.find('Remote Unit') == 0:
            if len(rhoa):
                Rhoa.append(rhoa)
                Phase.append(phase)
#            runo = int(sline[2])  # not being used so far
            rhoa = []
            phase = []
        elif line.find('Freq') >= 0:
            pass
        elif len(sline) > 1:
            rhoa.append(float(line[10:25]))  # sline[1]))
            phase.append(float(line[25:35]))  # sline[2]))
    Rhoa.append(rhoa)
    Phase.append(phase)
    RHOA.append(Rhoa)
    PHASE.append(Phase)
    if verbose:
        print('Reading ' + str(rdno) + ':' + str(len(Rhoa)) + ' RUs')

    return data, RHOA, PHASE, AB


def showAmplitudeSpectrum(ax, freq, amp, ylabel=r'$\rho_a$ in $\Omega$m',
                          **kwargs):
    """ show amplitude spectrum """
    ax.loglog(freq, amp, **kwargs)
    ax.set_xlabel('f in Hz')
    ax.set_ylabel(ylabel)


def showPhaseSpectrum(ax, freq, phi, ylabel=r'$\phi_a$ in mrad',
                      **kwargs):
    """ show amplitude spectrum """
    ax.loglog(freq, phi, **kwargs)
    ax.set_xlabel('f in Hz')
    ax.set_ylabel(ylabel)


class SIP():
    '''
    class for managing spectral induced polarisation data
    '''
    def __init__(self, filename=None, verbose=True, **kwargs):
        ''' init function with optional data load '''
        if filename is not None:
            self.load(filename)

    def __repr__(self):  # for print function
        out = ['SIP data: nf=' + str(self.nf) + ' nc=' + str(self.nc) +\
              ' nv=' + str(self.nv)]
        for key in self.header:
            val = self.header[key]
            if isinstance(val, int) or isinstance(val, float):
                out.append(key+' = '+str(val))
            else:
                out.append(key+' = array('+str(val.shape)+')')
        return "\n".join(out)

    def load(self, filename):
        ''' load data from file '''
        self.header, self.RHOA, self.PHASE, self.AB = readSIP256file(filename)
        self.resfile = filename
        self.nf = self.header['MinFrqInd']
        self.nc = self.header['Number_of_Readings']
        self.nv = self.header['Number_of_Remote_Units']

    def generatePDF(self, useall=False, maxphi=20.):
        """ make pdf file containing all spectra """
        freq = self.header['FrequencyParameter'][:, 0]
        pdf = PdfPages(self.resfile.replace('.RES', '.pdf'))
        fig, ax = plt.subplots(figsize=(8.5, 11), nrows=2)
        colors = 'bgrcmyk'
        markers = ('x', 'o', 'v', '^', 's', 'p', '>', '<', '+', 'd')
        for i in range(len(self.PHASE)):
            Phase = self.PHASE[i]
            Rhoa = self.RHOA[i]
            act = 0
            for j in range(len(Phase)):
                rhoa = np.array(Rhoa[j])
                phi = -np.array(Phase[j])
                useit = (j > self.AB[i][0] + 1) or useall
                if useit and np.isfinite(phi[0]) and np.isfinite(rhoa[0]):
                    col = colors[act % 7]
                    ma = markers[act / 7]
                    ax[0].loglog(freq[:len(rhoa)], rhoa, color=col,
                                 label=str(j), marker=ma)
                    ax[1].semilogx(freq[:len(phi)], phi, color=col,
                                   label=str(j), marker=ma)
                    act += 1

            if act:
                ax[0].legend(loc='upper right', numpoints=1, ncol=2)
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

if __name__ == "__main__":
    sip = SIP(sys.argv[1])
    print(sip)
    sip.generatePDF()
