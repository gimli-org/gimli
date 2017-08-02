#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Time Domain Electromagnetics (TDEM) functions and class"""
import sys
from math import pi
import numpy as np
import matplotlib.pyplot as plt


def rhoafromU(UbyI, t, Tx, Rx=None):
    """Apparent resistivity curve from classical TEM (U or dB/dt)

    rhoafromU(U/I, t, TXarea[, RXarea])
    .. math:: \rho_a = ( (A_{Rx} *A_{Tx} * \mu_0)/ (20 U/I) )^2/3*t^{-5/3}*4e-7
    """
    if Rx is None:
        Rx = Tx  # assume single/coincident loop

    mu0 = 4e-7 * pi
    rhoa = (Rx * Tx * mu0 / 20. / UbyI)**(2. / 3.) * \
        t**(-5. / 3.) * mu0 / pi
    return rhoa


def rhoafromB(B, t, Tx, I=1):
    """Apparent resistivity from B-field TEM

    .. math:: \rho_a = ( (A_{Tx}*I*\mu_0 ) / (30B) )^2/3 * 4e-7 / t
    """
    mu0 = 4e-7 * pi
    rhoa = (I * Tx * mu0 / 30. / B)**(2. / 3.) * mu0 / pi / t
    return rhoa


# TODO: better derive a class TEMsounding from dict and put functions in there
def TxArea(snd):
    """ return effective transmitter area """
    if isinstance(snd['LOOP_SIZE'], str):
        Tx = np.prod([float(a) for a in snd['LOOP_SIZE'].split()])
    else:
        Tx = snd['LOOP_SIZE']

    return Tx


def RxArea(snd):
    """Return effective receiver area."""
    Rx = 0  # just in case of error
    if 'COIL_SIZE' in snd:
        Rx = snd['COIL_SIZE']
        if Rx == 700.:
            Rx = 100.  # hack for wrong turns in NMR noise loop
    else:  # no coil size given ==> COI or SIN ==> take loop size
        Rx = TxArea(snd)
    return Rx


def get_rhoa(snd, cal=260e-9, corrramp=True):
    """Compute apparent resistivity from sounding (usf) dict."""
    Tx = TxArea(snd)
    Rx = RxArea(snd)
    if 'COIL_SIZE' in snd:
        Rx = snd['COIL_SIZE']
    else:
        Rx = Tx

    v = snd['VOLTAGE']
    istart, istop = 0, len(v)  # default: take all
    mav = np.arange(len(v))[v == max(v)]
    if len(mav) > 1:  # several equal big ones: start after
        istart = max(mav) + 1

    if min(v) < 0.0:  # negative values: stop at first
        istop = np.argmax(v[20:] < 0.0) + 20

    print(istart, istop)
    v = v[istart:istop]
    if 'ST_DEV' in snd:
        dv = snd['ST_DEV'][istart:istop]  # / snd['CURRENT']
    else:
        dv = v * 0.01

    t = snd['TIME'][istart:istop]
    if corrramp and 'RAMP_TIME' in snd:
        t = t - snd['RAMP_TIME']

    if Rx == 1:  # apparently B-field not dB/dt
        rhoa = rhoafromB(v * cal, t, Tx)
    else:
        rhoa = rhoafromU(v, t, Tx, Rx)

    rhoaerr = dv / v * (2. / 3.)
    return rhoa, t, rhoaerr


def readusffile(filename):
    """Read data from single USF (universal sounding file) file

    Examples
    --------
        DATA = readusffile(filename)
        DATA = readusffile(filename, DATA) will append to DATA
    """

    DATA = []
    columns = []
    nr = 0
    sounding = {}
    sounding['FILENAME'] = filename
    isdata = False
    fid = open(filename)
    for line in fid:
        zeile = line.rstrip('\n').replace(',', '')  # commas useless here
        if zeile:  # anything at all
            if zeile[0] == '/':  # comment-like
                if zeile[1:4] == 'END':  # end of a sounding
                    if isdata:  # already read some data
                        sounding['data'] = columns
                        for i, cn in enumerate(sounding['column_names']):
                            sounding[cn] = columns[:, i]

                        sounding['FILENAME'] = filename
                        if 'INSTRUMENT' in sounding and 'ST_DEV' in sounding:
                            if 'terraTEM' in sounding['INSTRUMENT']:
                                sounding['ST_DEV'] *= 0.01
                                print('taking default stdev')
                        DATA.append(sounding)
                        sounding = {}

                    isdata = not isdata  # turn off data mode
                elif zeile.find(':') > 0:  # key-value pair
                    key, value = zeile[1:].split(':')
                    try:
                        val = float(value)
                        sounding[key] = val
                    except:
                        sounding[key] = value

            else:
                if isdata:
                    values = zeile.split()
                    try:
                        for i, v in enumerate(values):
                            columns[nr, i] = float(v)

                        nr += 1
                    except:
                        sounding['column_names'] = values
                        columns = np.zeros((int(sounding['POINTS']),
                                            len(values)))
                        nr = 0

    fid.close()
    return DATA


def readusffiles(filenames):
    """Read all soundings data from a list of usf files

    Example
    -------
        DATA = readusffiles(filenames)
    """
    from glob import glob
    if isinstance(filenames, str):
        if filenames.find('*') >= 0:
            filenames = glob(filenames)
        else:
            filenames = [filenames]

    DATA = []
    for onefile in filenames:
        DATA.extend(readusffile(onefile))

    return DATA


def readTEMfastFile(temfile):
    """ReadTEMfastFile(filename) reads TEM-fast file into usf sounding."""
    snd = {}
    snd['FILENAME'] = temfile
    fid = open(temfile)
    for i in range(4):
        zeile = fid.readline()
    snd['STACK_SIZE'] = int(zeile.split()[3])
    snd['RAMP_TIME'] = float(zeile.split()[5])*1e-6
    snd['CURRENT'] = float(zeile.split()[7][2:])
    zeile = fid.readline()
    fid.close()
    snd['LOOP_SIZE'] = float(zeile.split()[2])**2
    snd['COIL_SIZE'] = float(zeile.split()[5])**2
    t, v, e, r = np.loadtxt(temfile, skiprows=8, usecols=(1, 2, 3, 4),
                            unpack=True)
    ind = np.nonzero((r > 0) * (v > 0) * (t > snd['RAMP_TIME']*1.2e6))  # us
    snd['TIME'] = t[ind] * 1e-6  # us
    snd['VOLTAGE'] = v[ind]
    snd['ST_DEV'] = e[ind]
    snd['RHOA'] = r[ind]

    return snd


def readUniKTEMData(filename):
    """Read TEM data format of University of Cologne."""
    if '*' in filename:
        from glob import glob
        allfiles = glob(filename)
    else:
        allfiles = [filename]

    DATA = []
    for filename in allfiles:
        snd = {}
        snd['FILENAME'] = filename
        A = np.loadtxt(filename)
        snd['TIME'] = A[:, 1]
        snd['VOLTAGE'] = A[:, 2]
        snd['ST_DEV'] = A[:, 4] / 100 * A[:, 2]
        DATA.append(snd)
    return DATA


def readSiroTEMData(fname):
    """Read TEM data from siroTEM instrument dump.

    Example
    -------
        DATA = readSiroTEMData(filename)
        .. list of soundings with USF and siro-specific keys
    """
    Time_ST = np.array([487., 887., 1287., 1687., 2087., 2687., 3487., 4287.,
                        5087., 5887., 7087., 8687., 10287., 11887., 13487.,
                        15887., 19087., 22287., 25487., 28687., 33487., 39887.,
                        46287., 52687., 59087., 68687., 81487., 94287.,
                        107090., 119890., 139090., 164690., 190290., 215890.,
                        241490., 279890., 331090., 382290., 433490., 484690.,
                        561490., 663890., 766290., 868690., 971090., 1124700.,
                        1329500., 1534300., 1739100., 1943900.])
    Time_ET = np.array([0.05, 0.1, 0.15, 0.25, 0.325, 0.425, 0.525, 0.625,
                        0.725, 0.875, 1.075, 1.275, 1.475, 1.675, 1.975,
                        2.375, 2.775, 3.175, 3.575, 4.175, 4.975, 5.775,
                        6.575, 7.375, 8.575, 10.175, 11.775, 13.375, 14.975,
                        17.375, 20.575, 23.775, 26.975, 30.175, 34.975, 41.375,
                        47.775, 54.175, 60.574, 70.175, 82.975, 95.775,
                        108.575, 121.375, 140.575, 166.175, 191.775, 217.375,
                        242.975, 281.375, 332.575])

    fid = open(fname)
    # read in file header until : sign
    line = 'a'
    while len(line) > 0 and line[0] != ':':
        line = fid.readline()

    DATA = []
    line = fid.readline()
    while line[0] != ';':
        header = line[1:-6].split(',')

        snd = {}  # dictionary, uppercase corresponds to USF format keys
        snd['INSTRUMENT'] = 'siroTEM'
        snd['dtype'] = int(header[3])
        dstring = header[1]
        snd['DATE'] = int('20' + dstring[6:8] + dstring[3:4] + dstring[0:1])
        snd['win0'], snd['win1'], ngain, snd['conf'], snd['nch'] = \
            [int(h) for h in header[5:10]]
        snd['SOUNDING_NUMBER'] = int(header[10])

        snd['GAIN_FACTOR'] = [0.1, 1.0, 10.0, 100.0][ngain]  # predefined gains
        snd['STACK_SIZE'] = int(header[14])
        snd['ttype'] = int(header[20])
        # 1=composite, 2=earlytime, 3=standard, 4=highresolution
        snd['CURRENT'] = float(header[17])
        snd['RAMP_TIME'] = float(header[18]) * 1e-6
        snd['TIME_DELAY'] = float(header[19])
        snd['LOOP_SIZE'] = float(header[21])
        snd['COIL_SIZE'] = float(header[22])

        fid.readline()
        data = []
        line = fid.readline()[:-1]  # trim CR+LF newline
        while len(line) > 0:
            while line[-1] == '/':
                line = line[:-1] + fid.readline()[:-1].replace('\t', '')
#                aline = line

            nums = [float(el[-7:-2]) * 10**(float(el[-2:])) for el in
                    line[1:-5].split(',')[1:]]
            data.append(np.array(nums))
            line = fid.readline().rstrip('\n').rstrip('\r')

        snd['VOLTAGE'] = data[0]
        if snd['ttype'] == 2:  # early time
            snd['TIME'] = Time_ET[snd['win0'] - 1:snd['win1']] * 1e-3
        if snd['ttype'] == 3:  # standard time
            snd['TIME'] = Time_ST[snd['win0'] - 1:snd['win1']] * 1e-6

        snd['ST_DEV'] = data[1]
        if snd['dtype'] > 0:  # normal measurement
            DATA.append(snd)

        line = fid.readline()

    fid.close()
#    DATA['FILENAME'] = fname  # makes no sense as DATA is an array->snd?
    return DATA


def getname(snd):
    """Generate label name from filename entry."""
    fname = snd['FILENAME']
    name = fname[fname.rfind('\\')+1:-4]
    if 'STACK_SIZE' in snd:
        name += '-' + str(int(snd['STACK_SIZE']))

    return name


class TDEM():
    """TEM class mainly for holding data etc."""

    def __init__(self, filename=None):
        """Initialize class and (optionally) load data"""
        self.DATA = []
        self.names = []

        if filename:
            self.load(filename)

    def load(self, filename):
        """Road data from usf, txt (siroTEM), tem (TEMfast) or UniK file."""
        if filename.lower().endswith('.usf'):
            self.DATA.extend(readusffiles(filename))
        elif filename.lower().endswith('.txt'):
            self.DATA = readSiroTEMData(filename)
        elif filename.lower().endswith('.tem'):
            self.DATA = readTEMfastFile(filename)
        elif filename.lower().endswith('.dat'):  # dangerous
            self.DATA = readUniKTEMData(filename)

    def __repr__(self):
        return "<TDEMdata: %d soundings>" % (len(self.DATA))

    def showInfos(self):  # only for old scripts using it
        print(self.__repr__)

    def plotTransients(self, ax=None, **kwargs):
        """Plot all transients into one window"""
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()

        cols = 'rgbmcyk'
        pl = []
        for i, data in enumerate(self.DATA):
            t = data['TIME']
            u = data['VOLTAGE'] / RxArea(data)
            col = cols[i % len(cols)]
            pl.append(ax.loglog(t, u, marker='.', label=getname(data),
                                color=col), **kwargs)
            if 'ST_DEV' in data:
                err = data['ST_DEV'] / RxArea(data)
                ax.errorbar(t, u, yerr=err, color=col)
#                uU = u + err
#                uL = u - err
#                ax.errorbar(t, u, yerr=[uL, uU], color=col)

            if 'RAMP_TIME' in data:
                ax.vlines(data['RAMP_TIME'], min(u), max(u), colors=col)

        ax.set_xlabel('t [s]')
        ax.set_ylabel('U/I [V/A]')
        ax.legend(loc='best')
#        xlim = [10e-6, 2e-3]
        ax.grid(True)
        return fig, ax

    def plotRhoa(self, ax=None, ploterror=False, **kwargs):
        """Plot all apparent resistivity curves into one window."""
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()

        cols = 'rgbmcyk'
        for i, data in enumerate(self.DATA):
            rhoa, t, err = get_rhoa(data)
            err[err > .99] = .99
            col = cols[i % len(cols)]
            ax.loglog(t, rhoa, marker='.', label=getname(data), color=col,
                      **kwargs)
            if ploterror:
                ax.errorbar(t, rhoa, yerr=rhoa * err, color=col)

        ax.set_xlabel('t [s]')
        ax.set_ylabel(r'$\rho_a$ [$\Omega$m]')
        ax.legend(loc='best')
        ax.grid(True)
        return fig, ax

if __name__ == '__main__':
    print("do some tests here")
    tem = TDEM(sys.argv[1])
    print(tem)
    tem.plotTransients()
    tem.plotRhoa()
