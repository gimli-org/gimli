#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import matplotlib.pyplot as plt

import pygimli as pg


def readusffile(filename, data=None):
    """Read data from single USF (universal sounding file) file.

    data = readusffile( filename )
    data = readusffile( filename, data ) will append to data
    """
    if not data:
        data = []

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
                    if isdata:  # have already read some data
                        sounding['data'] = columns
                        for i, cn in enumerate(sounding['column_names']):
                            sounding[cn] = columns[:, i]

                        data.append(sounding)
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
    return data


def readusffiles(filenames):
    """ read all soundings data from a list of usf files
        DATA = readusffiles( filenames ) """
    DATA = []
    for onefile in filenames:
        DATA = readusffile(onefile, DATA)

    return DATA


def importMaxminData(filename, verbose=False):
    """import function reading in positions, data, frequencies, geometry."""
    delim = None
    fid = open(filename)
    coilspacing = 0.
    freq = []
    for i, aline in enumerate(fid):
        if aline.split()[0][0].isdigit():  # number found
            break
        elif aline.upper().find('COIL') > 0:  # [:6] == '/ COIL':
            coilspacing = float(aline.split()[-2])
        elif aline.upper().find('FREQ') > 0:  # [:6] == '/ FREQ':
            freq = np.array([float(aa) for aa in aline[aline.find(
                ':') + 1:].replace(',', ' ').split() if aa[0].isdigit()])

    fid.close()

    if verbose:
        print("CS=", coilspacing, "F=", freq)
    if aline.find(',') > 0:
        delim = ','

    nf = len(freq)
    if verbose:
        print("delim=", delim, "nf=", nf)
    A = np.loadtxt(filename, skiprows=i, delimiter=delim).T
    x, IP, OP = A[0], A[2:nf * 2 + 2:2].T, A[3:nf * 2 + 2:2].T

    return x, freq, coilspacing, IP, OP

# Data file example
# / MAXMIN ELECTROMAGNETIC SURVEY
# / FILENAME: Muell150.xyz
# / PROJECT NUMBER:         1
# / OPERATOR NUMBER:        3
# / MAXMIN EQUIPMENT: MMI-9 S/N 3395
# / SLOPE METHOD: NO SLOPES
# / COIL SEPARATION:    150.0 METRES
# / STATION SPACING:     50.0 METRES
# / MODE: MAX1 (Horizontal Coplanar)
# / FREQUENCIES ON a.dat FILE: 110, 220, 440, 880, 1760, 3520, 7040, 14080 Hz
# LINE 1.0 Surveying: N ELEV 110 220 440 880 1760 3520 7040 14080 110C 220C 440C 880C 1760C 3520C 7040C 14080C BFC ERROR
# 	0	0	6.61	7.97	8.07	12.42	14.14	19.5	27.66	28.03	45.82	23.67	63.08	11.45	68.98	-8.62	58.82	-20.77
# 20      0       5.04    5.56    7.11    10.31   13.22   16.28   25.06
# 21.91   37.18   14.17   57.3    4.67    52.07   -17.81  42.18   -31.07

if __name__ == '__main__':
    print("print do some tests here")
