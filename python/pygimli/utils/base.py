# -*- coding: utf-8 -*-
"""Pygimli base functions.

Some needs to be sorted.
"""
import os
import time

import numpy as np

import pygimli as pg


def gmat2numpy(mat):
    """Convert pygimli matrix into numpy.array.

    TODO implement correct rval
    """
    nmat = np.zeros((len(mat), len(mat[0])))
    for i, row in enumerate(mat):
        nmat[i] = row
    return nmat


def numpy2gmat(nmat):
    """Convert numpy.array into pygimli RMatrix.

    TODO implement correct rval
    """
    gmat = pg.RMatrix()
    for arr in nmat:
        gmat.push_back(arr)  # pg.asvector(arr))
    return gmat


def rndig(a, ndig=3):
    """Round float using a number of counting digits."""
    if np.abs(a) < 1e-4:
        return a
    else:
        return np.around(a, ndig - int(np.ceil(np.log10(np.abs(a) + 1e-4))))


def num2str(a, fmtstr='%g'):
    """List of strings (deprecated, for backward-compatibility)."""
    return [fmtstr % rndig(ai) for ai in a]


def inthist(a, vals, bins=None, islog=False):
    """What is this good for?."""
    if bins is None:
        bins = np.min((np.round(len(a) / 20), 10))

    if islog:
        hists, edges = np.histogram(np.log(a), bins=bins)
    else:
        hists, edges = np.histogram(a, bins=bins)

    cums = np.cumsum(np.hstack((0., hists))) / np.sum(hists) * 100.
    out = np.interp(vals, cums, edges)
    if islog:
        return np.exp(out)
    else:
        return out


def interperc(a, trimval=3.0, islog=False, bins=None):
    """What is this good for?."""
    return inthist(
        a, np.array([trimval, 100. - trimval]), bins=bins, islog=islog)


def interpExtrap(x, xp, yp):
    """Like np.interp function with linear extrapolation."""
    y = np.interp(x, xp, yp)
    y = np.where(x < xp[0], yp[0]+(x-xp[0])*(yp[0]-yp[1])/(xp[0]-xp[1]), y)
    return np.where(x > xp[-1], yp[-1]+(x-xp[-1])*(yp[-1]-yp[-2]) /
                    (xp[-1]-xp[-2]), y)


def saveResult(fname, data, rrms=None, chi2=None, mode='w'):
    """Save rms/chi2 results into filename."""
    with open(fname, mode) as f:
        np.savetxt(f, data)
        if rrms is not None:
            f.write('\nrrms:{}\n'.format(rrms))
        if chi2 is not None:
            f.write('\nchi2:{}\n'.format(chi2))


def getSavePath(folder=None, subfolder='', now=None):
    """TODO."""
    if folder is None:
        path = createResultFolder(subfolder, now)
    else:
        path = createfolders([folder, subfolder])
    return path


def createResultFolder(subfolder, now=None):
    """Create a result Folder."""
    result = createDateTimeString(now)
    return createfolders(['./', result, subfolder])


def createDateTimeString(now=None):
    """Return datetime as string (e.g. for saving results)."""
    if now is None:
        now = time.localtime()
    return str(now.tm_year) + str(now.tm_mon).zfill(2) + \
        str(now.tm_mday).zfill(2) + '-' + \
        str(now.tm_hour).zfill(2) + '.' + \
        str(now.tm_min).zfill(2)


def createfolders(foldername_list):
    """Create the folder structure specified by the list."""
    path = ''

    for s in foldername_list:
        if s != '/':
            path = path + s + '/'

    try:
        os.makedirs(path)
    except OSError as e:
        if os.path.exists(path):
            print('Path "{}" already exists.'.format(path))
        else:
            print('Unable to create path "{}".'.format(path))
            raise e

    return path
