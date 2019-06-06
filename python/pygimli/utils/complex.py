#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Pygimli base functions to handle complex arrays"""

from math import pi

import numpy as np
import pygimli as pg


def isComplex(vals):
    """Check numpy or pg.Vector if have complex data type"""
    z = pg.CVector(2)
    if len(vals) > 0:
        if hasattr(vals, '__iter__'):
            if isinstance(vals[0], np.complex) or \
               isinstance(vals[0], complex):
             return True
    return False

def toComplex(amp, phi=None):
    """Convert real values into into complex valued array.
    
    If no phases phi are given assuming z = amp[0:N] + i amp[N:2N].

    If phi is given in (rad) complex values are generated:
    z = amp*(cos(phi) + i sin(phi))

    Parameters
    ----------
    amp: iterable (float)
        Amplitudes or real unsqueezed real valued array.
    phi: iterable (float)
        Phases in neg radiant

    Returns
    -------
    z: ndarray(dtype=np.complex)
        Complex values
    """
    if phi is not None:
        return amp.array() * (np.cos(phi) + 1j *np.sin(phi))
    N = len(amp) // 2 
    return np.array(amp[0:N]) + 1j * np.array(amp[N:])
    #return np.array(pg.toComplex(amp[0:N], amps[N:]))

def toPolar(z):
    """Convert complex values array into amplitude and phase in radiant

    If z is real valued we assume its squeezed.

    Parameters
    ----------
    z: iterable (floats, complex)
        If z contains of floats and squeezedComplex is assumed [real, imag]
        
    Returns
    -------
    amp, phi: ndarray
        Amplitude amp and phase angle phi in radiant.
    
    """
    if isComplex(z):
        return np.abs(z), np.angle(z)
    else:
        return toPolar(toComplex(z))

def squeezeComplex(z, polar=False):
    """Squeeze complex valued array into [real, imag] or [amp, -phase(rad)]"""
    if isinstance(z, pg.CSparseMapMatrix) or \
       isinstance(z, pg.CSparseMatrix) or isinstance(z, pg.CMatrix):
        return toRealMatrix(z)

    if isComplex(z):
        vals = np.array(z)
        if polar is True:
            vals = pg.cat(toPolar[:])
        else:
            vals = pg.cat(vals.real, vals.imag)
        return vals
    return z

def toRealMatrix(C):
    """Convert complex valued matrix into a real valued Blockmatrix"""
    R = pg.BlockMatrix()
    # we store the mats to keep the GC happy after leaving the scope
    Cr = pg.real(C)
    Ci = pg.imag(C)

    rId = R.addMatrix(Cr)
    iId = R.addMatrix(Ci)

    R.addMatrixEntry(rId, 0,         0,         scale=1.0)
    R.addMatrixEntry(rId, Cr.rows(), Cr.cols(), scale=1.0)
    R.addMatrixEntry(iId, 0,         Cr.cols(), scale=-1.0)
    R.addMatrixEntry(iId, Cr.rows(), 0,         scale=1.0)
    pg._y(Cr.rows(), Cr.cols(), R.rows(), R.cols())
    return R

def KramersKronig(f, re, im, usezero=False):
    """Return real/imaginary parts retrieved by Kramers-Kronig relations.

    formulas including singularity removal according to Boukamp (1993)
    """
    from scipy.integrate import simps

    x = 2. * pi * f
    im2 = np.zeros(im.shape)
    re2 = np.zeros(im.shape)
    re3 = np.zeros(im.shape)
    drdx = np.diff(re) / np.diff(x)
    didx = np.diff(im) / np.diff(x)
    dRedx = np.hstack((drdx[0], (drdx[:-1] + drdx[1:]) / 2, drdx[-1]))
    dImdx = np.hstack((didx[0], (didx[:-1] + didx[1:]) / 2, didx[-1]))
    for num, w in enumerate(x):
        x2w2 = x**2 - w**2
        x2w2[num] = 1e-12
        
        fun1 = (re - re[num]) / x2w2
        fun1[num] = dRedx[num] / 2 / w
        im2[num] = -2./pi * w * simps(fun1, x)
        
        if usezero:
            fun2 = (im * w / x - im[num]) / x2w2
            re2[num] = 2./pi * w * simps(fun2, x)  + re[0]
        else:
            fun2 = (im * x - im[num] * w) / x2w2
            fun2[num] = (im[num] / w + dImdx[num]) / 2
            re2[num] = 2./pi * simps(fun2, x) + re[-1]

        # re3 = re2

    return re2, im2

