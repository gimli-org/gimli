# -*- coding: utf-8 -*-
"""Collection of mathematical functions."""

import numpy as np
from .core import (angle, besselI0, besselI1, besselK0, besselK1, cos,
                   cot, det, dot, exp, exp10, imag, log, log10, max, median,
                   min, pow, rand, randn, real, rms, round, rrms, sign,
                   sin, sqrt, sum, toComplex, unique)
<<<<<<< HEAD
=======


def symlog(x, tol=1e-12, linearSpread=0):
    """Symmetric bi-logarithmic transformation (as used in matplotlib).

    Transforms a signed values in a logarithmic way preserving the sign.
    All absolute values below a certain threshold are treated zero or linearly
    distributed (if linearSpread>0).

    .. math::

        f(x) = sign(x) * (log10(1 + abs(x)/tol) + s/2)

    Parameters
    ----------
    x : iterable
        array to be transformed
    tol : float [None]
        tolerance for minimum values to be treated zero (or linear)
    linearSpread : float
        define how wide linear transformation is done (0-not, 1-one decade)

    Returns
    -------
    y : np.array
        transformed array of same size

    See Also
    --------
    pygimli.math.symlogInv
    """
    if tol is None:
        tol = np.min(np.abs(x))

    return np.sign(x) * (np.log10(1 + np.abs(x/tol))+linearSpread/2)


def symlogInv(y, tol=1e-12, linearSpread=0):
    """Inverse symlog transformation.

    .. math::

        f(y) = sign(y) * (10^(abs(y)-s/2) - 1) * tol

    Parameters
    ----------
    y : iterable
        array to be transformed
    tol : float [None]
        tolerance for minimum values to be treated zero (or linear)
    linearSpread : float
        define how wide linear transformation is done (0-not, 1-one decade)

    Returns
    -------
    x : np.array
        transformed array of same size

    See Also
    --------
    pygimli.math.symlog
    """
    return (10**(np.abs(y)-linearSpread/2) - 1.) * tol * np.sign(y)
>>>>>>> dev
