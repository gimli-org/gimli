# -*- coding: utf-8 -*-
"""Collection of mathematical functions."""

import numpy as np
from .core import (angle, besselI0, besselI1, besselK0, besselK1, cos,
                   cot, det, dot, exp, exp10, imag, log, log10, max, median,
                   min, pow, rand, randn, real, rms, round, rrms, sign,
                   sin, sqrt, sum, toComplex, unique)


def symlog(x, tol=None, linearSpread=0):
    """Symmetric bi-logarithmic transformation (as used in matplotlib).

    Transforms a signed values in a logarithmic way preserving the sign.
    All absolute values below a certain threshold are treated zero or linearly
    distributed (if linearSpread>0).

    Parameters
    ----------
    x : iterable
        array to be transformed
    tol : float [None]
        tolerance for minimum values to be treated zero (or linear)
    linearSpread : float
        define how wide linear transformation is done (0-not, 1-one decade)
    """
    if tol is None:
        tol = np.min(np.abs(x))

    return np.sign(x) * (np.log10(1 + np.abs(x/tol)+linearSpread/2))
