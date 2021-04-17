# -*- coding: utf-8 -*-
"""Curve-fitting with harmonic functions plus offset and drift."""
import numpy as np

import pygimli as pg


class HarmFunctor(object):
    """Functor for harmonic functions plus offset and drift."""

    def __init__(self, A, coeff, xmin, xSpan):
        """Initialize."""
        self.A_ = A
        self.coeff_ = coeff
        self.xmin_ = xmin
        self.xSpan_ = xSpan

    def __call__(self, x):
        """Yield function call."""
        nc = len(self.coeff_) / 2
        A = np.ones(nc * 2) * 0.5
        A[1] = 3. * x
        A[2::2] = np.sin(2.0 * np.pi * np.arange(1, nc) *
                         (x - self.xmin_) / self.xSpan_)
        A[3::2] = np.cos(2.0 * np.pi * np.arange(1, nc) *
                         (x - self.xmin_) / self.xSpan_)
        return sum(A * self.coeff_)


def harmfitNative(y, x=None, nc=None, xc=None, err=None):
    """Python-based curve-fit by harmonic functions.

    Parameters
    ----------
    y : iterable
        values of a curve to be fitted
    x : iterable
        abscissa, if none [0..len(y))
    nc : int
        number of coefficients
    err : iterable
        absolute data error
    xc : iterable
        abscissa to predict y on (otherwise equal to x)
    """
    y = np.asarray(y)

    if x is None:
        x = list(range(len(y)))

    x = np.asarray(x)

    if xc is not None:
        xc = np.asarray(xc)
    else:
        xc = x

    if err is None:
        err = np.ones((1, len(x)))  # ones(size(x)) end

    if nc is None or nc == 0:
        # nc=round(length(x)/30) % number of coefficients
        nc = int(len(x) / 30)

    xspan = max(x) - min(x)
    xmi = min(x)

    # A=ones(length(x),nc*2+2)/2 %nc*(sin+cos)+offset+drift
    A = np.ones((len(x), nc * 2 + 2)) / 2.0  # %nc*(sin+cos)+offset+drift
    # B=ones(length(xc),nc*2+2)/2
    B = np.ones((len(xc), nc * 2 + 2)) / 2.0

    A[:, 1] = x[:] * 3.0
    B[:, 1] = xc[:] * 3.0

    for i in range(1, nc + 1):
        # A(:,i*2+1)=sin(2*i*pi*(x-xmi)/xspan)
        # A(:,i*2+2)=cos(2*i*pi*(x-xmi)/xspan)
        # B(:,i*2+1)=sin(2*i*pi*(xc(:)-xmi)/xspan)
        # B(:,i*2+2)=cos(2*i*pi*(xc(:)-xmi)/xspan)

        A[:, i * 2 + 0] = np.sin(2.0 * i * np.pi * (x - xmi) / xspan)
        A[:, i * 2 + 1] = np.cos(2.0 * i * np.pi * (x - xmi) / xspan)
        B[:, i * 2 + 0] = np.sin(2.0 * i * np.pi * (xc - xmi) / xspan)
        B[:, i * 2 + 1] = np.cos(2.0 * i * np.pi * (xc - xmi) / xspan)

    # w = 1.0 / err w(~isfinite(w))=0
    w = 1.0 / err
    # w[(~isfinite(w)]=0

    # coeff=(spdiags(w,0,length(w),length(w))*A)\(y.*w)
    # WA = np.diag(w, 0).dot(A)
    # coeff, res, rank, s = np.linalg.lstsq(WA, y * w, rcond=None)
    coeff, res, rank, s = np.linalg.lstsq(np.diag(w, 0) * A, (y * w)[0, :],
                                          rcond=None)

    return sum((B * coeff).T), HarmFunctor(A, coeff, xmi, xspan)


def harmfit(y, x=None, error=None, nc=42, resample=None, lam=0.1,
            window=None, verbose=False, dosave=False,
            lineSearch=True, robust=False, maxiter=20):
    """GIMLi-based curve-fit by harmonic functions.

    Parameters
    ----------
    y : 1d-array - values to be fitted

    x : 1d-array(len(y)) - data abscissa data. default: [0 .. len(y))

    error : 1d-array(len(y)) error of y. default (absolute error = 0.01)

    nc : int - Number of harmonic coefficients

    resample : 1d-array - resample y to x using fitting coeffients

    window : int - just fit data inside window bounds

    Returns
    -------
    response : 1d-array(len(resample) or len(x)) - smoothed values

    coefficients : 1d-array - fitting coefficients
    """
    if x is None:
        x = np.arange(len(y))

    xToFit = None
    yToFit = None

    if window is not None:
        idx = pg.find((x >= window[0]) & (x < window[1]))
#        idx = getIndex(x , lambda v: v > window[0] and v < window[1])

        xToFit = x(idx)
        yToFit = y(idx)

        if error is not None:
            error = error(idx)
    else:
        xToFit = x
        yToFit = y

    fop = pg.core.HarmonicModelling(nc, xToFit, verbose)
    inv = pg.core.RInversion(yToFit, fop, verbose, dosave)
    # inv = pg.frameworks.MarquardtInversion(fop, verbose=verbose, debug=dosave)
    # inv.dataVals = yToFit

    if error is not None:
        inv.setAbsoluteError(error)
        # inv.errorVals = error
    else:
        inv.setRelativeError(0.01)
        # inv.errorVals =0.01

    if error is not None:
        inv.stopAtChi1(True)

    inv.setMarquardtScheme(0.9)
    inv.setLambda(lam)
    inv.setMaxIter(maxiter)
    inv.setLineSearch(lineSearch)
    inv.setRobustData(robust)
    # inv.setConstraintType(0)

    coeff = inv.run()

    if resample is not None:

        ret = fop.response(coeff, resample)

        if window is not None:
            # print pg.find((resample < window[0]) | (resample >= window[1]))
            ret.setVal(0.0, pg.find((resample < window[0]) |
                                    (resample >= window[1])))
#            idx = getIndex(resample,
#                           lambda v: v <= window[0] or v >= window[1])
#            for i in idx: ret[i] = 0.0
        return ret, coeff
    else:
        return inv.response(), coeff
# def harmfit(...)
