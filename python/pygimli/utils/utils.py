#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Collection of several utility functions."""

import sys
from math import sqrt, floor
import numpy as np
import pygimli as pg


class ProgressBar(object):
    """Animated text-based progressbar.

    Animated text-based progressbar for intensive loops. Should work in the
    console and in the IPython Notebook.

    Parameters
    ----------
    its : int
        Number of iterations of the process.
    width : int
        Width of the ProgressBar, default is 80.
    sign : str
        Sign used to fill the bar.

    Examples
    --------
    >>> from pygimli.utils import ProgressBar
    >>> pbar = ProgressBar(its=20, width=40, sign='+')
    >>> pbar.update(5)
    \r[+++++++++++      30%                  ]  6 of 20 complete
    """

    def __init__(self, its, width=80, sign=":"):
        """Constructor."""
        self.its = int(its)
        self.width = width
        self.sign = sign[0]  # take first character only if sign is longer
        self.pbar = '[]'
        self._amount(0)

    def update(self, iteration):
        """Update ProgressBar by iteration number starting at 0."""
        self._setbar(iteration + 1)
        print("\r" + self.pbar)
        sys.stdout.flush()

    def _setbar(self, elapsed_it):
        """Reset pbar based on current iteration number."""
        self._amount((elapsed_it / float(self.its)) * 100.0)
        self.pbar += '  %d of %s complete' % (elapsed_it, self.its)

    def _amount(self, new_amount):
        """Calculate amount by which to update the pbar."""
        pct_done = int(round((new_amount / 100.0) * 100.0))
        full_width = self.width - 2
        num_signs = int(round((pct_done / 100.0) * full_width))
        self.pbar = '[' + self.sign * num_signs + \
            ' ' * (full_width - num_signs) + ']'
        pct_place = (len(self.pbar) // 2) - len(str(pct_done))
        pct_string = '%d%%' % pct_done
        self.pbar = self.pbar[0:pct_place] + \
            (pct_string + self.pbar[pct_place + len(pct_string):])


def boxprint(s, width=80, sym="#"):
    """Print string centered in a box.

    Examples
    --------
    >>> from pygimli.utils import boxprint
    >>> boxprint("This is centered in a box.", width=40, sym='+')
    ++++++++++++++++++++++++++++++++++++++++
    +      This is centered in a box.      +
    ++++++++++++++++++++++++++++++++++++++++
    """
    row = sym * width
    centered = s.center(width - 2)
    print("\n".join((row, centered.join((sym, sym)), row)))


def trimDocString(docstring):
    """Return properly formatted docstring.

    From: https://www.python.org/dev/peps/pep-0257/

    Examples
    --------
    >>> from pygimli.utils import trimDocString
    >>> docstring = '    This is a string with indention and whitespace.   '
    >>> trimDocString(docstring).replace('with', 'without')
    'This is a string without indention and whitespace.'
    """
    if not docstring:
        return ''
    # Convert tabs to spaces (following the normal Python rules)
    # and split into a list of lines:
    lines = docstring.expandtabs().splitlines()
    # Determine minimum indentation (first line doesn't count):
    indent = 2**16 - 1
    for line in lines[1:]:
        stripped = line.lstrip()
        if stripped:
            indent = min(indent, len(line) - len(stripped))
    # Remove indentation (first line is special):
    trimmed = [lines[0].strip()]
    if indent < 2**16 - 1:
        for line in lines[1:]:
            trimmed.append(line[indent:].rstrip())
    # Strip off trailing and leading blank lines:
    while trimmed and not trimmed[-1]:
        trimmed.pop()
    while trimmed and not trimmed[0]:
        trimmed.pop(0)
    # Return a single string:
    return '\n'.join(trimmed)


def unicodeToAscii(text):
    """TODO DOCUMENTME."""
    if isinstance(text, str):
        return text.encode("iso-8859-1", "ignore")
    else:
        return text


def logDropTol(p, droptol=1e-3):
    """Create logarithmic scaled copy of p.

    Examples
    --------
    >>> from pygimli.utils import logDropTol
    >>> x = logDropTol((-10, -1, 0, 1, 100))
    >>> print(x.array())
    [-4. -3.  0.  3.  5.]
    """
    tmp = pg.RVector(p)

    tmp = pg.abs(tmp / droptol)
    tmp.setVal(1.0, pg.find(tmp < 1.0))

    tmp = pg.log10(tmp)
    tmp *= pg.sign(p)
    return tmp

def niceLogspace(vMin, vMax, nDec=10):
    """Create nice logarithmic space from the next decade lower to vMin to
    decade larger then vMax.

    Parameters
    ----------
    vMin : float
        lower limit need to be > 0
    vMax : float
        upper limit need to be >= vMin
    nDec : int
        Amount of logarithmic equidistant steps for one decade

    Examples
    --------
    >>> from pygimli.utils import niceLogspace
    >>> v1 = niceLogspace(vMin=0.1, vMax=0.1, nDec=1)
    >>> print(v1)
    [ 0.1  1. ]
    >>> v1 = niceLogspace(vMin=0.09, vMax=0.11, nDec=1)
    >>> print(v1)
    [ 0.01  0.1   1.  ]
    >>> v1 = niceLogspace(vMin=0.09, vMax=0.11, nDec=10)
    >>> print(len(v1))
    21
    >>> print(v1)
    [ 0.01        0.01258925  0.01584893  0.01995262  0.02511886  0.03162278
      0.03981072  0.05011872  0.06309573  0.07943282  0.1         0.12589254
      0.15848932  0.19952623  0.25118864  0.31622777  0.39810717  0.50118723
      0.63095734  0.79432823  1.        ]
    """
    if vMin > vMax or vMin < 1e-12:
        print("vMin:", vMin, "vMax", vMax)
        raise Exception('vMin > vMax or vMin <= 0.')

    vmin = 10**np.floor(np.log10(vMin))
    vmax = 10**np.ceil(np.log10(vMax))

    if vmax == vmin:
        vmax *= 10

    n = np.log10(vmax/vmin)*nDec + 1

    q = 10**(1/nDec)
    print(vmin, vmax, n)

    return vmin * q**np.arange(n)



def grange(start, end, dx=0, n=0, log=False):
    """Create array with possible increasing spacing.

    Create either array from start step-wise filled with dx until end reached
    [start, end] (like np.array with defined end) n or array filled from
    start to end with n steps. [start, end] (like np.linespace) n or an
    array with with logarithmic spacing if n is given, dx will be ignored.

    Parameters
    ----------
    start: float
        First value of the resulting array
    end: float
        Last value of the resulting array
    dx: float
        Linear step length, n will be ignored
    n: int
        Amount of steps
    log: bool

    Examples
    --------
    >>> from pygimli.utils import grange
    >>> v1 = grange(start=0, end=10, dx=3)
    >>> v2 = grange(start=0, end=10, n=3)
    >>> print(v1)
    <class 'pygimli.core._pygimli_.RVector'> 4 [0.0, 3.0, 6.0, 9.0]
    >>> print(v2)
    <class 'pygimli.core._pygimli_.RVector'> 3 [0.0, 5.0, 10.0]

    Returns
    -------
    ret: :gimliapi:`GIMLI::RVector`
        Return resulting array
    """
    s = float(start)
    e = float(end)
    d = float(dx)

    if dx != 0:
        if end < start and dx > 0:
            # print("grange: decreasing range but increasing dx, swap dx sign")
            d = -d
        if end > start and dx < 0:
            # print("grange: increasing range but decreasing dx, swap dx sign")
            d = -d
        ret = pg.RVector(range(int(floor(abs((e - s) / d)) + 1)))
        ret *= d
        ret += s
        return ret

    elif n > 0:
        if not log:
            return grange(start, end, dx=(e - s) / (n - 1))
        else:
            raise Exception('not yet implemented.')

    else:
        raise Exception('Either dx or n have to be given.')


def diff(v):
    """Calculate approximate derivative.

    Calculate approximate derivative from v as d = [v_1-v_0, v2-v_1, ...]

    Parameters
    ----------
    v : array(N) | pg.R3Vector(N)
        Array of double values or positions

    Returns
    -------
    d : [type(v)](N-1) |
        derivative array

    Examples
    --------
    >>> import pygimli as pg
    >>> from pygimli.utils import diff
    >>> p = pg.R3Vector(4)
    >>> p[0] = [0.0, 0.0]
    >>> p[1] = [0.0, 1.0]
    >>> print(diff(p)[0])
    RVector3: (0.0, 1.0, 0.0)
    >>> print(diff(p)[1])
    RVector3: (0.0, -1.0, 0.0)
    >>> print(diff(p)[2])
    RVector3: (0.0, 0.0, 0.0)
    >>> p = pg.RVector(3)
    >>> p[0] = 0.0
    >>> p[1] = 1.0
    >>> p[2] = 2.0
    >>> print(diff(p))
    <class 'pygimli.core._pygimli_.RVector'> 2 [1.0, 1.0]
    """
    d = None

    if isinstance(v, np.ndarray):
        if v.ndim == 2:
            v = pg.R3Vector(v)
    elif isinstance(v, list):
        v = pg.R3Vector(v)

    if isinstance(v, pg.R3Vector):
        d = pg.R3Vector(len(v) - 1)
    else:
        d = pg.RVector(len(v) - 1)

    for i, _ in enumerate(d):
        d[i] = v[i + 1] - v[i]
    return d


def dist(p, c=None):
    """Calculate the distance for each position in p relative to pos c(x,y,z).

    Parameters
    ----------
    p : ndarray(N,2) | ndarray(N,3) | pg.R3Vector

        Position array
    c : [x,y,z] [None]
        relative origin. default = [0, 0, 0]

    Returns
    -------
    d : ndarray(N)
        Distance array

    Examples
    --------
    >>> import pygimli as pg
    >>> from pygimli.utils import dist
    >>> import numpy as np
    >>> p = pg.R3Vector(4)
    >>> p[0] = [0.0, 0.0]
    >>> p[1] = [0.0, 1.0]
    >>> print(dist(p))
    [ 0.  1.  0.  0.]
    >>> x = pg.RVector(4, 0)
    >>> y = pg.RVector(4, 1)
    >>> print(dist(np.array([x, y]).T))
    [ 1.  1.  1.  1.]
    """
    if c is None:
        c = pg.RVector3(0.0, 0.0, 0.0)
    d = np.zeros(len(p))
    pI = None
    for i, _ in enumerate(p):
        if isinstance(p[i], pg.RVector3):
            pI = p[i]
        else:
            pI = pg.RVector3(p[i])
        d[i] = (pI-c).abs()

    return d


def cumDist(p):
    """The progressive i.e, cumulative length for the path p.

    d = [0.0, d[0]+ | p[1]-p[0] |, d[1] + | p[2]-p[1] | + ...]

    Parameters
    ----------
    p : ndarray(N,2) | ndarray(N,3) | pg.R3Vector

        Position array

    Returns
    -------
    d : ndarray(N)
        Distance array

    Examples
    --------
    >>> import pygimli as pg
    >>> from pygimli.utils import cumDist
    >>> import numpy as np
    >>> p = pg.R3Vector(4)
    >>> p[0] = [0.0, 0.0]
    >>> p[1] = [0.0, 1.0]
    >>> p[2] = [0.0, 1.0]
    >>> p[3] = [0.0, 0.0]
    >>> print(cumDist(p))
    [ 0.  1.  1.  2.]
    """
    d = np.zeros(len(p))
    d[1:] = np.cumsum(dist(diff(p)))
    return d


def chi2(a, b, err, trans=None):
    """Return chi square value."""
    if trans is None:
        trans = pg.RTrans()
    d = (trans(a) - trans(b)) / trans.error(a, err)
    return pg.dot(d, d) / len(d)


def randN(n, minVal=0.0, maxVal=1.0):
    """Create RVector of length n with normally distributed random numbers."""
    r = pg.RVector(n)
    pg.randn(r)
    r *= (maxVal-minVal)
    r += minVal
    return r


def rand(n, minVal=0.0, maxVal=1.0):
    """Create RVector of length n with normally distributed random numbers."""
    r = pg.RVector(n)
    pg.rand(r, minVal, maxVal)
    return r


def getIndex(seq, f):
    """TODO DOCUMENTME."""
    # DEPRECATED_SLOW
    idx = []
    if isinstance(seq, pg.RVector):
        for i, _ in enumerate(seq):
            v = seq[i]
            if f(v):
                idx.append(i)
    else:
        for i, d in enumerate(seq):
            if f(d):
                idx.append(i)
    return idx


def filterIndex(seq, idx):
    """TODO DOCUMENTME."""
    if isinstance(seq, pg.RVector):
        # return seq(idx)
        ret = pg.RVector(len(idx))
    else:
        ret = list(range(len(idx)))

    for i, ix in enumerate(idx):
        ret[i] = seq[ix]

    return ret


def findNearest(x, y, xp, yp, radius=-1):
    """TODO DOCUMENTME."""
    idx = 0
    minDist = 1e9
    startPointDist = pg.RVector(len(x))
    for i, _ in enumerate(x):
        startPointDist[i] = sqrt(
            (x[i] - xp) * (x[i] - xp) + (y[i] - yp) * (y[i] - yp))

        if startPointDist[i] < minDist and startPointDist[i] > radius:
            minDist = startPointDist[i]
            idx = i
    return idx, startPointDist[idx]


def unique_everseen(iterable, key=None):
    """Return iterator of unique elements ever seen with preserving order.

    Return iterator of unique elements ever seen with preserving order.

    From: http://docs.python.org/library/itertools.html#recipes

    Examples
    --------
    >>> from pygimli.utils import unique_everseen
    >>> s1 = 'AAAABBBCCDAABBB'
    >>> s2 = 'ABBCcAD'
    >>> list(unique_everseen(s1))
    ['A', 'B', 'C', 'D']
    >>> list(unique_everseen(s2, key=str.lower))
    ['A', 'B', 'C', 'D']

    See also
    --------
    unique, unique_rows
    """
    try:
        from itertools import ifilterfalse
    except BaseException as _:
        from itertools import filterfalse

    seen = set()
    seen_add = seen.add
    if key is None:
        try:
            for element in ifilterfalse(seen.__contains__, iterable):
                seen_add(element)
                yield element
        except BaseException as _:
            for element in filterfalse(seen.__contains__, iterable):
                seen_add(element)
                yield element
    else:
        for element in iterable:
            k = key(element)
            if k not in seen:
                seen_add(k)
                yield element


def unique(a):
    """Return list of unique elements ever seen with preserving order.

    Examples
    --------
    >>> from pygimli.utils import unique
    >>> unique((1,1,2,2,3,1))
    [1, 2, 3]

    See also
    --------
    unique_everseen, unique_rows
    """
    return list(unique_everseen(a))


def unique_rows(array):
    """Return unique rows in a 2D array.

    Examples
    --------
    >>> from pygimli.utils import unique_rows
    >>> import numpy as np
    >>> A = np.array(([1,2,3],[3,2,1],[1,2,3]))
    >>> unique_rows(A)
    array([[1, 2, 3],
           [3, 2, 1]])
    """
    b = array.ravel().view(np.dtype((np.void,
                                     array.dtype.itemsize*array.shape[1])))
    _, unique_idx = np.unique(b, return_index=True)

    return array[np.sort(unique_idx)]
    # A_1D = A.dot(np.append(A.max(0)[::-1].cumprod()[::-1][1:], 1))
    # sort_idx = A_1D.argsort()
    # mask = np.append(True, np.diff(A_1D[sort_idx]) !=0 )
    # return A[sort_idx[np.nonzero(mask)[0][np.bincount(mask.cumsum()-1)==1]]]
