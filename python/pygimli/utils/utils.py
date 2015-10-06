#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pygimli.utils - Collection of several utility functions.
"""
import sys
import pygimli as pg
from importlib import import_module
import numpy as np
from math import sqrt, floor

class ProgressBar(object):

    """
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
        self.its = int(its)
        self.width = width
        self.sign = sign[0]  # take first character only if sign is longer
        self.pbar = '[]'
        self._amount(0)

    def update(self, iteration):
        """Update ProgressBar by iteration number starting at 0."""
        self._setbar(iteration + 1)
        print("\r" + self.pbar, end='')
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


def opt_import(module, requiredTo="use the full functionality"):
    """
    Import and return module only if it exists.

    If `module` cannot be imported, a warning is printed followed by the
    `requiredFor` string. Otherwise, the imported `module` will be returned.
    This function should be used to import optional dependencies in order to
    avoid repeated try/except statements.

    Parameters
    ----------
    module : str
        Name of the module to be imported.
    requiredFor : str, optional
        Info string for the purpose of the dependency.

    Examples
    --------
    >>> from pygimli.utils import opt_import
    >>> pg = opt_import("pygimli")
    >>> pg.__name__
    'pygimli'
    >>> opt_import("doesNotExist", requiredTo="do something special")
    No module named 'doesNotExist'.
    You need to install this optional dependency to do something special.
    """

    # set default message for common imports
    if not requiredTo and "matplotlib" in module:
        requiredTo = "visualize 2D content"

    if module.count(".") > 2:
        raise ImportError("Can only import modules and sub-packages.")

    try:
        mod = import_module(module)
    except ImportError:
        msg = ("No module named \'%s\'.\nYou need to install this optional "
               "dependency to %s.")
        print(msg % (module, requiredTo))
        mod = None

    return mod


def trimDocString(docstring):
    """
    Return properly formatted docstring.

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
    if isinstance(text, str):
        return text.encode("iso-8859-1", "ignore")
    else:
        return text


def logDropTol(p, droptol=1e-3):
    """

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


def grange(start, end, dx=0, n=0, log=False, verbose=False):
    """
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
    <class 'pygimli._pygimli_.RVector'> 4 [0.0, 3.0, 6.0, 9.0]
    >>> print(v2)
    <class 'pygimli._pygimli_.RVector'> 3 [0.0, 5.0, 10.0]

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
    """Calculate approximate derivative from v as d = [v_1-v_0, v2-v_1, ...]
    
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
    >>> p = pg.R3Vector(4)
    >>> p[0] = [0.0, 0.0]
    >>> p[1] = [0.0, 1.0]
    >>> print(diff(p)[0], diff(p)[1], diff(p)[2])    
    >>> 
    >>> p = pg.RVector(3)
    >>> p[0] = 0.0
    >>> p[1] = 1.0
    >>> p[2] = 2.0
    >>> print(diff(p))    
    
    """
    d = None
    
    if isinstance(v, pg.R3Vector):
        d = pg.R3Vector(len(v) - 1)
    else:
        d = pg.RVector(len(v) - 1)
        
    for i in range(len(d)):
        d[i] = v[i + 1] - v[i]
    return d


def dist(p):
    """Calculate the distance for each position in p to the origin.
       
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
    >>> p = pg.R3Vector(4)
    >>> p[0] = [0.0, 0.0]
    >>> p[1] = [0.0, 1.0]
    >>> print(dist(p))
    >>> 
    >>> x = pg.RVector(4, 0)
    >>> y = pg.RVector(4, 1)
    >>> print(dist(np.array([x, y]).T))
    """
    
    d = np.zeros(len(p))
    pi = None
    for i in range(len(p)):
        if isinstance(pi, pg.RVector3):
            pi = p[i]
        else:
            pi = pg.RVector3(p[i])
        d[i] = pi.abs()
        
    return d

def xyToLength(x, y):
    raise("please use utils.distSum")
    
def distSum(p):
    """The progressive length for the path p.
    
    d = [0.0, d[0]+ |p[1]-p[0]|, d[1] + |p[2]-p[1]| + ...]
    
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
    >>> p = pg.R3Vector(4)
    >>> p[0] = [0.0, 0.0]
    >>> p[1] = [0.0, 1.0]
    >>> p[2] = [0.0, 1.0]
    >>> p[3] = [0.0, 0.0]
    >>> print(distSum(p), [0, 1, 1, 2])    
    
        
    """
    d = np.zeros(len(p))
    d[1:] = np.cumsum(dist(diff(p)))
    return d
    
def getIndex(seq, f):
    # DEPRECATED_SLOW
    idx = []
    if isinstance(seq, pg.RVector):
        for i in range(len(seq)):
            v = seq[i]
            if f(v):
                idx.append(i)
    else:
        for i, d in enumerate(seq):
            if f(d):
                idx.append(i)
    return idx


def filterIndex(seq, idx):
    if isinstance(seq, pg.RVector):
        # return seq(idx)
        ret = pg.RVector(len(idx))
    else:
        ret = list(range(len(idx)))

    for i, id in enumerate(idx):
        ret[i] = seq[id]
    return ret


def findNearest(x, y, xp, yp, radius=-1):
    idx = 0
    minDist = 1e9
    startPointDist = pg.RVector(len(x))
    for i in range(len(x)):
        startPointDist[i] = sqrt(
            (x[i] - xp) * (x[i] - xp) + (y[i] - yp) * (y[i] - yp))

        if startPointDist[i] < minDist and startPointDist[i] > radius:
            minDist = startPointDist[i]
            idx = i
    return idx, startPointDist[idx]


def unique_everseen(iterable, key=None):
    """
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
    unique
    """
    try:
        from itertools import ifilterfalse
    except:
        from itertools import filterfalse

    seen = set()
    seen_add = seen.add
    if key is None:
        try:
            for element in ifilterfalse(seen.__contains__, iterable):
                seen_add(element)
                yield element
        except:
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
    """
    Return list of unique elements ever seen with preserving order.

    Examples
    --------
    >>> from pygimli.utils import unique
    >>> unique((1,1,2,2,3,1))
    [1, 2, 3]

    See also
    --------
    unique_everseen
    """
    return list(unique_everseen(a))


def arrayToStdVectorUL(theArray):
    """Converts a 'ndarray' to pygimli.stdVectorUL."""
    vec = pg.stdVectorUL()
    for i in theArray:
        vec.append(int(i))
    return vec

if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    
    
    
