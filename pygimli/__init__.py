# -*- coding: utf-8 -*-
"""
pyGIMLi - An open-source library for modelling and inversion in geophysics
"""
import sys
import locale

from .core.decorators import (renamed, singleton)

### Import everything that should be accessible through main namespace.
from .core import (BVector, CVector, DataContainer, DataContainerERT,
                   IVector, Line, Mesh, Plane, Pos,
                   RVector3, Vector, abs, cat, center, exp, find,
                   interpolate, log, log10, logDropTol, max,
                   mean, median, min, search, setDebug, setThreadCount, sort,
                   Stopwatch, sum, trans, unique, versionStr, x, y, z, zero)

from .core import isScalar, isArray, isPos, isR3Array, isComplex, isMatrix

from .core import math # alias all from .core.math.* to pg.math.*
from .core import matrix # alias all from .core.matrix.* to pg.matrix.*

from .core.matrix import (BlockMatrix, Matrix, SparseMapMatrix, SparseMatrix)

from .core.logger import (_, _d, _g, _r, _y, _b, critical, d, debug, 
                          deprecated, renameKwarg,
                          error, info, setDebug, setLogLevel, setVerbose, v,
                          verbose, warn)

warning = warn # convenience

from .core.config import getConfigPath, rc, getCPUCount

from .meshtools import createGrid, interpolate
from .solver import solve
from .utils import boxprint, cache, cut, unique, unit, cmap, randn
from .utils import prettify as pf
from .viewer import plt, show, wait
from .frameworks import fit, Modelling, Inversion
from .testing import test#, setTestingMode, testingMode

from .core.load import getCachePath, getExampleFile
from .core.load import load, optImport

def checkAndFixLocaleDecimal_point(verbose=False):
    """
    """
    if locale.localeconv()['decimal_point'] == ',':
        if verbose:
            print("Found locale decimal_point ',' "
                  "and change it to: decimal point '.'")
    try:
        locale.localeconv()['decimal_point']
        locale.setlocale(locale.LC_NUMERIC, 'C')
    except Exception as e:
        print(e)
        print('cannot set locale to decimal point')

    # LC_CTYPE should be something with UTF-8
    # export LC_CTYPE="de_DE.UTF-8"
    # python -c 'import sys; print(sys.stdout.encoding)'

checkAndFixLocaleDecimal_point(verbose=True)
#print(locale.localeconv()['decimal_point'])
#if locale.localeconv()['decimal_point'] == ',':
  #print("Found locale decimal_point ',' and change it to: decimal point '.'")
#try:
   #locale.localeconv()['decimal_point']
   #locale.setlocale(locale.LC_NUMERIC, 'C')
#except:
   #print('cannot set locale to decimal point')


if '--debug' in sys.argv or '-d' in sys.argv:
    setDebug(True)
else:
    setDebug(False)

if '--verbose' in sys.argv or '-v' in sys.argv:
    setVerbose(True)
else:
    setVerbose(False)

# if '--test' in sys.argv or '-t' in sys.argv:
#     setTestingMode(True)
# else:
#     setTestingMode(False)


################################################################################
# Please leave this block here until the following issue is fixed:
# https://github.com/ContinuumIO/anaconda-issues/issues/1068
#if "conda" in __path__[0]:
#    try:
#        import PyQt5
#        import matplotlib
#        matplotlib.use("qt5agg", warn=False)
#    except ImportError:
#        pass
################################################################################
__version__ = "0"

def findVersion(cache=True):
    """
    Find current version either generated by versioneer or from local cache
    to avoid extensive git systemcalls.
    """
    import os
    global __version__

    setDebug(False)
    root = os.path.abspath(os.path.join(__file__, "../../"))
    gitPath = os.path.join(root, '.git')
    gitIndexFile = os.path.join(gitPath, 'index')
    versionCacheFile = os.path.join(getCachePath(), 'VERSION')
    versionPyFile = os.path.abspath(os.path.join(__file__, "_version.py"))

    loadCache = False

    if os.path.exists(versionCacheFile) and os.path.exists(gitIndexFile):
        # if git exists and cache is newer then load cache
        t1 = os.path.getmtime(versionCacheFile)
        t2 = os.path.getmtime(gitIndexFile)
        if t1 > t2:
            loadCache = True

    if os.path.exists(versionCacheFile) and os.path.exists(versionPyFile):
        # if git does not exists and cache is newer then _version.py load cache
        t1 = os.path.getmtime(versionCacheFile)
        t2 = os.path.getmtime(versionPyFile)
        if t1 > t2:
            loadCache = True

    if loadCache is True and cache is True:
        with open(versionCacheFile, 'r') as fi:
            __version__ = fi.read()
            debug('Loaded version info from cache.',
                    versionCacheFile, __version__)
        return __version__

    debug('Fetching version info.')
    from ._version import get_versions
    __versions__ = get_versions()
    __version__ = __versions__['version']

    def _get_branch():
        """Get current git branch."""
        from os.path import exists

        if exists(gitPath):
            from subprocess import check_output
            out = check_output(["git", "--git-dir", gitPath, "rev-parse",
                                "--abbrev-ref", "HEAD"]).decode("utf8")

            branch = out.split("\n")[0]
            if not "HEAD" in branch:
                return branch

        return None

    # def _get_latest_tag():
    #     from os.path import exists

    #     if exists(gitPath):
    #         from subprocess import check_output
    #         out = check_output(["git", "--git-dir", gitPath,
    #             "describe", "--tag"]).decode("utf8")

    #         tag = out.split("\n")[0].split('-')[0]
    #         return tag
    #     return None

    if __versions__["dirty"]:
        __version__ = __version__.replace(".dirty", " (with local changes")

        _branch = _get_branch()
        # print('######', _branch)
        # _tag = _get_latest_tag()
        # print('######', _tag)

        if _branch:
            __version__ += " on %s branch)" % _branch
        else:
            __version__ += ")"

    if not os.path.exists(versionCacheFile):
        os.makedirs(os.path.dirname(versionCacheFile), exist_ok=True)

    with open(versionCacheFile, 'w') as fi:
        fi.write(__version__)
        debug('Wrote version info to cache:', versionCacheFile, __version__)

# call once to get version from cache, setup or _version.py
findVersion()

def version(cache=True):
    """Shortcut to show and return current version."""
    findVersion(cache=cache)
    if cache is True:
        info('Version (cached): ' + __version__ + " core:" + versionStr())
    else:
        info('Version: ' + __version__ + " core:" + versionStr())
    return __version__

__swatch__ = dict()


def tic(msg=None, key=0):
    """Start global timer. Print elapsed time with `toc()`.

    You can start multiple stopwatches with optional identifier.

    Parameters
    ----------
    msg : string, optional
        Print message string just before starting the timer.
    key: dict key
        Identifier for your Stopwatch.
    """
    if msg:
        print(msg)
    try:
        __swatch__[key].start()
    except:
        __swatch__[key] = Stopwatch(start=True)


def toc(msg=None, box=False, reset=False, key=0):
    """Print elapsed time since global timer was started with `tic()`.

    Arguments
    ---------
    msg: string [None]
        Print message string just after printing the elapsed time. If box is
        True, then embed msg into its own box
    box: bool [False]
        Embed the time in an ascii box
    reset: bool [False]
        Reset the stopwatch.
    id: identifier
        Identifier for your Stopwatch.
    """
    if msg:
        if box is True:
            boxprint(msg)
        else:
            print(msg, end=' ')

    seconds = dur(reset=reset, key=key)
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    if h <= 0 and m <= 0:
        time = "%.2f" % s
    elif h <= 0:
        if m == 1.0:
            time = "%d minute and %.2f" % (m, s)
        else:
            time = "%d minutes and %.2f" % (m, s)
    elif h == 1.0:
        time = "%d hour, %d minutes and %.2f" % (h, m, s)
    else:
        time = "%d hours, %d minutes and %.2f" % (h, m, s)
    p = print if not box else boxprint

    if len(__swatch__.keys()) > 1:
        p("Elapsed time ({0}) is {1} seconds.".format(key, time))
    else:
        p("Elapsed time is {0} seconds.".format(time))


def dur(reset=False, key=0):
    """Return time in seconds since global timer was started with `tic()`."""
    if key in __swatch__:
        return __swatch__[key].duration(reset)
    else:
        warn("No stopwatch for id {0}".format(key))
        return 0.0
