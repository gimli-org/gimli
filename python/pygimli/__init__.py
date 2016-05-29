# -*- coding: utf-8 -*-

"""
pyGIMLi - Python package for GIMLi including bindings to C++ library.

Usage:

.. code-block:: python

    import pygimli as pg

"""

from . import core
from .core import *
#from .core._pygimli_ import *
    
import locale

def checkAndFixLocaleDecimal_point(verbose=False):
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

    #LC_CTYPE should be something with UTF-8
    # export LC_CTYPE="de_DE.UTF-8"
    #python -c 'import sys; print(sys.stdout.encoding)'

checkAndFixLocaleDecimal_point(verbose=True)
# print(locale.localeconv()['decimal_point'])
# if locale.localeconv()['decimal_point'] == ',':
#    print("Found locale decimal_point ',' and change it to: decimal point '.'")
# try:
#    locale.localeconv()['decimal_point']
#    locale.setlocale(locale.LC_NUMERIC, 'C')
# except:
#    print('cannot set locale to decimal point')

if '--debug' in sys.argv:
    print("set debug mode")
    core._pygimli_.setDebug(True)

def warnNonEmptyArgs(kwargs):
    if len(kwargs) > 0:
        print("Warning! unrecognized keyword arguments", kwargs)

def test(target=None, show=False, onlydoctests=False, coverage=False, htmlreport=False):
    """Run docstring examples and additional tests.

    Examples
    --------
    >>> from pygimli.utils import boxprint
    >>> test(target=boxprint)

    Parameters
    ----------
    target : function, optional
        Function or method to test. By default everything is tested.
    show : boolean, optional
        Show matplotlib windows during test run. They will be closed
        automatically.
    onlydoctests : boolean, optional
        Run test files in ../tests as well.
    coverage : boolean, optional
        Create a coverage report. Requires the pytest-cov plugin.
    htmlreport : str, optional
        Filename for HTML report such as www.pygimli.org/build_tests.html.
        Requires pytest-html plugin.
    """
    if target:
        import doctest
        doctest.run_docstring_examples(target, globals())
        return

    try:
        import pytest
    except ImportError:
        raise ImportError("pytest is required to run test suite. " + \
                          "Try 'sudo pip install pytest'.")

    from matplotlib import pyplot as plt
    from pygimli.utils import opt_import
    pc = opt_import("pytest_cov", "create a code coverage report")
    ph = opt_import("pytest_html", "create a html report")

    old_backend = plt.get_backend()
    if not show:
        plt.switch_backend("Agg")
    cwd = os.path.realpath(__path__[0])
    cfg = os.path.join(cwd, "../tests/setup.cfg")
    cmd = ""
    if os.path.exists(cfg):
        cmd += "-c %s " % cfg
    if pc and coverage:
        cmd += "--cov pygimli --cov-report term " + \
               "--cov-config %s " % cfg.replace("setup.cfg", ".coveragerc")
    if ph and htmlreport:
        cmd += "--html %s " % htmlreport
    cmd += "%s " % cwd
    if not onlydoctests and os.path.exists(cfg):
        cmd += os.path.join(cwd, "../tests")

    exitcode = pytest.main(cmd)
    plt.switch_backend(old_backend)
    plt.close('all')
    sys.exit(exitcode)

# provide __version__ string
from ._version import get_versions
__version__ = get_versions()['version']

def version():
    ''' Shortcut to show and return current version. '''
    print(__version__)
    return __version__

#__all__ = ['__version__']
#__all__.extend(core.__all__)  
