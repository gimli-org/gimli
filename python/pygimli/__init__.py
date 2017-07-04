# -*- coding: utf-8 -*-
"""pyGIMLi - Python package for GIMLi including bindings to C++ library.

Import convention:

.. code-block:: python

    import pygimli as pg
"""

# py 2.7 compatiblity
from __future__ import print_function, division

import sys
import locale
import os

from . import core
from ._version import get_versions
from .core import *
from .testing import test


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
# print(locale.localeconv()['decimal_point'])
# if locale.localeconv()['decimal_point'] == ',':
#   print("Found locale decimal_point ',' and change it to: decimal point '.'")
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


__version__ = get_versions()['version']


def version():
    """Shortcut to show and return current version."""
    print(__version__)
    return __version__

# __all__ = ['__version__']
# __all__.extend(core.__all__)
