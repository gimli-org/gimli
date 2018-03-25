# -*- coding: utf-8 -*-
"""pyGIMLi - An open-source library for modelling and inversion in geophysics

Import convention:

.. code-block:: python

    import pygimli as pg
"""

# py 2.7 compatiblity
from __future__ import division, print_function

import locale
import sys


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

################################################################################
# Please leave this block here until the following issue is fixed:
# https://github.com/ContinuumIO/anaconda-issues/issues/1068
if "conda" in __path__[0]:
    try:
        import PyQt5
        import matplotlib
        matplotlib.use("qt5agg", warn=False)
    except ImportError:
        pass
################################################################################
from . import core
from ._version import get_versions
from .core import *
from .testing import test

__version__ = get_versions()['version']
del get_versions

def version():
    """Shortcut to show and return current version."""
    logger.info('Version: ' + __version__ + " core:" + pg.versionStr())
    return __version__
