# -*- coding: utf-8 -*-

"""
Imports and extensions of the C++ bindings.
"""
import os
import sys
import traceback

import numpy as np


if sys.platform == 'win32':
    os.environ['PATH'] = __path__[0] + ';' + os.environ['PATH']

_pygimli_ = None

try:
    from . import _pygimli_
    from ._pygimli_ import *
except ImportError as e:
    print(e)
    traceback.print_exc(file=sys.stdout)
    sys.stderr.write("ERROR: cannot import the library '_pygimli_'.\n")
