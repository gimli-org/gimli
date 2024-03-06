# -*- coding: utf-8 -*-
"""
Basic import of gimli core extension.
"""
import sys

# if sys.platform == 'win32':
#     os.environ['PATH'] = os.path.abspath(__file__) + ';' + os.environ['PATH']

_pygimli_ = None

try:
    # from . import _pygimli_
    # from ._pygimli_ import *

    from .libs import _pygimli_ as pgcore  
    from .libs._pygimli_ import *  

except ImportError as e:
    import traceback
    print(e)
    traceback.print_exc(file=sys.stdout)
    sys.stderr.write("ERROR: cannot import the library '_pygimli_'.\n")

