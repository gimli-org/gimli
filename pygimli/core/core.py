# -*- coding: utf-8 -*-
import sys

if sys.platform == 'win32':
    os.environ['PATH'] = __path__[0] + ';' + os.environ['PATH']

#_pygimli_ = None
pgcore = None

try:
    from . import _pygimli_  # if it works: as pgcore, replace all _pygimli_
    # from . import _pygimli_ as pgcore  # if it works: as pgcore, replace all _pygimli_
    pgcore = _pygimli_
    dir(pgcore)
    from ._pygimli_ import *  # check if . can be omitted
except ImportError as e:
    print("did not find in-place pg core, try import pgcore")
    # import pgcore as _pygimli_
    import pgcore
    from pgcore import _pygimli_  # check version compatibility
    from pgcore import *
except ImportError as e:
    print(e)
    traceback.print_exc(file=sys.stdout)
    sys.stderr.write("ERROR: cannot import the library '_pygimli_'.\n")
