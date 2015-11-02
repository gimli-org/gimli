# -*- coding: utf-8 -*-
"""General physics independent solver."""

from . solver import *
from . green import *
from . solverFiniteVolume import *

try:
    from . solverDiscontinuousGalerkin import *
except:
    pass

class WorkSpace:
    pass
