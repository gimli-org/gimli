# -*- coding: utf-8 -*-
"""
    Utility functions for miscellaneous stuff.
"""

from .base import gmat2numpy
from .base import numpy2gmat
from .base import rndig
from .base import num2str
from .base import inthist
from .base import interperc
from .base import interpExtrap
from .base import saveResult
from .base import getSavePath
from .base import createResultFolder
from .base import createDateTimeString
from .base import createfolders

from .utils import ProgressBar
from .utils import boxprint
from .utils import trimDocString
from .utils import unicodeToAscii
from .utils import logDropTol
from .utils import grange
from .utils import niceLogspace
from .utils import diff
from .utils import dist
from .utils import cumDist
from .utils import chi2
from .utils import randN
from .utils import rand
from .utils import getIndex
from .utils import filterIndex
from .utils import findNearest
from .utils import unique_everseen
from .utils import unique
from .utils import unique_rows

from .postinversion import iterateBounds
from .postinversion import modCovar

from .geostatistics import covarianceMatrix
from .geostatistics import computeInverseRootMatrix

from . hankel import hankelFC

__all__ = [name for name in dir() if '_' not in name]
