# -*- coding: utf-8 -*-
"""
Useful utility functions.
"""

from .base import (createDateTimeString, createfolders, createResultFolder,
                   getSavePath, gmat2numpy, interperc, interpExtrap, inthist,
                   num2str, numpy2gmat, rndig, saveResult)
from .geostatistics import computeInverseRootMatrix, covarianceMatrix
from .hankel import hankelFC
from .postinversion import iterateBounds, modCovar
from .sparseMat2Numpy import convertCRSIndex2Map, sparseMatrix2Array
from .utils import (ProgressBar, boxprint, chi2, cumDist, diff, dist,
                    filterIndex, filterLinesByCommentStr, findNearest,
                    getIndex, grange, logDropTol, niceLogspace, rand, randN,
                    trimDocString, unicodeToAscii, unique, unique_everseen,
                    unique_rows, uniqueAndSum)

__all__ = [name for name in dir() if '_' not in name]
