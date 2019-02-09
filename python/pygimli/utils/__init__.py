# -*- coding: utf-8 -*-
"""
Useful utility functions.
"""

from .base import (rms, rmswitherr, nanrms, createDateTimeString, createfolders, createResultFolder,
                   getSavePath, gmat2numpy, interperc, interpExtrap, inthist,
                   num2str, numpy2gmat, rndig, saveResult)

from .geostatistics import (computeInverseRootMatrix, covarianceMatrix,
                            generateGeostatisticalModel)
from .hankel import hankelFC
from .postinversion import iterateBounds, modCovar
from .sparseMat2Numpy import (convertCRSIndex2Map, sparseMatrix2Array,
                              sparseMatrix2csr, sparseMatrix2coo)

from .utils import (ProgressBar, boxprint, chi2, cumDist, diff, dist,
                    filterIndex, filterLinesByCommentStr, findNearest,
                    getIndex, grange, logDropTol, niceLogspace, rand, randN,
                    trimDocString, unicodeToAscii, unique, unique_everseen,
                    unique_rows, uniqueAndSum, )

from .gps import readGPX, findUTMZone, getUTMProjection, getProjection, GKtoUTM

__all__ = [name for name in dir() if '_' not in name]
