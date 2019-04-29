# -*- coding: utf-8 -*-
"""
Useful utility functions.
"""

from .base import (rms, rmsWithErr, nanrms, createDateTimeString, createfolders, createResultFolder,
                   getSavePath, gmat2numpy, interperc, interpExtrap, inthist,
                   num2str, numpy2gmat, rndig, saveResult)

# backward compatibility
from .base import rmsWithErr as rmswitherr

from .geostatistics import (computeInverseRootMatrix, covarianceMatrix,
                            generateGeostatisticalModel)
from .hankel import hankelFC
from .postinversion import iterateBounds, modCovar
from .sparseMat2Numpy import (convertCRSIndex2Map, sparseMatrix2Array,
                              sparseMatrix2csr, sparseMatrix2coo)

from .utils import (ProgressBar, boxprint, 
                    prettyFloat, chi2, cumDist, diff, dist,
                    filterIndex, filterLinesByCommentStr, findNearest,
                    getIndex, grange, logDropTol, niceLogspace, rand, randN,
                    trimDocString, unicodeToAscii, unique, unique_everseen,
                    unique_rows, uniqueAndSum)

from .cache import cache

from .gps import readGPX, findUTMZone, getUTMProjection, getProjection, GKtoUTM

from .units import unit

__all__ = [name for name in dir() if '_' not in name]
