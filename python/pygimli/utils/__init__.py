# -*- coding: utf-8 -*-
"""
Useful utility functions.
"""

from .base import (rms, rmsWithErr, nanrms, createDateTimeString, 
                   createfolders, #remove me Nameing convention
                   createFolders, createResultFolder,
                   getSavePath, gmat2numpy, interperc, interpExtrap, inthist,
                   nanrms, num2str, numpy2gmat, rrms, rms, rndig, saveResult)

# compatibility for fw_cleaning
from .base import rmsWithErr as rmswitherr

from .cache import (cache, strHash)
from .geostatistics import (computeInverseRootMatrix, covarianceMatrix,
                            generateGeostatisticalModel)
from .gps import GKtoUTM, findUTMZone, getProjection, getUTMProjection, readGPX
from .hankel import hankelFC
from .postinversion import iterateBounds, modCovar
from .sparseMat2Numpy import (convertCRSIndex2Map, sparseMatrix2Array,
                              sparseMatrix2coo, sparseMatrix2csr)

from .units import (unit, cMap)

from .utils import (ProgressBar, boxprint, chi2, cumDist, diff, dist,
                    filterIndex, filterLinesByCommentStr, findNearest,
                    getIndex, grange, logDropTol, niceLogspace, prettyFloat,
                    rand, randN, trimDocString, unicodeToAscii, unique,
                    unique_everseen, unique_rows, uniqueAndSum)

__all__ = [name for name in dir() if '_' not in name]
