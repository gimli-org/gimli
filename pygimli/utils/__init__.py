# -*- coding: utf-8 -*-
"""
Useful utility functions.
"""
from .base import (
    rms,
    rmsWithErr,
    nanrms,
    createDateTimeString,
    createfolders,  # renamed 20200515
    createFolders,
    createResultFolder,  # renamed 20200515
    createResultPath,
    createPath,
    getSavePath,
    gmat2numpy,
    interperc,
    interpExtrap,
    inthist,
    nanrms,
    num2str,
    numpy2gmat,
    rrms,
    rms,
    rndig,
    saveResult,
    chi2,
)
# compatibility for dev #can be removed? (20200515)
from .base import rmsWithErr as rmswitherr
from .cache import cache, strHash, noCache
from .complex import (
    isComplex,
    toComplex,
    toPolar,
    squeezeComplex,
    toRealMatrix,
    KramersKronig,
)
from .dem import DEM
from .geostatistics import (
    computeInverseRootMatrix,
    covarianceMatrix,
    generateGeostatisticalModel,
)
from .gps import GKtoUTM, findUTMZone, getProjection, getUTMProjection, readGPX
from .hankel import hankelFC
from .postinversion import iterateBounds, modelCovariance, modelResolutionMatrix
from .sparseMat2Numpy import (
    convertCRSIndex2Map,
    sparseMatrix2Array,
    sparseMatrix2coo,
    sparseMatrix2csr,
    sparseMatrix2Dense,
    toSparseMatrix,
    toSparseMapMatrix,
    toCSR,
    toCOO,
)
from .streams import streamline, streamlineDir
from .units import cmap as cMap  # for compatibilty (will be removed)
from .units import unit, cmap
from .utils import (
    ProgressBar,
    boxprint,
    cumDist,
    cut,
    diff,
    dist,
    filterIndex,
    filterLinesByCommentStr,
    findNearest,
    getIndex,
    grange,
    logDropTol,
    niceLogspace,
    prettify,
    prettyFloat,
    prettyTime,
    randn,
    rand,
    trimDocString,
    unicodeToAscii,
    unique,
    unique_everseen,
    unique_rows,
    uniqueAndSum,
)

__all__ = [name for name in dir() if "_" not in name]
