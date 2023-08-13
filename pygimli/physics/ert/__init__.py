# -*- coding: utf-8 -*-
"""Direct current electromagnetics

    This package contains tools, modelling operators, and managers for
    Electrical Resistivity Tomography (ERT) & Induced polarization (IP)
"""

import pygimli as pg
from .ert import (simulate, estimateError,
                  createGeometricFactors, createInversionMesh)
from .ertManager import ERTManager
from .ertIPManager import ERTIPManager
from .ertModelling import ERTModelling, ERTModellingReference
from .ertScheme import createData
from pygimli.physics.ves import VESManager  # backward compatibility
from pygimli.physics.ves.vesModelling import VESModelling
# , VESCModelling, VESRhoModelling

from .visualization import showERTData, drawERTData, generateDataPDF
from .importData import load


@pg.renamed(createData, '1.3')  # 20210302
def createERTData(*args, **kwargs):
    pass

showData = showERTData
show = showERTData  # better create a function that can also handle mgr
pg.core.DataContainerERT.show.__doc__ = showERTData.__doc__
geometricFactor = pg.core.geometricFactors
geometricFactors = geometricFactor

# Module prototypes
DataContainer = pg.core.DataContainerERT
Manager = ERTManager
coverageERT = pg.core.coverageDCtrans
