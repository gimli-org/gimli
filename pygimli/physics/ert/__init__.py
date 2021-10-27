# -*- coding: utf-8 -*-
"""Direct current electromagnetics

    This package contains tools, modelling operators, and managers for:

    * Electrical Resistivity Tomography (ERT) / Induced polarization (IP)
    * Vertical Electric Sounding (VES)
"""

import pygimli as pg
from .ert import (simulate, estimateError,
                  createGeometricFactors, createInversionMesh)
from .ertManager import ERTManager
from .ertModelling import ERTModelling, ERTModellingReference
from .ertScheme import createData
from .ves import VESModelling, VESCModelling, VESManager

from .visualization import showERTData
from .importData import load

@pg.renamed(createData, '1.3')  # 20210302
def createERTData(*args, **kwargs):
    pass

showData = showERTData
show = showERTData  # better create a function that can also handle mgr
geometricFactor = pg.core.geometricFactor
geometricFactors = geometricFactor
