# -*- coding: utf-8 -*-
"""Direct current electromagnetics

    This package contains tools, modelling operators, and managers for:

    * Electrical Resistivity Tomography (ERT) / Induced polarization (IP)
    * Vertical Electric Sounding (VES)
"""

from .ert import (createERTData, simulate, estimateError,
                  createGeometricFactors, createInversionMesh)
from .ertManager import ERTManager
from .ertModelling import ERTModelling, ERTModellingReference
from .ves import VESModelling, VESCModelling, VESManager

from .visualization import showERTData
from .importData import load

show = showERTData
createData = createERTData
