# -*- coding: utf-8 -*-
"""Direct current electromagnetics

    This package contains tools, modelling operators, and managers for:

    * lightweight Electrical Resistivity Tomography (ERT)
    * Vertical Electric Sounding (VES)
"""

from .ert import (ERTManager, ERTModelling, createERTData, #
                 ERTModellingReference,
                 simulate, createGeometricFactors)
from .ves import VESModelling, VESCModelling, VESManager

from .visualization import showERTData 
show = showERTData
