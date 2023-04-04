# -*- coding: utf-8 -*-
"""Direct current electromagnetics

    This package contains tools, modelling operators, and managers for:

    * Electrical Resistivity Tomography (ERT) / Induced polarization (IP)
    * Vertical Electric Sounding (VES)
"""

import pygimli as pg
from .vesManager import VESManager
from .vesModelling import VESModelling, VESCModelling, VESRhoModelling


# Module prototypes
# DataContainer = pg.core.DataContainerERT
Manager = VESManager
