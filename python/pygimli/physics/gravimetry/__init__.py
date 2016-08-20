# -*- coding: utf-8 -*-
"""Solve gravimetric and magneto static problems in 2d and 3D analytical"""

from . gravMagModelling import BaZCylinderHoriz
from . gravMagModelling import BZPoly
from . gravMagModelling import BaZSphere

from . gravMagModelling import uSphere
from . gravMagModelling import gradUSphere
from . gravMagModelling import gradGZSphere

from . gravMagModelling import uCylinderHoriz
from . gravMagModelling import gradUCylinderHoriz
from . gravMagModelling import gradGZCylinderHoriz

from . gravMagModelling import gradUHalfPlateHoriz
from . gravMagModelling import gradGZHalfPlateHoriz

from . gravMagModelling import solveGravimetry

from . gravMagModelling import GravimetryModelling  # FOP class
# from . gravMagModelling import Gravimetry  # MethodManager class TODO

# __all__ = ['']
