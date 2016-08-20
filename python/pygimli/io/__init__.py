# -*- coding: utf-8 -*-
"""Input and Output.

GIMLi supports different data formats. This module provides all the necessary
import and export functionality including a general load and save function with
automatic filetype detection.
"""

from .load import load
from .load import opt_import
from .gps import underlayBKGMap, GKtoUTM
