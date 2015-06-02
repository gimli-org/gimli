# encoding: utf-8

"""
meshtools.

This package provides functions for mesh generation and modification.

.. note::
    Although we discriminate here between grids (structured meshes) and meshes
    (unstructured), both objects are treated the same internally.
"""

from .grid import *
from .mesh import *

from . polytools import createRectangle
from . polytools import createWorld
from . polytools import createCircle
from . polytools import createLine
from . polytools import mergePLC
from . polytools import readPLC

