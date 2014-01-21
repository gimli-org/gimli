#!/usr/bin/env python
# encoding: utf-8

"""
meshtools.

This package provides functions for mesh generation and modification.

.. note::
    Altough we discriminate here between grids (structured meshes) and meshes
    (unstructured), both objects are treated the same internally.
"""

from .grid import *
from .mesh import *
