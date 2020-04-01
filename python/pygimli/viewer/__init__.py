# -*- coding: utf-8 -*-
"""Interface for 2D and 3D visualizations."""

import matplotlib

from .mpl import hold, wait
from .showmesh import plt, show, showBoundaryNorm, showMesh

# Set global hold if mpl inline backend is used (as in Jupyter Notebooks)
if 'inline' in matplotlib.get_backend():
    hold(1)
