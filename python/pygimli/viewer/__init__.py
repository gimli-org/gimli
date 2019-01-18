# -*- coding: utf-8 -*-
"""Viewer interface for 2D visualizations."""

from ..mplviewer import hold
import matplotlib

# Set global hold if mpl inline backend is used (as in Jupyter Notebooks)
if 'inline' in matplotlib.get_backend():
    hold(1)

from ..mplviewer import wait

from .showmesh import plt, show, showBoundaryNorm, showMesh
