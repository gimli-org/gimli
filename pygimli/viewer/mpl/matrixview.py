#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Some data related viewer."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Rectangle, Wedge

import pygimli as pg

from .utils import updateAxes as updateAxes_


def drawMatrix(ax, mat, **kwargs):
    """Draw a view to a matrix into the axe."

    TODO
    ----
        * pg.core.BlockMatrix

    Parameters
    ----------
    ax : mpl axis instance, optional
        Axis instance where the matrix will be plotted.

    mat: obj
        obj can be so far:
        * pg.core.*Matrix
        * scipy.sparce

    Returns
    -------
    ax:
    """
    mat = pg.utils.sparseMatrix2coo(mat)
    ax.spy(mat)
    return ax

