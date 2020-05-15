#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions to draw various pygimli matrices with matplotlib."""

import matplotlib.pyplot as plt
import numpy as np

import pygimli as pg


def drawSparseMatrix(ax, mat, **kwargs):
    """Draw a view of a matrix into the axes.

    Parameters
    ----------
    ax : mpl axis instance, optional
        Axis instance where the matrix will be plotted.

    mat: pg.matrix.SparseMatrix or pg.matrix.SparseMapMatrix

    Returns
    -------
    ax:
    """
    mat = pg.utils.sparseMatrix2coo(mat)
    gci = ax.spy(mat)
    return gci

def drawBlockMatrix(ax, mat, **kwargs):
    """Draw a view of a matrix into the axes.

    Parameters
    ----------
    ax : mpl axis instance, optional
        Axis instance where the matrix will be plotted.

    mat: pg.Matrix.BlockMatrix

    Returns
    -------
    ax:
    """
    plcs = []
    for mid, col, row, scale in zip(mat.entries_matrixID(),
                                    mat.entries_colStart(),
                                    mat.entries_rowStart(),
                                    mat.entries_scale()):
        widthy = mat.matrices()[mid].rows()
        widthx = mat.matrices()[mid].cols()
        plc = pg.meshtools.createRectangle([col, row],
                                           [col + widthx, row+row+widthy],
                                           marker=mid)
        plcs.append(plc)

    bm = pg.meshtools.mergePLC(plcs)
    return pg.show(bm, markers=True, ax=ax)
