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
    for mid, col, row, scale, mati in zip(mat.entries_matrixID(),
                                          mat.entries_colStart(),
                                          mat.entries_rowStart(),
                                          mat.entries_scale(),
                                          mat.matrices()):
        widthy = mati.rows() - 0.1 # to make sure non-matrix regions are not connected in the plot
        widthx = mati.cols() - 0.1
        plc = pg.meshtools.createRectangle([col, row],
                                           [col + widthx, row + widthy],
                                           marker=mid)
        plcs.append(plc)

    bm = pg.meshtools.mergePLC(plcs)
    gci, cBar = pg.viewer.mpl.drawPLC(ax, bm, fitView=False)
    ax.invert_yaxis()
    cBar.set_label("Matrix ID")
    if len(mat.matrices()) > 10:
        gci.set_cmap("viridis")
    return gci, cBar
