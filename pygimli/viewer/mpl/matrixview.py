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
    sp = pg.optImport("scipy.sparse", "visualize Block Matrices")
    pmatrix = sp.lil_matrix((mat.rows(), mat.cols()))
    nmats = len(mat.matrices())
    pg.info(nmats)
    for mid, col, row, scale in zip(mat.entries_matrixID(),
                                    mat.entries_colStart(),
                                    mat.entries_rowStart(),
                                    mat.entries_scale()):
        widthy = mat.matrices()[mid].rows()
        widthx = mat.matrices()[mid].cols()
        pmatrix[row:row+widthy, col:col+widthx] = mid + 1
    plotmat = pmatrix.toarray()
    plotmat[plotmat == 0.0] = np.nan
    plotmat -= 1
    gci = ax.matshow(plotmat, vmin=-.5, vmax=nmats-.5)
    cMap = plt.cm.get_cmap("Set3", nmats)
    cMap.set_under
    gci.set_cmap(cMap)
    return gci
