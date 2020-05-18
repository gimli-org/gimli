#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions to draw various pygimli matrices with matplotlib."""

import numpy as np
import matplotlib as mpl

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
    row = kwargs.pop('rowOffset', 0)
    col = kwargs.pop('colOffset', 0)
    color = kwargs.pop('color', None)

    mat = pg.utils.sparseMatrix2coo(mat)
    mat.row += row
    mat.col += col
    gci = ax.spy(mat, color=color)

    ax.autoscale(enable=True, axis='both', tight=True)
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
    if kwargs.pop('spy', False):
        gci = []
        ids = mat.entries_matrixID()
        cMap = pg.plt.cm.get_cmap("Set3", len(pg.utils.unique(ids)))

        for mid, col, row, scale in zip(mat.entries_matrixID(),
                                        mat.entries_colStart(),
                                        mat.entries_rowStart(),
                                        mat.entries_scale()):
            gci.append(drawSparseMatrix(ax, mat.mat(mid),
                        rowOffset=row, colOffset=col, color=cMap(mid)))

        return gci, None
    else:

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
