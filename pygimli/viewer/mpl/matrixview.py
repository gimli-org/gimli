#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions to draw various pygimli matrices with matplotlib."""

import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg


def drawSparseMatrix(ax, mat, **kwargs):
    """Draw a view of a matrix into the axes.

    Parameters
    ----------
    ax : mpl axis instance, optional
        Axis instance where the matrix will be plotted.

    mat: pg.matrix.SparseMatrix or pg.matrix.SparseMapMatrix

    Examples
    --------
    >>> import pygimli as pg
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from pygimli.viewer.mpl import drawSparseMatrix
    >>> A = np.random.randn(10,10)
    >>> SM = pg.core.SparseMapMatrix()
    >>> for i in range(10):
    ...     SM.setVal(i, i, 5.0)
    >>> fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, sharex=True)
    >>> drawSparseMatrix(ax1, A, colOffset=5, rowOffset=5)
    >>> drawSparseMatrix(ax2, SM)

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

    Examples
    --------
    >>> import pygimli as pg
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> I = pg.core.IdentityMatrix(10)
    >>> SM = pg.core.SparseMapMatrix()
    >>> for i in range(10):
    ...     SM.setVal(i, 10 - i, 5.0)
    ...     SM.setVal(i,  i, 5.0)
    >>> B = pg.core.BlockMatrix()
    >>> B.add(I, 0, 0)
    >>> B.add(SM, 10, 10)
    >>> fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    >>> _ = pg.show(B, ax=ax1)
    >>> _ = pg.show(B, spy=True, ax=ax2)

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
            mati = mat.mat(mid)
            if isinstance(mati, pg.core.IdentityMatrix):
                mati = np.eye(mati.size())
            gci.append(drawSparseMatrix(ax, mati,
                       rowOffset=row, colOffset=col, color=cMap(mid)))

        return gci, None
    else:

        plcs = []
        for mid, col, row, scale in zip(mat.entries_matrixID(),
                                        mat.entries_colStart(),
                                        mat.entries_rowStart(),
                                        mat.entries_scale()):
            widthy = mat.mat(mid).rows() - 0.1 # to make sure non-matrix regions are not connected in the plot
            widthx = mat.mat(mid).cols() - 0.1
            plc = pg.meshtools.createRectangle([col, row],
                                               [col + widthx, row + widthy],
                                               marker=mid)
            plcs.append(plc)

        bm = pg.meshtools.mergePLC(plcs)
        gci, cBar = pg.viewer.mpl.drawPLC(ax, bm, fitView=False)
        ax.invert_yaxis()
        ax.xaxis.tick_top()
        cBar.set_label("Matrix ID")

        if len(mat.matrices()) > 10:
            gci.set_cmap("viridis")
        return gci, cBar
