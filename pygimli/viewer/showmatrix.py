# -*- coding: utf-8 -*-
"""Generic matrix visualization tools."""

import matplotlib as mpl
import numpy as np

import pygimli as pg

from .mpl import createColorBar, updateColorBar
from .mpl.matrixview import drawBlockMatrix, drawSparseMatrix


def showMatrix(mat, ax=None, **kwargs):
    """Show various pyGIMLi matrices using matplotlib.

    Args
    ----
    mat: matrix

    ax: mpl.axes

    Keyword Args
    ------------
        **kwargs : forwarded to mpl plotting commands

    Returns
    -------
        mpl.axes, Colorbar
    """
    if ax is None:
        print(ax)
        ax = pg.show()[0]

    try:
        from scipy.sparse import spmatrix

        if isinstance(mat, spmatrix):
            gci = drawSparseMatrix(ax, mat, **kwargs)
            return ax, None
    except ImportError:
        pass

    if isinstance(mat, (pg.core.RSparseMapMatrix, pg.core.RSparseMatrix)):
        gci = drawSparseMatrix(ax, mat, **kwargs)
        cBar = None
    elif isinstance(mat, pg.matrix.BlockMatrix):
        gci, cBar = drawBlockMatrix(ax, mat, **kwargs)

        if cBar is None:
            uniqueIDs = pg.unique([e.matrixID for e in mat.entries()])
            cMap = pg.plt.cm.get_cmap("Set3", len(uniqueIDs))
            sm = pg.plt.cm.ScalarMappable(cmap=cMap)

            cBar = createColorBar(sm, ax=ax, label="Matrix ID",
                                  cMin=-0.5, cMax=len(uniqueIDs)-0.5)
            ticks = np.arange(len(uniqueIDs))
            cBar.set_ticks(ticks)

            labels = []
            for ID in uniqueIDs:
                label = "{:d}".format(ID)
                labels.append(label)
            cBar.set_ticklabels(labels)

    else:
        pg.error("Matrix type not supported yet.")
    return ax, cBar
