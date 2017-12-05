# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 17:20:46 2017

@author: skibbe
"""
import pygimli as pg
import numpy as np


def convertCRSIndex2Map(rowIdx, colPtr):
    """Converts CRS indices to uncompressed indices (row, col)."""
    ii = []
    jj = []
    for i in range(len(colPtr) - 1):
        for j in range(colPtr[i], colPtr[i + 1]):
            ii.append(i)
            jj.append(rowIdx[j])
    return ii, jj


def sparseMatrix2Array(matrix, indices=True, getInCRS=True):
    """Extract indices and value from sparse matrix (SparseMap or CRS)

    Get python Arrays from SparseMatrix or SparseMapMatrix in either CRS
    convention (row index, column Start_End, values) or row index, column
    index, values.

    Parameters
    ----------

    matrix: pg.SparseMapMatrix or pg.SparseMatrix
        Input matrix to be transformed to numpy arrays.

    indices: boolean (True)
        Decides weather the indices of the matrix will be returned or not.

    getInCSR: boolean (True)
        If returned, the indices can have the format of a compressed row
        storage (CSR), the default or uncompressed lists with column and row
        indices.

    Returns
    -------

    vals: numpy.ndarray
        Entries of the matrix.

    indices: list, list
        Optional. Returns additional array with the indices for reconstructing
        the matrix in the defined format.

    """
    io_warn = 'Only working for pygimli.SparseMapMatrix or CSR shaped ' +\
              'pygimli.Sparse matrices. Import type is {}'
    assert isinstance(matrix, pg.SparseMapMatrix) or\
        isinstance(matrix, pg.SparseMatrix), io_warn.format(type(matrix))

    if not isinstance(matrix, pg.SparseMatrix):
        matrix = pg.RSparseMatrix(matrix)

    vals = np.array(matrix.vecVals())
    if indices is True:
        rows = list(matrix.vecRowIdx())
        cols = list(matrix.vecColPtr())
        if getInCRS:
            return list(rows), list(cols), vals
        else:
            rr, cc = convertCRSIndex2Map(rows, cols)
            return rr, cc, vals
    else:
        return vals

if __name__ == '__main__':
    pass
# The End
