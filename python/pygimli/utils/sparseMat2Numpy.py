#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Utility functions to convert pg SparseMatrices from and to numpy objects"""

import numpy as np
import pygimli as pg


def sparseMatrix2csr(A):
    """Convert SparseMatrix to scipy.csr_matrix.

    Compressed Sparse Row matrix, i.e., Compressed Row Storage (CRS)
    
    Parameters
    ----------
    A: pg.SparseMapMatrix | pg.SparseMatrix
        Matrix to convert from.

    Returns
    -------
    mat: scipy.csr_matrix
        Matrix to convert into.
    """
    #optImport(scipy.sparse, requiredFor="toCRS_matrix")
    from scipy.sparse import csr_matrix
    if isinstance(A, pg.SparseMapMatrix):
        C = pg.SparseMatrix(A)
        return csr_matrix((C.vecVals().array(),
                           C.vecRowIdx(),
                           C.vecColPtr()))
    elif isinstance(A, pg.SparseMatrix):
        return csr_matrix((A.vecVals().array(),
                           A.vecRowIdx(),
                           A.vecColPtr()))
    elif isinstance(A, pg.RBlockMatrix):
        M = A.sparseMapMatrix()
        return sparseMatrix2csr(M)

    return csr_matrix(A)


def sparseMatrix2coo(A, rowOffset=0, colOffset=0):
    """Convert SparseMatrix to scipy.coo_matrix.

    Parameters
    ----------
    A: pg.SparseMapMatrix | pg.SparseMatrix
        Matrix to convert from.

    Returns
    -------
    mat: scipy.coo_matrix
        Matrix to convert into.
    """
    from scipy.sparse import coo_matrix

    vals = pg.RVector()
    rows = pg.IndexArray([0])
    cols = pg.IndexArray([0])
    if isinstance(A, pg.SparseMatrix):
        C = pg.RSparseMapMatrix(A)
        C.fillArrays(vals=vals, rows=rows, cols=cols)
        rows += rowOffset
        cols += colOffset
        return coo_matrix((vals, (rows, cols)), shape=(A.rows(), A.cols()))
    elif isinstance(A, pg.SparseMapMatrix):
        A.fillArrays(vals, rows, cols)
        rows += rowOffset
        cols += colOffset

        return coo_matrix((vals, (rows, cols)), shape=(A.rows(), A.cols()))

    return coo_matrix(A)


def convertCRSIndex2Map(rowIdx, colPtr):
    """Converts CRS indices to uncompressed indices (row, col)."""
    ii = []
    jj = []
    for i in range(len(colPtr)-1):
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
        matrix = pg.SparseMatrix(matrix)

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
