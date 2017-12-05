#!/usr/bin/env python
# -*- coding: utf-8 -*-

# write a correct test!
import unittest

import pygimli as pg
import numpy as np

class TestSparseMatrix(unittest.TestCase):

    def test_Convert(self):
        """
        """
        colIds = range(10)
        rowIds = range(10)
        vals = np.ones(10)

        # Construct SparseMap Matrix from python arrays
        A = pg.SparseMapMatrix(colIds, rowIds, vals)

        # Construct SparseMap -> CRS (compressed row storage)
        S = pg.SparseMatrix(A)

        # Construct CRS -> SparseMap
        A2 = pg.SparseMapMatrix(S)

        # all should by identity matrix
        np.testing.assert_equal(A2.getVal(1,1), 1.0)
        np.testing.assert_equal(sum(S * np.ones(S.cols())), S.rows())
        np.testing.assert_equal(sum(A2 * np.ones(A2.cols())), A2.rows())

        MAP1 = pg.SparseMapMatrix(r=3, c=15)
        CSR = pg.SparseMatrix(MAP1)
        MAP2 = pg.SparseMapMatrix(CSR)

        v3 = pg.RVector(3)
        v15 = pg.RVector(15)

        np.testing.assert_equal((MAP1*v15).size(), 3)
        np.testing.assert_equal((MAP1.transMult(v3)).size(), 15)

        np.testing.assert_equal((CSR*v15).size(), 3)
        np.testing.assert_equal((CSR.transMult(v3)).size(), 15)

        np.testing.assert_equal(MAP1.cols(), MAP2.cols())
        np.testing.assert_equal(CSR.cols(), MAP1.cols())
        np.testing.assert_equal(CSR.rows(), MAP1.rows())
        np.testing.assert_equal(MAP1.rows(), MAP2.rows())

        # testing SparseMatrix to Numpy
        csr = pg.SparseMapMatrix(r=4, c=5)
        check_rows = [0, 0, 1, 2, 3]
        check_cols = [0, 1, 2, 3, 4]
        check_csr_rows = [0, 1, 2, 3, 4]
        check_col_s_e = [0, 2, 3, 4, 5, 5]
        check_vals = np.array([1.0, 3, np.pi, 1e-12, -1.12345e13])
        for i in range(len(check_rows)):
            csr.addVal(check_rows[i], check_cols[i], check_vals[i])

        r1, c1, v1 = pg.utils.sparseMatrix2Array(csr)
        np.testing.assert_allclose(r1, check_csr_rows)
        np.testing.assert_allclose(c1, check_col_s_e)
        np.testing.assert_allclose(v1, check_vals)

        r2, c2, v2 = pg.utils.sparseMatrix2Array(pg.SparseMatrix(csr),
                                                 getInCRS=False)
        np.testing.assert_allclose(r2, check_rows)
        np.testing.assert_allclose(c2, check_cols)
        np.testing.assert_allclose(v2, check_vals)


        A1 = pg.SparseMapMatrix(colIds, rowIds, vals)
        A2 = pg.SparseMapMatrix(colIds, rowIds, vals)
        A1 += A2

        sciA1 = pg.utils.sparseMatrix2csr(pg.SparseMatrix(csr))
        sciA2 = pg.utils.sparseMatrix2csr(csr)
        np.testing.assert_equal(len(sciA1.data), csr.size())
        np.testing.assert_equal(sciA1.data, sciA2.data)
        np.testing.assert_equal(sciA1.indices, sciA2.indices)
        np.testing.assert_equal(sciA1.indptr, sciA2.indptr)

        sciA1 = pg.utils.sparseMatrix2coo(pg.SparseMatrix(csr))
        sciA2 = pg.utils.sparseMatrix2coo(csr)
        np.testing.assert_equal(len(sciA1.data), csr.size())
        np.testing.assert_equal(sciA1.data, sciA2.data)
        np.testing.assert_equal(sciA1.row, sciA2.row)
        np.testing.assert_equal(sciA1.col, sciA2.col)



if __name__ == '__main__':
    unittest.main()
