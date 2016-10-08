# -*- coding: utf-8 -*-
"""
Test for interpolation matrix.
"""
import numpy as np

import unittest

import pygimli as pg


class TestInterpolate(unittest.TestCase):

    def test_Interpolate(self):
        grid = pg.createGrid(x=[0.0, 1.0], y=[0.0, 1.0])
        u = pg.RVector(grid.nodeCount(), 1.)

        # test with pg.interpolate
        queryPos = [0.2, 0.2]
        uI = pg.interpolate(srcMesh=grid,
                            inVec=u, destPos=[queryPos,queryPos] )

        np.testing.assert_allclose(uI[0], 1.)

        # test manual interpolation
        c = grid.findCell(queryPos)
        uI = c.pot(queryPos, u)
        np.testing.assert_allclose(uI, 1.)

        # test with manual interpolationMatrix generation
        I = pg.RSparseMapMatrix(1, grid.nodeCount())
        cI = c.N(c.shape().rst(queryPos))
        for i in range(c.nodeCount()):
            I.addVal(0, c.node(i).id(), cI[i])

        uI = I.mult(u)
        np.testing.assert_allclose(uI[0], 1)

        # test with automatic interpolationMatrix generation
        I = grid.interpolationMatrix([[0.0, 0.0],
                                      [1.0, 0.0],
                                      [1.0, 1.0],
                                      [0.0, 1.0]])
        uI = I * u
        np.testing.assert_allclose(uI, u)




if __name__ == '__main__':
    unittest.main()
