#!/usr/bin/env python
"""
"""
import unittest
import time

import numpy as np
import pygimli as pg
import pygimli.meshtools as mt

class TestMod(pg.ModellingBase):
    def __init__(self, mesh, verbose=True):
        super(TestMod, self).__init__(verbose)
        self.meshlist = []
        for i in range(2):
            for cell in mesh.cells():
                cell.setMarker(i + 1)
            self.meshlist.append(pg.Mesh(mesh))
            self.regionManager().addRegion(i + 1, self.meshlist[i])
            self.regionManager().region(i + 1).setConstraintType(1)

        self.mesh = self.meshlist[0]


class TestRM(unittest.TestCase):

    def test_Constraints(self):
        """ Test FOP """

        grid = pg.createGrid(x=[0., 1., 2.], y=[0., 1., 2.])
        fop = TestMod(grid)

        fop.createConstraints()
        #pg.solver.showSparseMatrix(fop.constraints(), full=True)

        fop.regionManager().setConstraintType(2)
        fop.createConstraints()
        #pg.solver.showSparseMatrix(fop.constraints(), full=True)

if __name__ == '__main__':

    unittest.main()
