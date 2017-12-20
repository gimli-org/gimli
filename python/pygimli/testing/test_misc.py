#!/usr/bin/env python
# -*- coding: utf-8 -*-

# write a correct test!
import unittest

import pygimli as pg
import numpy as np

class TestMisc(unittest.TestCase):

    def test_Trans(self):
        """
        """
        f = pg.RTrans()

        x = pg.RVector(3, 1.0)

        np.testing.assert_array_equal(f(x), x)
        np.testing.assert_array_equal(f.inv(x), x)
        np.testing.assert_array_equal(f.inv(f(x)), x)
        self.assertEqual(f.trans(x=1.0), 1.0)
        self.assertEqual(f(1.0), 1.0)
        self.assertEqual(f.inv(1.0), 1.0)

        f = pg.RTransLin(factor=2., offset=4.)
        np.testing.assert_array_equal(f(x), x*2. + 4.)
        np.testing.assert_array_equal(f.trans(x), x*2. + 4.)
        np.testing.assert_array_equal(f.inv(f(x)), x)
        self.assertEqual(f(1.0), 6.0)
        self.assertEqual(f.trans(1.0), 6.0)
        self.assertEqual(f.inv(6.0), 1.0)
        self.assertEqual(f.invTrans(6.0), 1.0)

    def test_DataContainerFilter(self):
        """
        """
        data = pg.DataContainer()
        data.resize(101)
        data.markInvalid([0, 100])

        x = np.array([0, 100], dtype="int")
        a = pg.IndexArray(x)
        data.markInvalid(a)

        data.markInvalid(pg.IndexArray(np.array([0, 100], dtype="int")))
        #data.markInvalid(np.array([0, 100], dtype="int")) # fails
        #data.markInvalid(np.array([0, 100], dtype=np.int64)) # fails

    def test_Operators(self):
        t = pg.RVector(10, 1.0)
        self.assertEqual(len(t == 1.0), len(t > 0))
        self.assertEqual(len(t == 1.0), len(t == 1))


if __name__ == '__main__':

    pg.setDebug(0)
    t = TestMisc()
    t.test_Operators()

    #unittest.main()


