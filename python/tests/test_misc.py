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
        
        
            
if __name__ == '__main__':
    from pygimli.physics.petro import *
    phi = [0.3]
    tFAPhi = transFwdArchiePhi(rFluid=20)
    r1 = tFAPhi.trans(phi)
    r2 = resistivityArchie(rFluid=20.0, porosity=phi, 
    a=1.0, m=2.0, S=1.0, n=2.0)
    (r1-r2 < 1e-12)

    #unittest.main()


