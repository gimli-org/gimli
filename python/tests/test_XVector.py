#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as pg

import numpy as np

import traceback
import unittest

class TestRVectorMethods(unittest.TestCase):

    def test_XVectorBasics(self):

        def testVector(v):

            for t in v:
                self.assertEqual(t, 1.0)

            self.assertEqual(sum(v), 5.0)
            self.assertFalse(pg.haveInfNaN(v))

            v[1] = 0
            self.assertEqual(v[1], 0)                         
                
            v[1] = 1
            self.assertEqual(v[1], 1)                         
            #print(v/v)

        testVector(pg.RVector(5, 1.0))
        testVector(pg.CVector(5, 1.0))
        testVector(pg.BVector(5, True))
        testVector(pg.IVector(5, 1))

    def test_XVectorSetVal(self):
        """
        """
        # 
        vec = pg.RVector(5, 1.0)
        self.assertEqual(vec[0], 1.0)
        
        vec[0] += 1.0
        self.assertEqual(vec[0], 2.0)
        
        vec[0] = vec[0] + 1.0
        self.assertEqual(vec[0], 3.0)
        
        vec.setVal(4.0, 0)
        self.assertEqual(vec[0], 4.0)
        
        vec.setVal(5.0, 0, 1)
        self.assertEqual(vec[0], 5.0)
    
        vec.setVal(6.0)
        self.assertEqual(vec[0], 6.0)
        
        vec.setVal(7.0, vec > 0.0)
        self.assertEqual(vec[0], 7.0)

        
    def test_RVectorOP(self):
        v = pg.RVector(5, 1.0)

        self.assertEqual(sum(v + 1), 10)
        self.assertEqual(sum(v - 2), -5)
        self.assertEqual(sum(v * 2), 10)
        self.assertEqual(sum(v / 1), 5)
        self.assertEqual(sum(1 + v), 10)
        self.assertEqual(sum(-1 - v), -10)
        self.assertEqual(sum(2 * v), 10)
        self.assertEqual(sum(1 / v), 5)
        self.assertEqual(sum(v + 1.0), 10)
        self.assertEqual(sum(v - 2.0), -5)
        self.assertEqual(sum(v * 2.0), 10)
        self.assertEqual(sum(v / 1.0), 5)
        self.assertEqual(sum(1.0 + v), 10)
        self.assertEqual(sum(-1.0 - v), -10)
        self.assertEqual(sum(2.0 * v), 10)
        self.assertEqual(sum(1.0 / v), 5)

        #v = pg.CVector(5, 1.0)

        #self.assertEqual(sum(v + 1), 10)
        #self.assertEqual(sum(v - 1), 0)
        #self.assertEqual(sum(v * 1), 5)
        #self.assertEqual(sum(v / 1), 5)
        #self.assertEqual(sum(1 + v), 10)
        #self.assertEqual(sum(1 - v), 0)
        #self.assertEqual(sum(1 * v), 5)
        #self.assertEqual(sum(1 / v), 5)
        #self.assertEqual(sum(v + 1.0), 10)
        #self.assertEqual(sum(v - 1.0), 0)
        #self.assertEqual(sum(v * 1.0), 5)
        #self.assertEqual(sum(v / 1.0), 5)
        #self.assertEqual(sum(1.0 + v), 10)
        #self.assertEqual(sum(1.0 - v), 0)
        #self.assertEqual(sum(1.0 * v), 5)
        #self.assertEqual(sum(1.0 / v), 5)
     
    def test_RVectorIndexRW(self):
        v = pg.RVector(5, 2.0)
        
        self.assertEqual(v[0], 2.0)
        v += 1.0
        self.assertEqual(v[0], 3.0)         
        v += 1
        self.assertEqual(v[0], 4.0)         
        v[1] = 1.0
        self.assertEqual(v[1], 1.0)
        v[1] += 1.0
        self.assertEqual(v[1], 2.0)  
        v[[1,2]] = 2.0
        self.assertEqual(v[1], 2.0)

        v[pg.IVector(1,3)] = 3.0
        self.assertEqual(v[3], 3.0)
        
        v.setVal(1.0, 1)
        self.assertEqual(v[1], 1.0)
         
    def test_RVectorFuncts(self):
        v = pg.RVector(5, 2.0)
        self.assertEqual(sum(pg.pow(v, 2)), 20)
        self.assertEqual(sum(pg.pow(v, 2.0)), 20)
        self.assertEqual(sum(pg.pow(v, v)), 20)
        
    def test_Slices(self):
        a = pg.RVector(np.arange(10.))

        np.testing.assert_array_equal(a[:], np.arange(10.)[:])
        np.testing.assert_array_equal(a[::], np.arange(10.)[::])
        np.testing.assert_array_equal(a[::1], np.arange(10.)[::1])
        np.testing.assert_array_equal(a[::-1], np.arange(10.)[::-1])
        
        np.testing.assert_array_equal(a[0:3:1], np.arange(10.)[0:3:1])
        np.testing.assert_array_equal(a[0:3:2], np.arange(10.)[0:3:2])
        
        np.testing.assert_array_equal(a[3:0:-1], np.arange(10.)[3:0:-1])
        np.testing.assert_array_equal(a[3:0:-2], np.arange(10.)[3:0:-2])
        
        np.testing.assert_array_equal(a[0:3:-1], np.arange(10.)[0:3:-1])
        np.testing.assert_array_equal(a[0:3:-2], np.arange(10.)[0:3:-2])
        

    #def test_isupper(self):
        #self.assertTrue('FOO'.isupper())
        #self.assertFalse('Foo'.isupper())

    #def test_split(self):
        #s = 'hello world'
        #self.assertEqual(s.split(), ['hello', 'world'])
        ## check that s.split fails when the separator is not a string
        #with self.assertRaises(TypeError):
            #s.split(2)

    def test_IndexAccess(self):
        # (double) array/vector
        an = np.arange(10.)
        ag = pg.RVector(an)

        # bn = nd.array(bool)
        bn = (an > 4.)
        self.assertEqual(type(bn), np.ndarray)
        self.assertEqual(bn.dtype, 'bool')
        self.assertEqual(sum(bn), 5)

        # bg = BVector
        bg = (ag > 4.)
        self.assertEqual(type(bg), pg.BVector)
        self.assertEqual(sum(bg), 5)

        # BVector(nd.array(bool))
        self.assertEqual(len(bg), len(pg.BVector(bn)))
        self.assertEqual(sum(bg), sum(pg.BVector(bn)))
        self.assertEqual(bg[0], pg.BVector(bn)[0])
        np.testing.assert_array_equal(bg, pg.BVector(bn))
         
        # In = nd.array(int)
        In = np.nonzero(bn)[0]
        self.assertEqual(type(In), np.ndarray)
        self.assertEqual(In.dtype, 'int')
        self.assertEqual(len(In), 5)
        self.assertEqual(In[0], 5)
        
        # np.nonzero(bg)
        np.testing.assert_array_equal(In, np.nonzero(bg)[0])

        # Ig = IndexArray
        Ig = pg.find(bg)
        self.assertEqual(type(Ig), pg.IndexArray)
        self.assertEqual(len(Ig), 5)
        self.assertEqual(Ig[0], 5)
        
        # pg.find(nd.array(bool))
        np.testing.assert_array_equal(Ig, pg.find(bn))
        
        ## Indexoperators ##
        # ndarray [nd.array(bool)] == ndarray [nd.array(int)]
        np.testing.assert_equal(an[bn], an[In])
        self.assertEqual(len(an[bn]), 5)
        self.assertEqual(an[bn][0], 5)
        
        # ndarray[IndexArray] == ndarray [nd.array(int)]
        np.testing.assert_equal(an[Ig], an[In])
        
        # ndarray[BVector] == ndarray [nd.array(bool)]
        np.testing.assert_array_equal(an[np.array(bg, dtype='bool')], an[bn])
        np.testing.assert_array_equal(an[np.array(bg)], an[bn])
        np.testing.assert_array_equal(an[bg.array()], an[bn])
        ## this fails because it is interpreted as an[[0,0,0,1,1,1]] .. 
        #np.testing.assert_equal(an[bg], an[bn])
        

        # RVector [BVector] == RVector [IndexArray]
        np.testing.assert_array_equal(ag[bg], ag[Ig])
        self.assertEqual(len(ag[bg]), 5)
        self.assertEqual(ag[bg][0], 5)

        # RVector [BVector] ==  RVector [nd.array(bool)]
        np.testing.assert_array_equal(ag[bg], ag[bn])
        np.testing.assert_equal(sum(ag[bg]), sum(ag[bn]))
        
        # RVector [IndexArray] ==  RVector [nd.array(int)]
        np.testing.assert_array_equal(ag[Ig], ag[In])
                
        # RVector(BVector) ==  RVector(nd.array(bool))
        
        # RVector(IndexArray) ==  RVector(nd.array(int))


if __name__ == '__main__':
    
    unittest.main()
   
    #suite = unittest.TestSuite()
    
    ##suite.addTest(TestRVectorMethods("test_RVectorOP"))
    #suite.addTest(TestRVectorMethods("test_IndexAccess"))
    ##suite.addTest(TestRVectorMethods("test_Slices"))
        
    #runner = unittest.TextTestRunner()
    #runner.run(suite)
    
    