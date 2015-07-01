#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as pg

import traceback
import unittest

class TestRVectorMethods(unittest.TestCase):

    def test_XVectorBasics(self):
        
        def testVector(v):
            
            for t in v:
                self.assertEqual(t, 1.0)
                            
            self.assertEqual(sum(v), 5.0)
            self.assertFalse(pg.haveInfNaN(v))
                
            #print(v/v)
        
        testVector(pg.RVector(5, 1.0))
        testVector(pg.CVector(5, 1.0))
        testVector(pg.BVector(5, True))
        testVector(pg.IVector(5, 1))

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
        
    #def test_isupper(self):
        #self.assertTrue('FOO'.isupper())
        #self.assertFalse('Foo'.isupper())

    #def test_split(self):
        #s = 'hello world'
        #self.assertEqual(s.split(), ['hello', 'world'])
        ## check that s.split fails when the separator is not a string
        #with self.assertRaises(TypeError):
            #s.split(2)

if __name__ == '__main__':
    
    
    
    unittest.main()
   