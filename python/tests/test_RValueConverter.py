#!/usr/bin/env python
# -*- coding: utf-8 -*-
import unittest
import numpy as np
import pygimli as pg

class TestRVectorMethods(unittest.TestCase):

    def test_RVector(self):
        a = pg.RVector(10)
        self.assertEqual(a.size(), 10.0)
        self.assertEqual(sum(a), 0.0)

    def test_ListToRVector3(self):
        '''
            custom_rvalue.cpp
        '''
        x = [0.0, 1.0, 0.0]
        p = pg.RVector3(x)
        self.assertEqual(p.dist(x), 0.0)
        self.assertEqual(p.dist([1.0, 1.0]), 1.0)

        p = pg.RVector3((0.0, 1.0, 0.0))
        self.assertEqual(p.dist([0.0, 1.0, 0.0]), 0.0)
        
    def test_NumpyToRVector3(self):
        '''
            custom_rvalue.cpp
        '''
        x = np.array([0.0, 1.0, 0.0])
        p = pg.RVector3(x)
        self.assertEqual(p.dist(x), 0.0)
        self.assertEqual(p.dist([1.0, 1.0]), 1.0)

        x = np.array([0.0, 1.0])
        p = pg.RVector3(x)
        self.assertEqual(p.dist([0.0, 1.0, 0.0]), 0.0)

    def test_ListToIndexArray(self):
        '''
            custom_rvalue.cpp
        '''
        idx = [0, 1, 1, 0]

        I = pg.IndexArray(idx)
        self.assertEqual(pg.sum(I), sum(idx))

        bn = (np.array(idx) > 0)  # numpy bool
        idx = np.nonzero(bn)[0]  # numpy int64

        # numyp int64 -> IndexArray
        I = pg.IndexArray(idx)

        self.assertEqual(I.size(), 2)
        self.assertEqual(pg.sum(I), sum(idx))

    def test_ListToRVector(self):
        '''
            custom_rvalue.cpp
        '''
        l = [1.0, 2.0, 3.0, 4.0]
        a = pg.RVector(l)
        self.assertEqual(a.size(), len(l))
        self.assertEqual(pg.sum(a), sum(l))

        l = (0.2, 0.3, 0.4, 0.5, 0.6)
        x = pg.RVector(l)
        self.assertEqual(x.size(), len(l))

    def test_ListToR3Vector(self):
        '''
            custom_rvalue.cpp
        '''
        x = [0.0, 1.0, 0.0]
        p = pg.RVector3(x)
        pl = [p, p, p]
        t = pg.R3Vector(pl)
        self.assertEqual(t.size(), len(pl))

    def test_NumpyToRVector(self):
        '''
            custom_rvalue.cpp
        '''
        x = np.arange(0, 1., 0.2)

        a = pg.RVector(x)
        self.assertEqual(a.size(), len(x))
        self.assertEqual(pg.sum(a), sum(x))

        x = np.arange(0, 1., 0.2, dtype=np.float64)
        a = pg.RVector(x)
        self.assertEqual(a.size(), len(x))
        self.assertEqual(pg.sum(a), sum(x))

    def test_NumpyToIndexArray(self):
        '''
            custom_rvalue.cpp
        '''
        x = np.array(range(10))

        a = pg.IndexArray(x)
        print(a)
        self.assertEqual(a.size(), len(x))
        self.assertEqual(pg.sum(a), sum(x))

        x = np.arange(0, 10, dtype=np.int64)
        a = pg.IndexArray(x)
        print(a)
        self.assertEqual(a.size(), len(x))
        self.assertEqual(pg.sum(a), sum(x))

    def test_RVectorToNumpy(self):
        '''
            implemented through hand_made_wrapper.py
        '''
        # check ob wirklich from array genommen wird!
        v = pg.RVector(10, 1.1)

        a = np.asarray(v)
        self.assertEqual(type(a), np.ndarray)
        self.assertEqual(len(a), 10)

        a = np.array(v)
        self.assertEqual(type(a), np.ndarray)
        self.assertEqual(len(a), 10)

    def test_BVectorToNumpy(self):
        '''
            implemented through hand_made_wrapper.py
        '''
        # check ob wirklich from array genommen wird!
        # wird es noch nicht .. siehe __init__.py:__BVectorArrayCall__
        v = pg.RVector(10, 1.1)
        b = (v == 1.1)

        self.assertEqual(type(b), pg.BVector)

        a = np.asarray(b)
        self.assertEqual(type(a), np.ndarray)
        self.assertEqual(a.dtype, 'bool')
        self.assertEqual(len(a), 10)
        self.assertEqual(sum(a), 10)

        a = np.array(b)
        self.assertEqual(type(a), np.ndarray)
        self.assertEqual(len(a), 10)
        self.assertEqual(sum(a), 10)

    def test_IndexArrayToNumpy(self):
        '''
            implemented through hand_made_wrapper.py
        '''
        #check ob wirklich from array genommen wird!
        # wird es noch nicht .. siehe __init__.py:__BVectorArrayCall__
        v = pg.IndexArray(10, 2)
        self.assertEqual(type(v), pg.IndexArray)
        # print(type(v[0]))
        # print(pg.showSizes())
        a = np.asarray(v)
        self.assertEqual(type(a), np.ndarray)
        #self.assertEqual(a.dtype, 'int64')
        self.assertEqual(len(a), 10)
        self.assertEqual(sum(a), 20)

        a = np.array(v)
        self.assertEqual(type(a), np.ndarray)
        self.assertEqual(len(a), 10)
        self.assertEqual(sum(a), 20)

    def test_RVector3ToNumpy(self):
        '''
            implemented through hand_made_wrapper.py
        '''
        v = pg.RVector3()

        a = np.array(v)
        self.assertEqual(type(a), np.ndarray)
        self.assertEqual(len(a), 3)

    def test_R3VectorToNumpy(self):
        '''
            implemented through hand_made_wrapper.py
        '''
        mesh = pg.createGrid(x=[0, 1, 2], y=[0, 1, 2])

        v = np.asarray(mesh.nodeCenters())

        self.assertEqual(type(v), np.ndarray)
        self.assertEqual(len(v), mesh.nodeCount())

        a = np.array(mesh.cellCenter())
        self.assertEqual(type(a), np.ndarray)
        self.assertEqual(len(a), mesh.cellCount())

if __name__ == '__main__':
    pg.setDebug(0)
    unittest.main()
    # do we need Implicit converter .. currently deactivated in vector.h

    #suite = unittest.TestSuite()
#
   #suite.addTest(TestRVectorMethods("test_ListToR3Vector"))
    #suite.addTest(TestRVectorMethods("test_NumpyToIndexArray"))
#
#    suite.addTest(TestRVectorMethods("test_BVectorToNumpy"))
#    suite.addTest(TestRVectorMethods("test_IndexArrayToNumpy"))
#    suite.addTest(TestRVectorMethods("test_ListToIndexArray"))
#    suite.addTest(TestRVectorMethods("test_ListToRVector"))
#    suite.addTest(TestRVectorMethods("test_NumpyToRVector"))
#
    
    #unittest.TextTestRunner().run(suite)
    
