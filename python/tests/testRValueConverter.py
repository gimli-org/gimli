#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pygimli as pg

import numpy as np
import sys
import gc # for reference counter

def testRVector():
    a=pg.RVector(10)

def testRValSeqRVector3():
    x = [0.0, 1.0, 0.0]
    p = pg.RVector3(x)
    print(p, p.dist(x))

def testRValSequenz():
    a=pg.RVector([1.0, 2.0, 3.0, 4.0])
    print(a)

def testRVecIter():
    a = pg.RVector(10, 1.1)
    print(sum(a))
    
def testFunct(a):
    if a[0] != 1.1:
        print('#'*40  + ' Fail: check __RVectorArrayCall__')
    
    print("testFunct", a)
    
def testNumpyFromRVec():
    a = pg.RVector(10, 1.1)
    print(type(a), sys.getrefcount(a), a)
    a = np.asarray(a)
    print(type(a), sys.getrefcount(a), a)
    testFunct(a)
    print(type(a), sys.getrefcount(a), a)
    
    a = pg.RVector(10, 1.1)
    a = np.array(a)
    testFunct(a)
    
    
    print(sys.getrefcount(a))
    print(a)
    x = np.array(a)
    print(type(x))
    print(x)
    print(sys.getrefcount(x))
    print(sys.getrefcount(a))
    
def testNumpyFromRVec3():
    a = pg.RVector3()
    print(sys.getrefcount(a))
    print(a.array())
    x = np.array(a)
    print(type(x))
    print(x)
    print(sys.getrefcount(x))
    print(sys.getrefcount(a))
    
def testRValNumpyArray():
    x = np.arange(0, 1., 0.2 )
    print(sys.getrefcount(x))
    
    a=pg.RVector(x)
    print(a)
    # should return 2 (self & counter) since the counter is not increased while conversion
    print(sys.getrefcount(x))
    
    x = np.arange(0, 1., 0.2, dtype=np.float64)
    a=pg.RVector(x)
    print('pg.RVector(x):', a)
    
    x = np.array(a)
    a=pg.RVector(x)
    print('pg.RVector(array):', a)
    
    a=pg.RVector([0.2, 0.3, 0.4, 0.5, 0.6])
    print('pg.RVector(list[float]):', a)
    
    a=pg.RVector((0.2, 0.3, 0.4, 0.5, 0.6))
    print('pg.RVector(tuple(float,)):', a)
    
    a=pg.RVector(np.arange(0, 1., 0.2 ))
    print(a)
    a=pg.RVector(np.arange(10.))
    print(a)
    print(pg.norm(np.arange(10.)))
    
def testSlices():
    a=pg.RVector(np.arange(10.))
    print(a[0:3:1], np.arange(10.)[0:3:1])
    print(np.arange(10.)[0:3:-1])
    try:
        print(a[0:3:-1])
    except:
        print("oki")
        pass
    print(a[3:0:-2])
    print(pg.norm(a[0:3:1]))


if __name__ == '__main__':
    print("start tests")
    #testRVector()
    #testRValSeqRVector3()
    #testRValSequenz()
    #testNumpyFromRVec()
    #testNumpyFromRVec3()
    testRValNumpyArray()
    #testSlices()
    