#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as pg


def testVector(v):
    print(v)
    for t in v:
        print(t)
    print(sum(v))
    
testVector(pg.RVector(5, 1.0))
testVector(pg.CVector(5, 1.0))
testVector(pg.BVector(5, True))
testVector(pg.IVector(5, 1))

