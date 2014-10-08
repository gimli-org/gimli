#! /usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pygimli as pg

# (double) array/vector
an = np.arange(10.)
ap = pg.RVector(an)
# boolean array/vector
bn = (an>4.)
bp = (ap>4.)
# index vectors from bools
fn = np.nonzero(bn)[0]
fp = pg.find(bp)
# pure numpy indexing
print(an[bn])
print(an[fn])
# pure pygimli indexing
print(ap(bp))
print(ap(fp))
#%% mixed indexing
print(np.nonzero(bp)) #  yields empty result
print(pg.find(bn)) # breaks
#%% mixed access
print(ap(bn)) # would be nice
print(an[bp]) # could be possible as well


