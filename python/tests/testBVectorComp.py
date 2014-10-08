#! /usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pygimli as pg

# (double) array/vector
an = np.arange(10.)
ag = pg.RVector(an)

# boolean array/vector
bn = (an > 4.)
bg = (ag > 4.)

# index vectors from bools
fn = np.nonzero(bn)[0]
fg = pg.find(bg)

# pure numpy indexing
print(an[bn])
print(an[fn])

# pure pygimli indexing
print(ag(bg))
print(ag(fg))

# work
print(ag[bg])
print(ag[fg])
print(ag[fn])
print(ag[bn])

#%% mixed access
print(ag(bn)) # would be nice ........... TODO need entry in custom_rvalue.cpp
print(an[bg]) # could be possible as well .......... TODO need entry in hand_made_wrappers.py

#%% mixed indexing
print(np.nonzero(bg)) #  yields empty result
print(pg.find(bn)) # breaks

