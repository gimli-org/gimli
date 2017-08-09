#!/usr/bin/env python
# encoding: utf-8

r"""
Author: carsten-forty2

Datahandling
============

This tutorial shows ... 
"""

import pygimli as pg
import numpy as np

###############################################################################
#

data = pg.DataContainer()
print(data)
data.setSensorPosition(0, (0.0, 0.0, 0))
print(data)
#define global datasize
data.resize(4)
data.set('a', [1.0, 2.0, 3.0, 4.0])


data.save("dummy.dat")
data = pg.load("dummy.dat")
print(data)

filename = "dummy.dat"
data = pg.load(filename)
print(data)

print(data('a'))


