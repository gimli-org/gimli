#! /usr/bin/python
# -*- coding: utf-8 -*-

import sys

sys.path.append('../')

import pygimli as pg

print(pg.version())
print("Lin solver autodetection chose: %s " % pg.LinSolver().solverName())
print("Number of virtual cpu (c++): ", pg.numberOfCPU())
