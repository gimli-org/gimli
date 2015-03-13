#! /usr/bin/python
# -*- coding: utf-8 -*-

import sys

sys.path.append('../')

import pygimli as pg

print pg.__version__
print "Lin solver autodetection chose: %s " % pg.LinSolver().solverName()
print "Number of virtual cpu (c++): ", pg.numberOfCPU()
print "Number of virtual cpu (python): "
