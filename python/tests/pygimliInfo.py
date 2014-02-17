#! /usr/bin/python
# -*- coding: utf-8 -*-

import sys

sys.path.append( '../' )

import pygimli as g

print g.__version__
print "Lin solver autodetection chose: %s " % g.LinSolver( ).solverName()
print "Number of virtual cpu (c++): ", g.numberOfCPU()
print "Number of virtual cpu (python): "
