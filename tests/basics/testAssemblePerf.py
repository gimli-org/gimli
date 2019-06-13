#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pygimli integration function
"""

import pygimli as pg
import numpy as np

def test(N):
    
    x = np.linspace(0, 1, N)

    pg.tic()
    mesh = pg.createGrid(x, x, x)
    print(mesh)
    pg.toc()
    
    A = pg.matrix.SparseMatrix()
    A.fillStiffnessMatrix(mesh)
    pg.toc()

test(20)
#for N in range(10, 1000, 100):
