# -*- coding: utf-8 -*-
"""
Created on Thu Feb 02 14:15:13 2012

@author: Guenther.T
"""

import pygimli as g
import numpy as N

def iterateBounds( inv, dchi2 = 0.5, maxiter = 100, change = 1.02 ):
    '''
    iterateBounds - find parameter bounds by iterating model parameter 
    until error bound is reached
    '''
    f = inv.forwardOperator()

    model = inv.model()
    resp  = inv.response()

    nd, nm = len( resp ), len( model )
    modelU = N.zeros( nm )
    modelL = N.zeros( nm )
    maxchi2 = inv.chi2() + dchi2

    for im in range( nm ):
        model1 = g.RVector( model )
        chi2, iter = 0., 0
        while ( chi2 < maxchi2 ) & ( iter < maxiter ) :
            iter += 1
            model1[ im ] *= change
            resp1 = f( model1 )
            chi2 = inv.getPhiD( resp1 ) / nd

        #if iter < maxiter - 1:
        modelU[ im ] = model1[ im ]

        model2 = g.RVector( model )
        chi2, iter = 0., 0
        while ( chi2 < maxchi2 ) & ( iter < maxiter ) :
            iter += 1
            model2[ im ] /= change
            resp2 = f( model2 )
            chi2 = inv.getPhiD( resp2 ) / nd

        #if iter < maxiter - 1:
        modelL[ im ] = model2[ im ]

    return modelL, modelU