# -*- coding: utf-8 -*-
"""
Created on Thu Feb 02 14:15:13 2012

@author: Guenther.T
"""

import pygimli as pg
from pygimli.utils import gmat2numpy
import numpy as np
from scipy.linalg import inv

def iterateBounds( inv, dchi2 = 0.5, maxiter = 100, change = 1.02 ):
    '''
    iterateBounds - find parameter bounds by iterating model parameter 
    until error bound is reached
    '''
    f = inv.forwardOperator()

    model = inv.model()
    resp  = inv.response()

    nd, nm = len( resp ), len( model )
    modelU = np.zeros( nm )
    modelL = np.zeros( nm )
    maxchi2 = inv.chi2() + dchi2

    for im in range( nm ):
        model1 = pg.RVector( model )
        chi2, iter = 0., 0
        while ( chi2 < maxchi2 ) & ( iter < maxiter ) :
            iter += 1
            model1[ im ] *= change
            resp1 = f( model1 )
            chi2 = inv.getPhiD( resp1 ) / nd

        #if iter < maxiter - 1:
        modelU[ im ] = model1[ im ]

        model2 = pg.RVector( model )
        chi2, iter = 0., 0
        while ( chi2 < maxchi2 ) & ( iter < maxiter ) :
            iter += 1
            model2[ im ] /= change
            resp2 = f( model2 )
            chi2 = inv.getPhiD( resp2 ) / nd

        #if iter < maxiter - 1:
        modelL[ im ] = model2[ im ]

    return modelL, modelU

def modCovar( inv ):
    td = np.asarray( inv.transData().deriv( inv.response() ) )
    tm = np.asarray( inv.transModel().deriv( inv.model() ) )
    J = td.reshape(len(td),1) * gmat2numpy( inv.forwardOperator().jacobian() ) * (1./tm)
    d = 1. / np.asarray( inv.transData().error( inv.response(), inv.error() ) )
    DJ = d.reshape(len(d),1) * J
    JTJ = DJ.T.dot( DJ )
    MCM = inv( JTJ )   # model covariance matrix
    varVG = np.sqrt( np.diag( MCM ) ) # standard deviations from main diagonal
    di = ( 1. / varVG )  # variances as column vector
    MCMs = di.reshape(len(di),1) * MCM * di  # scaled model covariance (=correlation) matrix
    return varVG, MCMs

def print1dBlockVar( var, thk, xpos=None ):
    if xpos is None:
        xpos = plt.xlim()[0]
    
    nlay = len(thk)+1
    zl  = np.cumsum(thk)
    zvec = np.hstack((zl,zl-thk/2,zl[-1]+thk[-1]/2))
    for j in range(nlay*2-1):
        v = np.log(1.+var[j])
        if j<nlay-1:
            plt.text( xpos ,zvec[j],'$\delta$='+str(np.round_( v*thk[j], 1 ))+'m')            
        else:
            plt.text( xpos, zvec[j],'$\delta$='+str(np.round_( v*100., 1 ))+'$\%$')

