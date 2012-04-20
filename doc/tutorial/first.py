#!/usr/bin/env python

'''

*This introductory sentence should state the intent and goal of the tutorial. Keep it brief.
This next block should state any assumptions that you the writer are making. Present them in list form.*

Let us start with the very simple inverse problem of fitting a polynomial curve of degree :math:`P`

.. math::

    f(x) = p_0 + p_1 x + \ldots + p_P x^P = \sum\limits_{i=0}^{P} p_i x^i
'''

__docformat__ = 'restructuredtext'

import pygimli as g
import numpy as np

class FunctionModelling( g.ModellingBase ):
    def __init__( self, nc, xvec, verbose = False  ):
        g.ModellingBase.__init__( self, verbose )
        self.x_ = xvec
        self.nc_ = nc
        self.regionManager().setParameterCount( nc )
    
    def response( self, par ):
        y = g.RVector( self.x_.size(), par[ 0 ] )
        
        for i in range( 1, self.nc_ ):
            y += g.pow( self.x_, i ) * par[ i ];
        return y;
    
    def startModel( self ): 
        return g.RVector( self.nc_, 0.5 )


#!this is doku within


# evaluate f(x) = 1.1 + 2.1 * x
x = g.asvector( np.arange( 0., 10., 1 ) )

'''
this is doku within
'''

y = 1.1 + 2.1 * x

'''
this is doku within
'''

print x, y

nP = 3
# two coefficients and x-vector (first data column)
fop = FunctionModelling( nP, x )

# initialize inversion with data and forward operator and set options
inv = g.RInversion( y, fop );

# constant absolute error of 0.01 (not necessary, only for chi^2)
inv.setAbsoluteError( 0.01 );

# the problem is well-posed and does not need regularization
inv.setLambda( 0 );

# actual inversion run yielding coefficient model
coeff = inv.run();

import pylab as P
P.plot( x, y, 'rx', x, inv.response(), 'b-' )
P.show()
