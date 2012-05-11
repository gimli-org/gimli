# -*- coding: utf-8 -*-

import math
import numpy as np

import pygimli as g
from pygimli.utils import getIndex, filterIndex

#from numpy import ones, array, sin, cos, pi, linalg, diag

class harmFunctor():
    def __init__( self, A, coeff, xmin, xSpan ):
        self.A_ = A
        self.coeff_ = coeff
        self.xmin_ = xmin
        self.xSpan_ = xSpan
        
    def __call__( self, x ):
        nc = len( self.coeff_ )/2
        A = np.ones( nc * 2 ) * 0.5
        A[ 1 ] = 3. * x;
        A[ 2::2 ] = np.sin( 2.0 * np.pi * np.arange( 1, nc ) * ( x - self.xmin_ ) / self.xSpan_ )
        A[ 3::2 ] = np.cos( 2.0 * np.pi * np.arange( 1, nc ) * ( x - self.xmin_ ) / self.xSpan_ )
        return sum( A * self.coeff_ )
    
def harmfitNative( y, x = None, nc = None, xc = None, err = None ):
    '''
        python based curve-fit by harmonic functions
        yc = harmfitNativ(x,y[,nc,xc,err]);
        y   .. values of a curve to be fitted
        x   .. abscissa, if none [0..len(y)) 
        nc  .. number of coefficients
        xc  .. abscissa to fit on (otherwise equal to x)
        err .. data error
    '''
    y = np.asarray( y )
    
    if x is None:
       x = range( len(y) )
       
    x = np.asarray( x )
    
    if xc is not None:
        xc = np.asarray( xc )
    else:
        xc = x
    
    if err is None:
        err = np.ones(  ( 1, len( x ) ) ) # ones(size(x)); end

    if nc is None or nc == 0:
    #nc=round(length(x)/30); % number of coefficients
        nc = int( len( x ) / 30 )

    xspan = max( x ) - min( x );
    xmi = min( x );

    #A=ones(length(x),nc*2+2)/2; %nc*(sin+cos)+offset+drift
    A = np.ones( ( len( x ), nc * 2 + 2 ) ) / 2.0 #; %nc*(sin+cos)+offset+drift
    #B=ones(length(xc),nc*2+2)/2;
    B = np.ones( ( len( xc ), nc * 2 + 2 ) ) / 2.0;

    A[:,1] = x[:] * 3.0;
    B[:,1] = xc[:] * 3.0;

    for i in range( 1, nc +1):
        #A(:,i*2+1)=sin(2*i*pi*(x-xmi)/xspan);
        #A(:,i*2+2)=cos(2*i*pi*(x-xmi)/xspan);
        #B(:,i*2+1)=sin(2*i*pi*(xc(:)-xmi)/xspan);
        #B(:,i*2+2)=cos(2*i*pi*(xc(:)-xmi)/xspan);

        A[ :, i * 2 + 0 ] = np.sin( 2.0 * i * np.pi * ( x - xmi ) / xspan );
        A[ :, i * 2 + 1 ] = np.cos( 2.0 * i * np.pi * ( x - xmi ) / xspan );
        B[ :, i * 2 + 0 ] = np.sin( 2.0 * i * np.pi * ( xc - xmi ) / xspan );
        B[ :, i * 2 + 1 ] = np.cos( 2.0 * i * np.pi * ( xc - xmi ) / xspan );

    #w = 1.0 / err; w(~isfinite(w))=0;
    w = 1.0 / err; 
    #w[ (~isfinite(w) ]=0;
    
    #coeff=(spdiags(w,0,length(w),length(w))*A)\(y.*w);
    coeff, res, rank, s = np.linalg.lstsq( np.diag( w, 0) * A, (y * w)[0,:] )
    
    return sum( (B * coeff).T ), harmFunctor( A, coeff, xmi, xspan )

def harmfit( y, x = None, error = None, nCoefficients = 42, resample = None
            , window = None, verbose = False, dosave = False
            , lineSearch = True, robust = False, maxiter = 20 ):
    '''
        HARMFIT - GIMLi based curve-fit by harmonic functions
        y .. 1d-array( len(y) ); values to be fitted
        x .. 1d-array( len(y) ); abscissa data spacing. if not given: 1 * [ 0 .. len( y ) ) 
        error .. 1d-array( len(y) );  data error of y (fit y into the range of error). if not given (absolute err = 1%)
        nCoefficients .. int; Number of harmonic coefficients
        resample .. 1d-array( len( resample ) ); resample y based on fitting coeffients
        window .. just fit data inside window bounds
        return response, coefficients
        
        response .. 1d-array( len( y ) ) if no resample given, else 1d-array( len( resample ) )
        coefficients .. coefficients for harmic functions that fit y
    '''    
    if x is None:
        x = g.asvector( np.arange( len( y ) ) )
    else:
        if not isinstance( x, g.RVector ):
            x = g.asvector( x )

    if not isinstance( y, g.RVector ):
        y = g.asvector( y )
        
    xToFit = None
    yToFit = None
    
    if window is not None:
        idx = g.find( ( x >= window[ 0 ] ) & ( x < window[ 1 ] ) )
        #idx = getIndex( x , lambda v: v > window[0] and v < window[1] )
        
        xToFit = x( idx );
        yToFit = y( idx );
        
        if error is not None: error = error( idx );
    else:
        xToFit = x
        yToFit = y
            
#    print xToFit
#    print yToFit
    fop = g.HarmonicModelling( nCoefficients, xToFit, verbose )
    inv = g.RInversion( yToFit, fop, verbose, dosave );
    if error is not None:
        if not isinstance( error, g.RVector ):
            error = g.asvector( error )
        
        inv.setRelativeError( error )
    else:
        inv.setAbsoluteError( 0.01 )
        
    inv.setLambda( 000.0 )
    inv.setLocalRegularization( True )
    inv.setMaxIter( maxiter )
    inv.setLineSearch( lineSearch )
    inv.setRobustData( robust );
    #inv.setConstraintType( 0 );
    
    coeff = inv.run();
    
    if resample is not None:
        if not isinstance( resample, g.RVector ):
            resample = g.asvector( resample )
            
        ret = fop.response( coeff, resample )
        
#        print ret
        
        if window is not None:
#            print resample
#            print window[0], window[1]
#            print g.find( ( resample < window[ 0 ] ) | ( resample >= window[ 1 ] ) )
            ret.setVal( 0.0, g.find( ( resample < window[ 0 ] ) | ( resample >= window[ 1 ] ) ) )
#            idx = getIndex( resample, lambda v: v <= window[0] or v >= window[1] )
#            for i in idx: ret[ i ] = 0.0
#        print ret
        #sys.exit
        return ret, coeff
    else:
        return inv.response(), coeff
# def harmfit( ... )