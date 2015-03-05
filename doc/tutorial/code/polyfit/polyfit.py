# -*- coding: utf-8 -*-
    # -*- coding: iso-8859-1 -*-

import sys
import pygimli as g
import pylab as P

class FunctionModelling( g.ModellingBase ):
    """
        new modelling operator returning f(x) derived from modelling base class
    """
    def __init__( self, nc, xvec, verbose = False  ):
        """
            constructor, nc: number of coefficients, xvec: abscissa, */
        """
        g.ModellingBase.__init__( self, verbose )
        self.x_ = xvec
        self.nc_ = nc
        self.regionManager().setParameterCount( nc )

    def response( self, par ):
        """
           the main thing - the forward operator: return f(x)
        """
        y = g.RVector( self.x_.size(), par[ 0 ] )
        for i in range( 1, self.nc_ ):
            y += g.pow( self.x_, i ) * par[ i ];
        return y;

    def startModel( self ):
        """
            define the startmodel
        """
        return g.RVector( self.nc_, 0.5 )


def main( argv ):
    from optparse import OptionParser

    parser = OptionParser( "Curvefit - fits data in datafile with different curves\n usage: %prog [options] Datafile" )
    parser.add_option("-v", "--verbose", dest = "verbose", action = "store_true"
                            , help = "be verbose", default = False )
    parser.add_option("-n", "--np", dest = "np", type = "int", default = 1
                            , help = "Number of polynomials" )

    (options, args) = parser.parse_args()

    if len( args ) == 0:
        parser.print_help()
        print("Please add a datafile.")
        sys.exit( 2 )
    else:
        datafile = args[ 0 ];

    xy = g.RMatrix()
    g.loadMatrixCol( xy, datafile );

    if options.verbose:
        print("data:", xy)

    # two coefficients and x-vector (first data column)
    f = FunctionModelling( options.np + 1, xy[ 0 ] )

    # initialize inversion with data and forward operator and set options
    inv = g.RInversion( xy[ 1 ], f, True, True );

    # constant absolute error of 0.01 (not necessary, only for chi^2)
    inv.setAbsoluteError( 0.01 );

    # the problem is well-posed and does not need regularization
    inv.setLambda( 0 );

    # actual inversion run yielding coefficient model
    coeff = inv.run();

    # get actual response and write to file.
    g.save( inv.response(), "resp.out" );

    # print result to screen and save coefficient vector to file
    s = "y = " + str( round( coeff[ 0 ] * 1000 ) /1000 )

    for i in range( 1, options.np+1 ):
        s = s + " + " + str( round( coeff[ i ] *1000 ) / 1000 ) + " x^" + str( i )

    print(s)

    g.save( coeff, "out.vec" );

    P.plot( xy[0], xy[1], 'rx', xy[0], inv.response(), 'b-' )
    P.title(s)
    P.xlabel("x");
    P.ylabel("y");
    P.legend(("measured", "fitted"),loc="upper left");
    P.show()

if __name__ == "__main__":
    main( sys.argv[ 1: ] )
