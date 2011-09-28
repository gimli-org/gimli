/***************************************************************************
 * This file belongs to GIMLi (Geophysical Inversion and Modelling) library*
 * and was created for documentation purposes only.                        *
 * The underlying application can be found in gimli/apps/fit/polyfit.cpp   *
 * POLYFIT0 is the easiest version of polynomial curve fitting.            *
 * It sets up the modelling class, reads data and does the inversion.      *
 * Run polyfit0 with the synthetic data using: polyfit0 -n 1 y_2.1x+1.1.dat*
 ***************************************************************************/

#include <optionmap.h>
#include <inversion.h>
#include <string>
using namespace GIMLI;

/*! new modelling operator returning f(x) derived from modelling base class */
class FunctionModelling : public ModellingBase {
public:
    /*! constructor, nc: number of coefficients, xvec: abscissa, */
    FunctionModelling( size_t nc, const RVector & xvec, bool verbose = false  )
        : ModellingBase( verbose ), x_( xvec ), nc_( nc ){ 
        this->regionManager().setParameterCount( nc ); 
    }

    /*! the main thing - the forward operator: return f(x) */
    RVector response( const RVector & par ){
        RVector y( x_.size(), par[ 0 ] );
        for ( size_t i = 1; i < nc_; i ++ ) y += pow( x_, i ) * par[ i ];
        return y;
    }

    /*! define the startmodel */
    RVector startModel( ){ return RVector( nc_, 0.5 ); }

protected:
    RVector x_; //! abscissa vector x
    size_t nc_;
};

int main( int argc, char *argv [] ){
    /*! parse command line: data file and number of coefficents */
    std::string datafile;
    int np = 1;
    /*! Parse option map using longoption map */
    OptionMap oMap;
    oMap.setDescription("Curvefit - fits data in datafile with different curves");
    oMap.addLastArg( datafile, "Datafile" );
    oMap.add( np, "n:", "np", "Number of polynomials" );
    oMap.parse( argc, argv );

    /*! load two-column matrix from file (input argument) holding x and y */
    RMatrix xy; loadMatrixCol( xy, datafile );

    /*! initialize modelling operator */
    FunctionModelling f( np + 1, xy[ 0 ] ); //! two coefficients and x-vector (first data column)

    /*! initialize inversion with data and forward operator and set options */
    RInversion inv( xy[ 1 ], f );

    /*! constant absolute error of 0.01 (not necessary, only for chi^2) */
    inv.setAbsoluteError( 0.01 );

     /*! the problem is well-posed and does not need regularization */
    inv.setLambda( 0 );

    /*! actual inversion run yielding coefficient model */
    RVector coeff( inv.run() );

    /*! get actual response and write to file */
    save( inv.response(), "resp.out" );

    /*! print result to screen and save coefficient vector to file */
    std::cout << "y = " << coeff[ 0 ];
    for( int i = 1 ; i <= np ; i++ ) std::cout << " + " << coeff[ i ] << " x^" << i;
    std::cout << std::endl;
    save( coeff, "out.vec" );
    /*! exit programm legally */
    return EXIT_SUCCESS;
}
