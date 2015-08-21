/***************************************************************************
 *   Copyright (C) 2009-2012 by the resistivity.net development team       *
 *   Thomas Günther thomas@resistivity.net                                 *
 *   Carsten Rücker carsten@resistivity.net                                *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <optionmap.h>
#include <inversion.h>
#include <string>
using namespace GIMLI;

/*! new modelling operator returning f(x) derived from modelling base class */
class PolynominalModelling : public ModellingBase {
public:
    /*! Constructor, nc: number of coefficients, xvec: abscissa, */
    PolynominalModelling( size_t nc, const RVector & xvec, bool verbose = false  )
        : ModellingBase( verbose ), x_( xvec ), nc_( nc ){ 
        regionManager().setParameterCount( nc ); 
    }

    /*! The main thing - the forward operator: return f(x) */
    RVector response( const RVector & par ){
        RVector y( x_.size(), par[ 0 ] );
        for ( size_t i = 1; i < nc_; i ++ ) y += pow( x_, (double)i ) * par[ i ];
        return y;
    }

    /*! Optional: generation of jacobian matrix, uncomment for default behaviour (brute force) */
    void createJacobian( const RVector & model ) {
        RMatrix *J = dynamic_cast < RMatrix * >( jacobian_ );
        if ( jacobian_->rows() != x_.size() || jacobian_->cols() != nc_ ) J->resize( x_.size(), nc_ );
        
        for ( size_t i = 0 ; i < nc_ ; i++ )
            for ( size_t j = 0 ; j < x_.size() ; j++ )
                (*J)[ j ][ i ] = pow( x_[ j ], (double)i );
    }

    /*! define the startmodel */
    RVector startModel( ){ return RVector( nc_, 0.0 ); }

protected:
    RVector x_; //! abscissa vector x
    size_t nc_;
};

int main( int argc, char *argv [] ){
    /*! parse command line: data file and number of coefficents */
    std::cout << "size_t" << sizeof(size_t) << std::endl;
    std::cout << "ssize_t" << sizeof(ssize_t) << std::endl;
    std::cout << "Index" << sizeof(GIMLI::Index) << std::endl;
    std::cout << "Sindex" << sizeof(GIMLI::SIndex) << std::endl;
   
    std::cout << "8" << sizeof(GIMLI::int8) << std::endl;
    std::cout << "16" << sizeof(GIMLI::int16) << std::endl;
    std::cout << "32" << sizeof(GIMLI::int32) << std::endl;
    std::cout << "64" << sizeof(GIMLI::int64) << std::endl;
    
    std::cout << "8" << sizeof(GIMLI::uint8) << std::endl;
    std::cout << "16" << sizeof(GIMLI::uint16) << std::endl;
    std::cout << "32" << sizeof(GIMLI::uint32) << std::endl;
    std::cout << "64" << sizeof(GIMLI::uint64) << std::endl;
    std::string datafile;
    int np = 1;
    double lambda = 0.001;
    bool verbose = true;
    /*! Parse option map using longoption map */
    OptionMap oMap;
    oMap.setDescription("Fits data in datafile with polynominals.");
    oMap.addLastArg( datafile, "Datafile" );
    oMap.add( np    , "n:", "np", "Number of polynomials" );
    oMap.add( lambda, "l:", "lambda", "Regularization" );
    oMap.parse( argc, argv );

    /*! load two-column matrix from file (input argument) holding x and y */
    RMatrix xy; loadMatrixCol( xy, datafile );

    /*! initialize modelling operator */
    PolynominalModelling f( np + 1, xy[ 0 ] ); //! two coefficients and x-vector (first data column)

    /*! initialize inversion with data and forward operator and set options */
    RInversion inv( xy[ 1 ], f, verbose, false );

    /*! maximum iterations to 1 due to linear problem (not necessary) */
    //inv.setMaxIter( 1 );

    /*! constant absolute error of 0.01 (not necessary, only for chi^2) */
    inv.setAbsoluteError( 0.4 );

     /*! the problem is well-posed and does not need regularization */
    inv.setLambda( lambda );

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
