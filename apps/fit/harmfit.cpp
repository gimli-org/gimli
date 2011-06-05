/***************************************************************************
 *   Copyright (C) 2009-2011 by the resistivity.net development team       *
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

#include <longoptions.h>
#include <inversion.h>
#include <string>
using namespace GIMLI;

/*! new modelling operator returning f(x) derived from modelling base class */
class HarmonicModelling : public ModellingBase {
public:
    /*! constructor, nh: number of coefficients, xvec: abscissa, */
    HarmonicModelling( size_t nh, const RVector & tvec, bool verbose = false  )
        : ModellingBase( verbose ), t_( tvec ), nh_( nh ), np_( 2 * nh + 2 ), tMin_( min( tvec ) ), tMax_( max( tvec ) ){
            regionManager_->setParameterCount( np_ );
//            A_.resize( tvec.size(), nh * 2 + 2 );
            A_.clear();
            nt_ = tvec.size();
            RVector one( nt_, 1.0 ); //! constant offset
            A_.push_back( one );
            RVector tOne( ( t_ - tMin_ ) / ( tMax_ - tMin_ ) );
            A_.push_back( tOne );
            for ( size_t j = 1 ; j <= nh_ ; j++ ){
                one = cos( tOne * 2.0 * PI * j );
//                for ( size_t i = 0 ; i < nt_ ; i++ ) one[ i ] = cos( tOne * 2 * PI * j ); //! cosine
                A_.push_back( one );
                one = sin( tOne * 2.0 * PI * j );
//                for ( size_t i = 0 ; i < nt_ ; i++ ) one[ i ] = sin( tOne * 2 * PI * j ); //! sine
                A_.push_back( one );
            }
    }


    /*! the main thing - the forward operator: return f(x) */
    RVector response( const RVector & par ){
        return transMult ( A_ , par );
    }

    /*! an additional forward operator for another time basis */
    RVector response( const RVector & par, const RVector tvec ){
        RVector result( tvec.size(), par[ 0 ] );
        RVector tOne( ( tvec - tMin_ ) / ( tMax_ - tMin_ ) );
        result += tOne * par[ 1 ];
        for ( size_t j = 1 ; j <= nh_ ; j++ ){
            result += cos( tOne * 2 * PI * j ) * par[ j * 2 ];
            result += sin( tOne * 2 * PI * j ) * par[ j * 2 + 1 ];
        }
        return result;
    }

    /*! optional: generation of jacobian matrix, uncomment for default behaviour (brute force) */
    void createJacobian( RMatrix & jacobian, const RVector & model ) {
        if ( jacobian.rows() != nt_ || jacobian.cols() != np_ ) jacobian.resize( nt_, np_ );
        for ( size_t i = 0 ; i < np_ ; i++ )
            for ( size_t j = 0 ; j < nt_ ; j++ )
                jacobian[ j ][ i ] = A_[ i ][ j ];
    }

    /*! define the startmodel */
    RVector startModel( ){ return RVector( np_, 0.0 ); }

protected:
    RVector t_; //! abscissa vector x
    RMatrix A_; //! function matrix
    size_t nh_, nt_, np_;
    double tMin_, tMax_;
};

int main( int argc, char *argv [] ){
    /*! parse command line: data file and number of coefficents */
    std::string datafile = NOT_DEFINED, tfile = NOT_DEFINED;
    int nh = 1;
    double lambda = 0.0;
    bool verbose = false;
    /*! Parse option map using longoption map */
    OptionMap oMap;
    oMap.setDescription("Curvefit - fits data in datafile with different curves");
    oMap.addLastArg( datafile, "Datafile" );
    oMap.add( verbose,"v", "verbose", "Verbose output" );
    oMap.add( nh    , "n:", "nh", "Number of harmonic pairs" );
    oMap.add( lambda, "l:", "lambda", "Regularization" );
    oMap.add( tfile,  "t:", "tfile", "Time file for output" );
    oMap.parse( argc, argv );

    /*! load two-column matrix from file (input argument) holding x and y */
    RMatrix xy; loadMatrixCol( xy, datafile );

    /*! initialize modelling operator */
    HarmonicModelling f( nh, xy[ 0 ] ); //! two coefficients and x-vector (first data column)

    /*! initialize inversion with data and forward operator and set options */
    RInversion inv( xy[ 1 ], f );

    /*! maximum iterations to 1 due to linear problem (not necessary) */
    //inv.setMaxIter( 1 );

    /*! constant absolute error of 0.01 (not necessary, only for chi^2) */
    inv.setAbsoluteError( 0.01 );

     /*! the problem is well-posed and does not need regularization */
    inv.setLambda( lambda );
    inv.setLocalRegularization( true );

    /*! actual inversion run yielding coefficient model */
    RVector coeff( inv.run() );
    save( coeff, "coeff.vec" );

    /*! get actual response and write to file */
    save( inv.response(), "resp.out" );
    if ( tfile != NOT_DEFINED ) {
        RVector newt; load( newt, tfile );
        save( f.response( coeff, newt ), "vout.vec" );
    }

    /*! print result to screen and save coefficient vector to file */
    std::cout << "rms/rrms( in, out ) = " << rms( xy[ 1 ], inv.response() ) << " / " << rrms( xy[ 1 ], inv.response() ) << "%" << std::endl;
    save( coeff, "out.vec" );
    /*! exit programm legally */
    return EXIT_SUCCESS;
}
