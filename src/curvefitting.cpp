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

#include "curvefitting.h"

#include "regionManager.h"
#include "vectortemplates.h"

namespace GIMLI {

HarmonicFunction::HarmonicFunction( const RVector & coeff, double xmin, double xmax )
: xMin_( xmin ), xMax_( xmax ){
    setCoefficients( coeff );
}

HarmonicFunction::~HarmonicFunction( ){
}

void HarmonicFunction::setCoefficients( const RVector & coeff ){

    nHarmonic_ = coeff.size() / 2;
    if ( ( (double)coeff.size() / 2.0 - nHarmonic_ ) > TOLERANCE ){
        throwError( 1, WHERE_AM_I + " coefficients size is uneven" + str( coeff.size() ) );
    }
    coeff_ = coeff;
}

double HarmonicFunction::getValue( const double & arg ) const {
    double ret = coeff_[ 0 ];

    double tOne = ( arg - xMin_ ) / ( xMax_ - xMin_ );

    ret += tOne * coeff_[ 1 ];

    for ( size_t j = 1 ; j < nHarmonic_; j++ ){
        ret += ::cos( tOne * PI2 * j ) * coeff_[ j * 2 ];
        ret += ::sin( tOne * PI2 * j ) * coeff_[ j * 2 + 1 ];
    }

    return ret;
}

RVector HarmonicFunction::getValue( const RVector & arg ) const {
    RVector ret( arg.size(), coeff_[ 0 ] );

    RVector tOne( ( arg - xMin_ ) / ( xMax_ - xMin_ ) );

    ret += tOne * coeff_[ 1 ];

    for ( size_t j = 1 ; j < nHarmonic_; j++ ){
        ret += cos( tOne * PI2 * j ) * coeff_[ j * 2 ];
        ret += sin( tOne * PI2 * j ) * coeff_[ j * 2 + 1 ];
    }

    return ret;
}

void HarmonicFunction::copy_( const HarmonicFunction & funct ){
    xMin_ = funct.xMin();
    xMax_ = funct.xMax();
    this->setCoefficients( funct.coefficients() );

}


HarmonicModelling::HarmonicModelling( size_t nh, const RVector & tvec, bool verbose )
: ModellingBase( verbose ),
    t_( tvec ), tMin_( min( tvec ) ), tMax_( max( tvec ) ), nh_( nh ), np_( 2 * nh + 2 ) {

    regionManager_->setParameterCount( np_ );
    A_.clear();
    nt_ = tvec.size();

    //! constant vector of 1 -- offset
    RVector one( nt_, 1.0 );
    A_.push_back( one ); //! const

    //! vector linearly ascending from 0 (tmin) to 1 (tmax) -- drift
    double tMin = min( tvec ), tMax = max( tvec );
    RVector tOne( ( t_ - tMin ) / ( tMax - tMin ) ); //** wieso nicht so:
//    RVector tOne( ( t_ - tMin_ ) / ( tMax_ - tMin_ ) );

    A_.push_back( tOne );

    //! harmonic functions cos/sin( n pi t )
    for ( size_t j = 1 ; j <= nh_ ; j++ ){
        one = cos( tOne * PI2 * j );
        A_.push_back( one );
        one = sin( tOne * PI2 * j );
        A_.push_back( one );
    }
}

RVector HarmonicModelling::response( const RVector & par ){
    return transMult ( A_ , par );
}

RVector HarmonicModelling::response( const RVector & par, const RVector tvec ){
    RVector ret( tvec.size(), par[ 0 ] );

    RVector tOne( ( tvec - tMin_ ) / ( tMax_ - tMin_ ) );

    ret += tOne * par[ 1 ];

    for ( size_t j = 1 ; j <= nh_ ; j++ ){
        ret += cos( tOne * PI2 * j ) * par[ j * 2 ];
        ret += sin( tOne * PI2 * j ) * par[ j * 2 + 1 ];
    }
    return ret;
}

void HarmonicModelling::createJacobian( RMatrix & jacobian, const RVector & model ) {
    //!! jacobian = transpose( A );
    if ( jacobian.rows() != nt_ || jacobian.cols() != np_ ) {
        jacobian.resize( nt_, np_ );

        for ( size_t i = 0 ; i < np_ ; i++ ){
            for ( size_t j = 0 ; j < nt_ ; j++ ){
                jacobian[ j ][ i ] = A_[ i ][ j ];
            }
        }
    }
}

} // namespace GIMLI{
