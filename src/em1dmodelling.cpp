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

#include "gimli.h"
#include "em1dmodelling.h"
#include "meshgenerators.h"

namespace GIMLI {

RVector MT1dModelling::rhoaphi( const RVector & rho, const RVector & thk ) { // after mtmod.c by R.-U. Börner
    size_t nperiods = periods_.size();
    RVector rhoa( nperiods ), phi( nperiods );
    static double my0 = PI * 4e-7;
    Complex i_unit( 0.0 , 1.0 ), adm, alpha, tanalpha;
    CVector z( nlay_ );
    for ( size_t i = 0 ; i < nperiods ; i++ ) {
        double omega = 2.0 * PI / periods_[ i ];
        z[ nlay_ - 1 ] = sqrt( i_unit * omega * rho[ nlay_ - 1 ] / my0 );
        for ( int k = nlay_ - 2 ; k >= 0 ; k-- ) {
            adm = sqrt( my0 / ( rho[ k ] * i_unit * omega ) );
            alpha = thk[ k ] * sqrt( i_unit * my0 * omega / rho[k] );
            tanalpha = sinh( alpha ) / cosh( alpha );
            z[ k ] = ( adm * z[ k + 1 ] + tanalpha ) / ( adm * z[ k + 1 ] * tanalpha + 1.0 );
            z[ k ] /= adm;
        }
        rhoa[i] = abs( z[ 0 ] ) * abs( z[ 0 ] ) * my0 / omega;
        phi[i] = std::atan( imag( z[ 0 ] ) / real( z[ 0 ] ) );
    }
    return cat( rhoa, phi );
}

RVector MT1dModelling::rhoa( const RVector & model ){ //! app. res. for thk/res vector
    if ( model.size() != nlay_ * 2 - 1 ) return EXIT_VECTOR_SIZE_INVALID;
    RVector thk( model, 0, nlay_ - 1 ), rho( model, nlay_ - 1, 2 * nlay_ - 1 );
    return rhoa( rho, thk );
}

    /*! the actual (full) forward operator returning app.res.+phase for thickness+resistivity */
RVector MT1dModelling::response( const RVector & model ) {
    if ( model.size() != nlay_ * 2 - 1 ) return EXIT_VECTOR_SIZE_INVALID;
    RVector thk( model, 0, nlay_ - 1 ), rho( model, nlay_ - 1, 2 * nlay_ - 1 );
    return rhoaphi( rho, thk );
}

RVector MT1dModelling::rhoa( const RVector & rho, const RVector & thk ) { // after mtmod.c by R.-U. Börner
    size_t nperiods = periods_.size();
    RVector rhoa( nperiods );
    static double my0 = PI * 4e-7;
    Complex i_unit( 0.0 , 1.0 ), adm, alpha, tanalpha;
    CVector z( nlay_ );
    for ( size_t i = 0 ; i < nperiods ; i++ ) {
        double omega = 2.0 * PI / periods_[ i ];
        z[ nlay_ - 1 ] = sqrt( i_unit * omega * rho[ nlay_ - 1 ] / my0 );
        for ( int k = nlay_ - 2 ; k >= 0 ; k-- ) {
            adm = sqrt( my0 / ( rho[ k ] * i_unit * omega ) );
            alpha = thk[ k ] * sqrt( i_unit * my0 * omega / rho[k] );
            tanalpha = sinh( alpha ) / cosh( alpha );
            z[ k ] = ( adm * z[ k + 1 ] + tanalpha ) / ( adm * z[ k + 1 ] * tanalpha + 1.0 );
            z[ k ] /= adm;
        }
        rhoa[i] = abs( z[ 0 ] ) * abs( z[ 0 ] ) * my0 / omega;
    }
    return rhoa;
}

RVector MRSModelling::response( const RVector & model ) { 
    RVector outreal( *KR_ * model );
    RVector outimag( *KI_ * model );
    return RVector( sqrt( outreal * outreal + outimag * outimag ) ); 
}
    
void MRSModelling::createJacobian( RMatrix & jacobian, const RVector & model ) {
    RVector ddr( *KR_ * model );
    RVector ddi( *KI_ * model );
    RVector dda( sqrt( ddr * ddr + ddi * ddi ) );
    
    jacobian.resize( dda.size(), model.size() );
    
    for ( size_t i = 0 ; i < KR_->rows() ; i++ ) {
        jacobian[ i ] = ( (*KR_)[ i ] * ddr[i] + (*KI_)[ i ] * ddi[i] ) / dda[i];         
    }
}

RVector MRS1dBlockModelling::response( const RVector & model ){
        //! extract water content and thickness from model vector
        RVector wc( model, nlay_ - 1 , nlay_ * 2 - 1 );
        RVector thk( model, 0 , nlay_ - 1 );
        //! fill vector of original size wit last layer water content
        RVector wcvec( nvec_, wc[ nlay_ - 1 ] );
        size_t iz1 = 0, iz2 = 0;
        double zthk = 0;
        //! run through layers and fill water content
        for ( size_t i = 0 ; i < nlay_ - 1 ; i++ ){
            zthk += thk[ i ];
            iz2 = 0;
            while ( iz2 < zvec_.size() && zvec_[ iz2 ] < zthk ) iz2++;
            if ( iz2 > nvec_ ) iz2 = nvec_;
            for ( size_t j = iz1 ; j < iz2 ; j++ ) wcvec[ j ] = wc[ i ];
            if ( iz2 + 1 >= zvec_.size() ) break; // end reached
            wcvec[ iz2 ] = ( ( zthk - zvec_[ iz2 ] ) * wc[ i ] 
                           + ( zvec_[ iz2 + 1 ] - zthk ) * wc[ i + 1 ] ) 
                         / ( zvec_[ iz2 + 1 ] - zvec_[ iz2 ] );
            iz1 = iz2 + 1;
        }
        
        if ( verbose_ ) save( wcvec, "wctmp.vec" );
        //! call original forward response and return;
        return MRSModelling::response( wcvec );
    }
    
} // namespace GIMLI{
