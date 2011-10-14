/***************************************************************************
 *   Copyright (C) 2009 by the resistivity.net development team            *
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

#ifndef _GIMLI_DC1DMODELLING__H
#define _GIMLI_DC1DMODELLING__H

#include "gimli.h"
#include "mesh.h"
#include "meshgenerators.h"
#include "modellingbase.h"
#include "vectortemplates.h"

namespace GIMLI{

/*! Syntactic suger for calculating the response of a model file and output to file. Better in bert Misc /Tools?, pygimli??*/
//DLLEXPORT void calculateDC1D( DataContainer & data, const std::string & modelFile, const std::string & outfile );


//! DC (direct current) 1D modelling 
/*! Classical DC 1D forward operator for given resistivities and thicknesses 
    DC1dModelling( nlayers, ab2, mn2, verbose )
    DC1dModelling( nlayers, am, an, bm, bn, verbose )
    DC1dModelling( nlayers, dataContainer, verbose ) */
class DLLEXPORT DC1dModelling : public ModellingBase {
public:
//    /*! normal constructor DC1dModelling( dataContainer, nlayers, verbose ) deactivated (DataContainer) */
//    DC1dModelling( size_t nlayers, DataContainer & data, bool verbose = false );
    
    /*! constructor for classical Schlumberger sounding */
    DC1dModelling( size_t nlayers, RVector & ab2, RVector & mn2, bool verbose = false );

    /*! general constructor using AM, AN, BM, BN distances */
    DC1dModelling( size_t nlayers, RVector & am, RVector & an, RVector & bm, RVector & bn, bool verbose = false );
    
    /*! constructor using a data container */
//    DC1dModelling( size_t nlayers, DataContainer & data, bool verbose = false );

    virtual ~DC1dModelling() { }

    /*! Returns an RVector of the 1dc response for model = [thickness_ 0, ..., thickness_(n-1), rho_0 .. rho_n]. For n = nlayers. */
    RVector response( const RVector & model );
    
    RVector rhoa( const RVector & rho, const RVector & thk );
    
    RVector kern1d( const RVector & lam, const RVector & rho, const RVector & h );
    
    RVector pot1d( const RVector & R, const RVector & rho, const RVector & thk );

    template < class Vec > Vec rhoaT( const Vec & rho, const RVector & thk ){
        Vec tmp;
        tmp  = pot1dT<Vec>( am_, rho, thk );
        tmp -= pot1dT<Vec>( an_, rho, thk );
        tmp -= pot1dT<Vec>( bm_, rho, thk );
        tmp += pot1dT<Vec>( bn_, rho, thk );
        return tmp * k_ + rho[ 0 ];    
    }
    template < class Vec > Vec kern1dT( const RVector & lam, const Vec & rho, const RVector & h ){
        size_t nr = rho.size();
        size_t nl = lam.size();
        Vec z( nl, rho[ nr - 1 ] );
        Vec p( nl );
        RVector th( nl );
        for ( int i = nr - 2; i >= 0; i-- ) {
            p = ( z - rho[ i ] ) / ( z + rho[ i ] );
            th = tanh( lam * h[ i ] );
            z = ( z + toComplex( th ) * rho[ i ] ) / ( z * th + rho[ i ] ) * rho[ i ];
        }

        Vec ehl( p * RVector( exp( -2.0 * lam * h[ 0 ] ) ) );
        return ehl / ( 1.0 - ehl ) * rho[ 0 ] / 2.0 / PI ;    
    }
    template < class Vec > Vec pot1dT( const RVector & R, const Vec & rho, const RVector & thk ){
        Vec z0( R.size() );
        //double rabs;
        RVector rabs( abs(R) );
        for ( size_t i = 0; i < R.size(); i++ ) {
            //rabs = std::fabs( R[ i ] );
            z0[ i ] = sum( myw_ * kern1dT<Vec>( myx_ / rabs[i], rho, thk ) * 2.0 ) / rabs[i];
        }
        return z0;
    }

    RVector createDefaultStartModel( ) {
        RVector mod( nlayers_ * 2 - 1, meanrhoa_ );
        for ( size_t i = 0; i < nlayers_ -1; i++ ) mod[ i ] = std::pow( 2.0, 1.0 + i );
        return mod;
    }
    

protected:

    /*! init myw and myx */
    void init_();

    void postprocess_();

    size_t nlayers_;
    double meanrhoa_;
    RVector am_;
    RVector an_;
    RVector bm_;
    RVector bn_;
    RVector k_;
    RVector tmp_;

    RVector myx_;
    RVector myw_;
};

/*! DC (direct current) 1D modelling for complex resistivity */
class DLLEXPORT DC1dModellingC : public DC1dModelling {
public:    
//    /*! normal constructor DC1dModelling( dataContainer, nlayers, verbose ) */
//    DC1dModellingC( size_t nlayers, DataContainer & data, bool verbose = false ) :
//        DC1dModelling( nlayers, data, verbose ){}

    /*! constructor for classical Schlumberger sounding */
    DC1dModellingC( size_t nlayers, RVector & ab2, RVector & mn2, bool verbose = false ) :
        DC1dModelling( nlayers, ab2, mn2, verbose ){
            setMesh( createMesh1DBlock( nlayers, 2 ) );
        }

    /*! general constructor using AM, AN, BM, BN distances */
    DC1dModellingC( size_t nlayers, RVector & am, RVector & an, RVector & bm, RVector & bn, bool verbose = false ) :
        DC1dModelling( nlayers, am, an, bm, bn, verbose ){}

    virtual ~DC1dModellingC() { }

    RVector response( const RVector & model );
};

/*! DC1dRhoModelling - Variant of DC 1D modelling with fixed parameterization (occam inversion) */
/*! DC1dRhoModelling( mesh, dataContainer, thicknesses, verbose ) */
class DLLEXPORT DC1dRhoModelling : public DC1dModelling {
public:
//    DC1dRhoModelling( Mesh & mesh, DataContainer & data, RVector & thk, bool verbose = false )
//            : DC1dModelling( mesh, data, thk.size(), verbose ), thk_( thk ) {}
//    DC1dRhoModelling( DataContainer & data, RVector & thk, bool verbose = false )
//            : DC1dModelling( data, thk.size(), verbose ), thk_( thk ) { 
//                mesh_ = createMesh1D( thk.size() + 1, 1 );
//                setMesh( mesh_ );
//            }
    DC1dRhoModelling( RVector & thk, RVector & am, RVector & an, RVector & bm, RVector & bn, bool verbose = false )
            : DC1dModelling( thk.size(), am, bm, am, an, verbose ), thk_( thk ) {
                setMesh( createMesh1D( thk.size() + 1, 1 ) );
            }
    DC1dRhoModelling( RVector & thk, RVector & ab2, RVector & mn2, bool verbose = false )
            : DC1dModelling( thk.size(), ab2, mn2, verbose ), thk_( thk ) { 
                setMesh( createMesh1D( thk.size() + 1, 1 ) );
            }
    virtual ~DC1dRhoModelling() { }

    RVector response( const RVector & rho ) {  return rhoa( rho, thk_ ); }
    
    RVector createDefaultStartModel( ) {  return RVector( thk_.size() + 1, meanrhoa_ ); }
    
protected:
    RVector thk_;
};

} // namespace GIMLI{

#endif // _GIMLI_DC1DMODELLING__H
