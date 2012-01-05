/***************************************************************************
 *   Copyright (C) 2006-2012 by the resistivity.net development team       *
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

#ifndef _GIMLI_ELEMENTMATRIX__H
#define _GIMLI_ELEMENTMATRIX__H

#include "gimli.h"
#include "vector.h"
#include "matrix.h"

namespace GIMLI{

DLLEXPORT std::ostream & operator << ( std::ostream & str, const ElementMatrix< double > & pos );

class DLLEXPORT IntegrationRules{
public:
    IntegrationRules();

    inline const std::vector < RVector3 > & gauAbscissa( uint order ) const { return gauAbscissa_[ order ]; }
    inline const RVector & gauWeights( uint order ) const { return gauWeights_[ order ]; }

    inline const std::vector < RVector3 > & edgAbscissa( uint order ) const { return edgAbscissa_[ order ]; }
    inline const RVector & edgWeights( uint order ) const { return edgWeights_[ order ]; }

    inline const std::vector < RVector3 > & triAbscissa( uint order ) const { return triAbscissa_[ order ]; }
    inline const RVector & triWeights( uint order ) const { return triWeights_[ order ]; }

    inline const std::vector < RVector3 > & tetAbscissa( uint order ) const { return tetAbscissa_[ order ]; }
    inline const RVector & tetWeights( uint order ) const { return tetWeights_[ order ]; }

    inline const std::vector < RVector3 > & quaAbscissa( uint order ) const { return quaAbscissa_[ order ]; }
    inline const RVector & quaWeights( uint order ) const { return quaWeights_[ order ]; }

    inline const std::vector < RVector3 > & hexAbscissa( uint order ) const { return hexAbscissa_[ order ]; }
    inline const RVector & hexWeights( uint order ) const { return hexWeights_[ order ]; }

protected:
    void initGau_();
    void initEdg_();
    void initTri_();
    void initTet_();
    void initQua_();
    void initHex_();

    std::vector < std::vector < RVector3 > > gauAbscissa_;
    std::vector < RVector > gauWeights_;
    std::vector < std::vector < RVector3 > > edgAbscissa_;
    std::vector < RVector > edgWeights_;
    std::vector < std::vector < RVector3 > > triAbscissa_;
    std::vector < RVector > triWeights_;
    std::vector < std::vector < RVector3 > > tetAbscissa_;
    std::vector < RVector > tetWeights_;
    std::vector < std::vector < RVector3 > > quaAbscissa_;
    std::vector < RVector > quaWeights_;
    std::vector < std::vector < RVector3 > > hexAbscissa_;
    std::vector < RVector > hexWeights_;
};

template < class T > class DLLEXPORT ElementMatrix {
public:
    ElementMatrix( ) { }

    ElementMatrix( const ElementMatrix < T > & E ) {
        std::cout << "ElementMatrix( const ElementMatrix < T > & E ) " << std::endl;
        THROW_TO_IMPL
//     this->resize( E.size() );
//         for ( uint i = 0; i < E.size(); i ++ ) mat_[ i ] = E.mat( i );
//         idx_ = E.idx();
//         initBaseMatricies();
    }

    ElementMatrix < T > & operator = ( const ElementMatrix < T > & E ) {
        std::cout << "ElementMatrix::operator = (" << std::endl;
        THROW_TO_IMPL
        if ( this != & E ){
//             this->resize( E.size() );
//             for ( uint i = 0; i < E.size(); i ++ ) mat_[ i ] = E.row( i );
//             idx_ = E.idx();
        } return *this;
    }

    ~ElementMatrix() {}

    inline const RVector & operator[]( uint row ) const { return mat_[ row ]; }

    void resize( uint newSize ) {
        idx_.resize( newSize );
        mat_.resize( newSize, newSize );
    }

    ElementMatrix < T > & operator += ( const ElementMatrix < T > & E ){
        for ( uint i = 0; i < size(); i ++ ){ mat_[ i ] += E.row( i ); } return *this;
    }
  /*ElementMatrix < T > & operator += ( T a ){
    for ( uint i = 0; i < size(); i ++ ){ mat_[ i ] += a; } return *this; }
  ElementMatrix < T > & operator *= ( T a ){
    for ( uint i = 0; i < size(); i ++ ){ mat_[ i ] *= a; } return *this; }
  */
// ElementMatrix < T > & operator OP##= ( T val ) {
//    for ( register size_t i = 0; i < size(); i ++ ) mat_[ i ] OP##= val; return *this; }

    #define DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__( OP )                   \
        void operator OP##= ( T val ) { \
            for ( register size_t i = 0; i < size(); i ++ ) mat_[ i ] OP##= val; \
        } \

        DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__( + )
        DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__( - )
        DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__( / )
        DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__( * )

    #undef DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__

    inline uint idx( uint i ) const { return idx_[ i ]; }
    inline uint size() const { return mat_.rows(); }
    inline const T getVal( uint i, uint j ) const { return mat_[ i ][ j ]; }

    inline const RVector & row( uint i ) const { return mat_[ i ]; }
    inline const RMatrix & mat() const { return mat_; }
    inline const std::vector < uint > & idx() const { return idx_; }

    ElementMatrix < T > & u( const MeshEntity & ent );

    ElementMatrix < T > & u2( const MeshEntity & ent );

    ElementMatrix < T > & ux2uy2uz2( const Cell & cell );

    ElementMatrix < T > & u( const MeshEntity & ent, const RVector & w, const std::vector < RVector3 > & x, bool verbose = false );
    ElementMatrix < T > & u2( const MeshEntity & ent, const RVector & w, const std::vector < RVector3 > & x, bool verbose = false );
    ElementMatrix < T > & ux2( const MeshEntity & ent, const RVector & w, const std::vector < RVector3 > & x, bool verbose = false );
    ElementMatrix < T > & ux2uy2( const MeshEntity & ent, const RVector & w, const std::vector < RVector3 > & x, bool verbose = false );
    ElementMatrix < T > & ux2uy2uz2( const MeshEntity & ent, const RVector & w, const std::vector < RVector3 > & x, bool verbose = false );

protected:
    RMatrix mat_;
    std::vector < uint > idx_;

    RMatrix functx_;
    RMatrix functy_;
    RMatrix functz_;

    std::map< uint, RVector > uCache_;
    std::map< uint, RMatrix > u2Cache_;
    
    RMatrix dNdr_;
    RMatrix dNds_;
    RMatrix dNdt_;

    IntegrationRules intRules_;
};

} // namespace GIMLI{

#endif // _GIMLI_ELEMENTMATRIX__H
