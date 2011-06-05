/***************************************************************************
 *   Copyright (C) 2006-2011 by the resistivity.net development team       *
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

#include "elementmatrix.h"

namespace GIMLI{

// static double Triangle6_S1[ 6 ][ 6 ] = {
//     {  3.0,  1.0, 0.0, -4.0,  0.0,  0.0 },
//     {  1.0,  3.0, 0.0, -4.0,  0.0,  0.0 },
//     {  0.0,  0.0, 0.0,  0.0,  0.0,  0.0 },
//     { -4.0, -4.0, 0.0,  8.0,  0.0,  0.0 },
//     {  0.0,  0.0, 0.0,  0.0,  8.0, -8.0 },
//     {  0.0,  0.0, 0.0,  0.0, -8.0,  8.0 }
//   };
//   static double Triangle6_S2[ 6 ][ 6 ] = {
//     {  6.0,  1.0,  1.0, -4.0,  0.0, -4.0 },
//     {  1.0,  0.0, -1.0, -4.0,  4.0,  0.0 },
//     {  1.0, -1.0,  0.0,  0.0,  4.0, -4.0 },
//     { -4.0, -4.0,  0.0,  8.0, -8.0,  8.0 },
//     {  0.0,  4.0,  4.0, -8.0,  8.0, -8.0 },
//     { -4.0,  0.0, -4.0,  8.0, -8.0,  8.0 }
//   };
//   static double Triangle6_S3[ 6 ][ 6 ] = {
//     {  3.0,  0.0,  1.0,  0.0,  0.0, -4.0 },
//     {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },
//     {  1.0,  0.0,  3.0,  0.0,  0.0, -4.0 },
//     {  0.0,  0.0,  0.0,  8.0, -8.0,  0.0 },
//     {  0.0,  0.0,  0.0, -8.0,  8.0,  0.0 },
//     { -4.0,  0.0, -4.0,  0.0,  0.0,  8.0 }
//   };

IntegrationRules::IntegrationRules(){
    initGau_();
    initEdg_();
    initTri_();
    initTet_();
    initQua_();
    initHex_();
}

void IntegrationRules::initGau_(){
    //** 0.Order, n=1 -- just placeholder
    gauAbscissa_.push_back( std::vector < RVector3 > ( 0 ) );
    gauWeights_.push_back( RVector( 0, 0.0 ) );

    //** 1.Order, n=1
    gauAbscissa_.push_back( std::vector < RVector3 > ( 1 ) );
    gauAbscissa_.back()[ 0 ] = RVector3(  0.0,0.0);
    gauWeights_.push_back( RVector( 1, 2.0 ) );

    //** 2.Order, n=2
    gauAbscissa_.push_back( std::vector < RVector3 > ( 2 ) );
    gauAbscissa_.back()[ 0 ] = RVector3( -0.577350269189626,0.0 );
    gauAbscissa_.back()[ 1 ] = RVector3( 0.577350269189626,0.0 );
    gauWeights_.push_back( RVector( 2, 1.0 ) );

    //** 3.Order, n=3
    gauAbscissa_.push_back( std::vector < RVector3 > ( 3 ) );
    gauAbscissa_.back()[ 0 ] = RVector3(  0.0,0.0 );
    gauAbscissa_.back()[ 1 ] = RVector3( -0.774596669241483,0.0 );
    gauAbscissa_.back()[ 2 ] = RVector3(  0.774596669241483,0.0 );
    gauWeights_.push_back( RVector( 3, 0.555555555555556 ) );
    gauWeights_.back()[ 0 ] = 0.888888888888889;

    //** 4.Order, n=4
    gauAbscissa_.push_back( std::vector < RVector3 > ( 4 ) );
    gauAbscissa_.back()[ 0 ] = RVector3(  0.339981043584856,0.0 );
    gauAbscissa_.back()[ 1 ] = RVector3( -0.339981043584856,0.0 );
    gauAbscissa_.back()[ 2 ] = RVector3(  0.861136311594053,0.0 );
    gauAbscissa_.back()[ 3 ] = RVector3( -0.861136311594053,0.0 );
    gauWeights_.push_back( RVector( 4, 0.0 ) );
    gauWeights_.back()[ 0 ] = 0.652145154862546;
    gauWeights_.back()[ 1 ] = 0.652145154862546;
    gauWeights_.back()[ 2 ] = 0.347854845137454;
    gauWeights_.back()[ 3 ] = 0.347854845137454;

    //** 5.Order, n=5
    gauAbscissa_.push_back( std::vector < RVector3 > ( 5 ) );
    gauAbscissa_.back()[ 0 ] = RVector3(  0.000000000000000,0.0 );
    gauAbscissa_.back()[ 1 ] = RVector3( -0.538469310105683,0.0 );
    gauAbscissa_.back()[ 2 ] = RVector3(  0.538469310105683,0.0 );
    gauAbscissa_.back()[ 3 ] = RVector3( -0.906179845938664,0.0 );
    gauAbscissa_.back()[ 4 ] = RVector3(  0.906179845938664,0.0 );
    gauWeights_.push_back( RVector( 5, 0.0 ) );
    gauWeights_.back()[ 0 ] = 0.568888888888889;
    gauWeights_.back()[ 1 ] = 0.478628670499366;
    gauWeights_.back()[ 2 ] = 0.478628670499366;
    gauWeights_.back()[ 3 ] = 0.236926885056189;
    gauWeights_.back()[ 4 ] = 0.236926885056189;

    //** 6.Order, n=6
    gauAbscissa_.push_back( std::vector < RVector3 > ( 6 ) );
    gauAbscissa_.back()[ 0 ] = RVector3(  0.238619186083197,0.0 );
    gauAbscissa_.back()[ 1 ] = RVector3( -0.238619186083197,0.0 );
    gauAbscissa_.back()[ 2 ] = RVector3(  0.661209386466265,0.0 );
    gauAbscissa_.back()[ 3 ] = RVector3( -0.661209386466265,0.0 );
    gauAbscissa_.back()[ 4 ] = RVector3(  0.932469514203152,0.0 );
    gauAbscissa_.back()[ 5 ] = RVector3( -0.932469514203152,0.0 );
    gauWeights_.push_back( RVector( 6, 0.0 ) );
    gauWeights_.back()[ 0 ] = 0.467913934572691;
    gauWeights_.back()[ 1 ] = 0.467913934572691;
    gauWeights_.back()[ 2 ] = 0.360761573048139;
    gauWeights_.back()[ 3 ] = 0.360761573048139;
    gauWeights_.back()[ 4 ] = 0.171324492379170;
    gauWeights_.back()[ 5 ] = 0.171324492379170;
}

void IntegrationRules::initEdg_(){
    /* just transform \int_1^� to \int_0� */
    edgAbscissa_.push_back( std::vector < RVector3 > ( 0 ) );
    edgWeights_.push_back( RVector( 0, 0.0 ) );

    for ( uint i = 1; i < gauAbscissa_.size(); i ++ ){
        edgAbscissa_.push_back( std::vector < RVector3 > ( gauAbscissa_[ i ].size() ) );
        edgWeights_.push_back( 0.5 * gauWeights_[ i ] );

        for ( uint j = 0; j < gauAbscissa_[ i ].size(); j ++ ){
            edgAbscissa_[ i ][ j ] = gauAbscissa_[ i ][ j ] / 2.0 + RVector3( 0.5,0.0 );
        }
    }
}

void IntegrationRules::initTri_(){
    //** 0.Order, n=1, Error: O(h�) -- just placeholder
    triAbscissa_.push_back( std::vector< RVector3 >( 0 ) );
    triWeights_.push_back( RVector( 0, 0.0 ) );

    //** 1.Order, n=1, Error: O(h�)
    triAbscissa_.push_back( std::vector< RVector3 >( 1 ) );
    triAbscissa_.back()[ 0 ] = RVector3( 1.0/3.0, 1.0/3.0  );
    triWeights_.push_back( RVector( 1, 1.0 ) );

    //** 2.Order, n=3, Error: O(h�)
    triAbscissa_.push_back( std::vector< RVector3 >( 3 ) );
    triAbscissa_.back()[ 0 ] = RVector3( 0.5, 0.0 );
    triAbscissa_.back()[ 1 ] = RVector3( 0.5, 0.5 );
    triAbscissa_.back()[ 2 ] = RVector3( 0.0, 0.5 );
    triWeights_.push_back( RVector( 3, 1.0/3.0 ) );

    //** 3.Order, n=4, Error: O(h4)
    triAbscissa_.push_back( std::vector< RVector3 >( 4 ) );
    triAbscissa_.back()[ 0 ] = RVector3( 1.0/3.0, 1.0/3.0 ) ;
    triAbscissa_.back()[ 1 ] = RVector3( 0.2, 0.2 );
    triAbscissa_.back()[ 2 ] = RVector3( 0.6, 0.2 );
    triAbscissa_.back()[ 3 ] = RVector3( 0.2, 0.6 );
    triWeights_.push_back( RVector( 4, 25.0 / 48.0 ) );
    triWeights_.back()[ 0 ] = -27.0/48.0;

    //** 4.Order, n=6, Error: O(h?)
    //**    Joseph E. Flaherty -- Finite Element Analysis CSCI-6860 / MATH-6860
    triAbscissa_.push_back( std::vector< RVector3 >( 6 ) );
    double a = 0.816847572980459, b = 0.091576213509771;
    triAbscissa_.back()[ 0 ] = RVector3( b, b );
    triAbscissa_.back()[ 1 ] = RVector3( a, b );
    triAbscissa_.back()[ 2 ] = RVector3( b, a );
    a = 0.108103018168070, b = 0.445948490915965;
    triAbscissa_.back()[ 3 ] = RVector3( b, b );
    triAbscissa_.back()[ 4 ] = RVector3( a, b );
    triAbscissa_.back()[ 5 ] = RVector3( b, a );
    triWeights_.push_back( RVector( 6, 0.109951743655322 ) );
    triWeights_.back()[ 3 ] = triWeights_.back()[ 4 ] = triWeights_.back()[ 5 ] = 0.223381589678011;

    //** 5.Order, n=7, Error: O(h6)
    triAbscissa_.push_back( std::vector< RVector3 >( 7 ) );
    double sqrt15 = std::sqrt( 15.0 );
    triAbscissa_.back()[ 0 ] = RVector3( 1.0/3.0, 1.0/3.0 );
    b = 2.0 / 7.0 + sqrt15 / 21.0, a = 1.0 - 2.0 * b;
    triAbscissa_.back()[ 1 ] = RVector3( b, b );
    triAbscissa_.back()[ 2 ] = RVector3( a, b );
    triAbscissa_.back()[ 3 ] = RVector3( b, a );
    b = 2.0 / 7.0 + sqrt15 / 21.0, a = 1.0 - 2.0 * b;
    triAbscissa_.back()[ 4 ] = RVector3( b, b );
    triAbscissa_.back()[ 5 ] = RVector3( a, b );
    triAbscissa_.back()[ 6 ] = RVector3( b, a );

    triWeights_.push_back( RVector( 7, 270.0 ) );
    triWeights_.back()[ 1 ] = triWeights_.back()[ 2 ] = triWeights_.back()[ 3 ] = 155.0 + sqrt15;
    triWeights_.back()[ 4 ] = triWeights_.back()[ 5 ] = triWeights_.back()[ 6 ] = 155.0 - sqrt15;
    triWeights_.back() /= 1200.0;
}

void IntegrationRules::initTet_(){
    //**    Joseph E. Flaherty -- Finite Element Analysis CSCI-6860 / MATH-6860
    //** 0.Order, n=1, Error: O(h0) -- just placeholder
    tetAbscissa_.push_back( std::vector< RVector3 >( 0 ) );
    tetWeights_.push_back( RVector( 0, 0.0 ) );

    //** 1.Order, n=1, Error: O(h2)
    tetAbscissa_.push_back( std::vector< RVector3 >( 1 ) );
    tetAbscissa_.back()[ 0 ] = RVector3( 0.25, 0.25, 0.25 );
    tetWeights_.push_back( RVector( 1, 1.0 ) );

    //** 2.Order, n=4, Error: O(h3)
    tetAbscissa_.push_back( std::vector< RVector3 >( 4 ) );
    double a = 0.585410196624969, b = 0.138196601125011;
    tetAbscissa_.back()[ 0 ] = RVector3( b, b, b );
    tetAbscissa_.back()[ 1 ] = RVector3( a, b, b );
    tetAbscissa_.back()[ 2 ] = RVector3( b, a, b );
    tetAbscissa_.back()[ 3 ] = RVector3( b, b, a );
    tetWeights_.push_back( RVector( 4, 0.25 ) );

    //** 3.Order, n=5, Error: O(h4)
    tetAbscissa_.push_back( std::vector< RVector3 >( 5 ) );
    tetAbscissa_.back()[ 0 ] = RVector3( 0.25, 0.25, 0.25 );
    a = 0.5, b = 1.0 / 6.0;
    tetAbscissa_.back()[ 1 ] = RVector3( b, b, b );
    tetAbscissa_.back()[ 2 ] = RVector3( a, b, b );
    tetAbscissa_.back()[ 3 ] = RVector3( b, a, b );
    tetAbscissa_.back()[ 4 ] = RVector3( b, b, a );
    tetWeights_.push_back( RVector( 5, 9.0 / 20.0 ) );
    tetWeights_.back( )[ 0 ] = -4.0 / 5.0;

    //** 4.Order, n=11, Error: O(h?)
    tetAbscissa_.push_back( std::vector< RVector3 >( 11 ) );
    tetAbscissa_.back()[ 0 ] = RVector3( 0.25, 0.25, 0.25 );
    a = 0.785714285714286; b = 0.071428571428571;
    tetAbscissa_.back()[ 1 ] = RVector3( b, b, b );
    tetAbscissa_.back()[ 2 ] = RVector3( a, b, b );
    tetAbscissa_.back()[ 3 ] = RVector3( b, a, b );
    tetAbscissa_.back()[ 4 ] = RVector3( b, b, a );
    a = 0.399403576166799; b = 0.100596423833201;
    tetAbscissa_.back()[ 5 ] = RVector3( a, b, b );
    tetAbscissa_.back()[ 6 ] = RVector3( b, b, a );
    tetAbscissa_.back()[ 7 ] = RVector3( b, a, a );
    tetAbscissa_.back()[ 8 ] = RVector3( a, a, b );
    tetAbscissa_.back()[ 9 ] = RVector3( b, a, b );
    tetAbscissa_.back()[10 ] = RVector3( a, b, a );

    tetWeights_.push_back( RVector( 11, -0.01315555555555555555555556 ) );
    for ( uint i = 1; i < 5; i ++ ) tetWeights_.back( )[ i ] = 0.007622222222222222222;
    for ( uint i = 5; i < 11; i ++ ) tetWeights_.back( )[ i ] = 0.02488888888888888888888;
    tetWeights_.back() *= 6.0;

    //** 5.Order, n=15, Error: O(h?)
    tetAbscissa_.push_back( std::vector< RVector3 >( 15 ) );
    tetAbscissa_.back()[ 0 ] = RVector3( 0.25, 0.25, 0.25 );
    a = 0.0; b = 1.0 / 3.0;
    tetAbscissa_.back()[ 1 ] = RVector3( b, b, b );
    tetAbscissa_.back()[ 2 ] = RVector3( a, b, b );
    tetAbscissa_.back()[ 3 ] = RVector3( b, a, b );
    tetAbscissa_.back()[ 4 ] = RVector3( b, b, a );
    a = 0.7272727272727272727; b = 0.090909090909090909;
    tetAbscissa_.back()[ 5 ] = RVector3( b, b, b );
    tetAbscissa_.back()[ 6 ] = RVector3( a, b, b );
    tetAbscissa_.back()[ 7 ] = RVector3( b, a, b );
    tetAbscissa_.back()[ 8 ] = RVector3( b, b, a );
    a = 0.066550153573664; b = 0.433449846426336;
    tetAbscissa_.back()[ 9 ] = RVector3( a, b, b );
    tetAbscissa_.back()[10 ] = RVector3( b, b, a );
    tetAbscissa_.back()[11 ] = RVector3( b, a, a );
    tetAbscissa_.back()[12 ] = RVector3( a, a, b );
    tetAbscissa_.back()[13 ] = RVector3( b, a, b );
    tetAbscissa_.back()[14 ] = RVector3( a, b, a );

    tetWeights_.push_back( RVector( 15, 0.030283678097089 ) );
    for ( uint i = 1; i < 5; i ++ ) tetWeights_.back( )[ i ] = 0.006026785714286;
    for ( uint i = 5; i < 9; i ++ ) tetWeights_.back( )[ i ] = 0.011645249086029;
    for ( uint i = 9; i < 15; i ++ ) tetWeights_.back( )[ i ] = 0.010949141561386;
    tetWeights_.back() *= 6.0;
}

void IntegrationRules::initQua_(){

    quaAbscissa_.push_back( std::vector < RVector3 > ( 0 ) );
    quaWeights_.push_back( RVector( 0, 0.0 ) );

    for ( uint order = 1; order < gauAbscissa_.size(); order ++ ){
        uint nK = gauAbscissa_[ order ].size();
        quaAbscissa_.push_back( std::vector < RVector3 > ( nK * nK ) );
        quaWeights_.push_back( RVector( gauAbscissa_[ order ].size() * gauAbscissa_[ order ].size() ) );

        for ( uint i = 0; i < nK; i ++ ){
            for ( uint j = 0; j < nK; j ++ ){
                uint k = i * nK + j;
                quaAbscissa_[ order ][ k ] = RVector3( edgAbscissa_[ order ][ i ][ 0 ], edgAbscissa_[ order ][ j ][ 0 ] );
                quaWeights_[ order ][ k ] = edgWeights_[ order ][ i ] * edgWeights_[ order ][ j ];
            }
        }
    }
}

void IntegrationRules::initHex_(){
}

std::ostream & operator << ( std::ostream & str, const ElementMatrix< double > & e ){
    for ( uint i = 0; i < e.idx().size(); i ++ ) str << e.idx(i) << " " ;

    str << std::endl;
    for ( uint i = 0; i < e.size(); i ++ ){
        str << e.idx( i ) << "\t: ";
        for ( uint j = 0; j < e.size(); j ++ ){
            str << e.getVal( i , j ) << " ";
        }
        str << std::endl;
    }
    return str;
}

template < > ElementMatrix < double > & ElementMatrix < double >::u( const MeshEntity & ent,
                                const RVector & w, const std::vector < RVector3 > & x, bool verbose ){

    uint nVerts = ent.nodeCount();
    std::map< uint, RVector >::const_iterator it = uCache_.find( ent.rtti() );

    if ( it == uCache_.end() ) {
        uint nRules = w.size();

        RVector u( nVerts );
        RMatrix N( nVerts, nRules );

        RVector tmp;
        for ( uint i = 0; i < nRules; i ++ ){
            ent.shapeFunctionsL( x[ i ], tmp );
            N.setCol( i, tmp );
        }
        for ( uint i = 0; i < nVerts; i ++ ){
            u[ i ] = sum( w * N[ i ] );
        }
        uCache_[ ent.rtti() ] = u;
        it = uCache_.find( ent.rtti() );
    }

    double A = ent.shape().domainSize();
    for ( uint i = 0; i < nVerts; i ++ ){
        mat_[ 0 ][ i ] = A * it->second[ i ];
    }

    if ( verbose ) std::cout << "int u " << *this << std::endl;
    return *this;
}

template < > ElementMatrix < double > & ElementMatrix < double >::u2( const MeshEntity & ent,
                                const RVector & w, const std::vector < RVector3 > & x, bool verbose ){

    uint nVerts = ent.nodeCount();
    std::map< uint, RMatrix>::const_iterator it = u2Cache_.find( ent.rtti() );

    if ( it == u2Cache_.end() ) {
        uint nRules = w.size();

        RMatrix u2( nVerts, nVerts );
        RMatrix N( nVerts, nRules );

        RVector tmp;
        for ( uint i = 0; i < nRules; i ++ ){
            ent.shapeFunctionsL( x[ i ], tmp );
            N.setCol( i, tmp );
        }
        for ( uint i = 0; i < nVerts; i ++ ){
            for ( uint j = i; j < nVerts; j ++ ){
                u2[ i ][ j ] = sum( w * N[ j ] * N[ i ] );
                u2[ j ][ i ] = u2[ i ][ j ];
            }
        }
        u2Cache_[ ent.rtti() ] = u2;
        it = u2Cache_.find( ent.rtti() );
    }

    double A = ent.shape().domainSize();
//** very slow yet, pimp this with expressions ( matrix::operator = ( matExpression & ) )
//    mat_ = it->second * A;
    for ( uint i = 0; i < nVerts; i ++ ){
        for ( uint j = 0; j < nVerts; j ++ ){
            mat_[ i ][ j ] = A * it->second[ i ][ j ];
        }
    }

    if ( verbose ) std::cout << "int u2 " << *this << std::endl;

    return *this;
}

template < > ElementMatrix < double > & ElementMatrix < double >::ux2( const MeshEntity & ent,
                                const RVector & w, const std::vector < RVector3 > & x, bool verbose ){
                                 uint nVerts = ent.nodeCount();
    uint nRules = w.size();

    if ( dNdr_.rows() != nVerts ){
        dNdr_.resize( nVerts, nRules );
        for ( uint i = 0; i < nRules; i ++ ){
            dNdr_.setCol( i, ent.deriveNdL( x[ i ], 0 ) );
        }
    }

    double drdx = ent.shape().deriveCoordinates( 1, 0 );

    double A = ent.shape().domainSize();
    for ( uint i = 0; i < nVerts; i ++ ){
        for ( uint j = i; j < nVerts; j ++ ){
            mat_[ i ][ j ] = A * sum( w * ( drdx * dNdr_[ i ] * drdx * dNdr_[ j ] ) );
            mat_[ j ][ i ] = mat_[ i ][ j ];
        }
    }
    if ( verbose ) std::cout << "int ux2uy2 " << *this << std::endl;
    return *this;
}

template < > ElementMatrix < double > & ElementMatrix < double >::ux2uy2( const MeshEntity & ent,
                                const RVector & w, const std::vector < RVector3 > & x, bool verbose ){

    uint nVerts = ent.nodeCount();
    uint nRules = w.size();

    if ( dNdr_.rows() != nVerts ){
        dNdr_.resize( nVerts, nRules );
        dNds_.resize( nVerts, nRules );

        for ( uint i = 0; i < nRules; i ++ ){
            dNdr_.setCol( i, ent.deriveNdL( x[ i ], 0 ) );
            dNds_.setCol( i, ent.deriveNdL( x[ i ], 1 ) );
        }
    }

    double drdx = ent.shape().deriveCoordinates( 0, 0 );
    double dsdx = ent.shape().deriveCoordinates( 1, 0 );
    double drdy = ent.shape().deriveCoordinates( 0, 1 );
    double dsdy = ent.shape().deriveCoordinates( 1, 1 );

    double A = ent.shape().domainSize();
    for ( uint i = 0; i < nVerts; i ++ ){
        for ( uint j = i; j < nVerts; j ++ ){
            mat_[ i ][ j ] = A * sum( w * ( ( drdx * dNdr_[ i ] + dsdx * dNds_[ i ] ) * ( drdx * dNdr_[ j ] + dsdx * dNds_[ j ] ) +
                                            ( drdy * dNdr_[ i ] + dsdy * dNds_[ i ] ) * ( drdy * dNdr_[ j ] + dsdy * dNds_[ j ] ) ) );
            mat_[ j ][ i ] = mat_[ i ][ j ];
        }
    }
    if ( verbose ) std::cout << "int ux2uy2 " << *this << std::endl;
    return *this;
}

template < > ElementMatrix < double > & ElementMatrix < double >::ux2uy2uz2( const MeshEntity & ent,
                                const RVector & w, const std::vector < RVector3 > & x, bool verbose ){

    uint nVerts = ent.nodeCount();
    uint nRules = w.size();

    if ( dNdr_.rows() != nVerts ){
        dNdr_.resize( nVerts, nRules );
        dNds_.resize( nVerts, nRules );
        dNdt_.resize( nVerts, nRules );

        for ( uint i = 0; i < nRules; i ++ ){
            dNdr_.setCol( i, ent.deriveNdL( x[ i ], 0 ) );
            dNds_.setCol( i, ent.deriveNdL( x[ i ], 1 ) );
            dNdt_.setCol( i, ent.deriveNdL( x[ i ], 2 ) );
        }
    }

    double drdx = ent.shape().deriveCoordinates( 0, 0 );
    double dsdx = ent.shape().deriveCoordinates( 1, 0 );
    double dtdx = ent.shape().deriveCoordinates( 2, 0 );
    double drdy = ent.shape().deriveCoordinates( 0, 1 );
    double dsdy = ent.shape().deriveCoordinates( 1, 1 );
    double dtdy = ent.shape().deriveCoordinates( 2, 1 );
    double drdz = ent.shape().deriveCoordinates( 0, 2 );
    double dsdz = ent.shape().deriveCoordinates( 1, 2 );
    double dtdz = ent.shape().deriveCoordinates( 2, 2 );

    double A = ent.shape().domainSize();
    for ( uint i = 0; i < nVerts; i ++ ){
        for ( uint j = i; j < nVerts; j ++ ){
            mat_[ i ][ j ] = A * sum( w * ( ( drdx * dNdr_[ i ] + dsdx * dNds_[ i ] + dtdx * dNdt_[ i ] ) *
                                            ( drdx * dNdr_[ j ] + dsdx * dNds_[ j ] + dtdx * dNdt_[ j ] ) +
                                            ( drdy * dNdr_[ i ] + dsdy * dNds_[ i ] + dtdy * dNdt_[ i ] ) *
                                            ( drdy * dNdr_[ j ] + dsdy * dNds_[ j ] + dtdy * dNdt_[ j ] ) +
                                            ( drdz * dNdr_[ i ] + dsdz * dNds_[ i ] + dtdz * dNdt_[ i ] ) *
                                            ( drdz * dNdr_[ j ] + dsdz * dNds_[ j ] + dtdz * dNdt_[ j ] )
                                        )
                                    );

            mat_[ j ][ i ] = mat_[ i ][ j ];
        }
    }
    if ( verbose ) std::cout << "int ux2uy2uz2 " << *this << std::endl;
    return *this;
}

template < > ElementMatrix < double > & ElementMatrix < double >::u( const MeshEntity & ent ){
    uint nVerts = ent.nodeCount();
    if ( size() != nVerts ) resize( nVerts );

    for ( uint i = 0; i < nVerts; i ++ ) idx_[ i ] = ent.node( i ).id();

    *this *= 0.0;

    switch( ent.rtti() ){
        case MESH_BOUNDARY_NODE_RTTI:
            mat_[ 0 ][ 0 ] = 1.0;
        break;
        case MESH_EDGE_CELL_RTTI:
        case MESH_EDGE_RTTI:
        case MESH_EDGE3_CELL_RTTI:
        case MESH_EDGE3_RTTI:
            return u( ent, intRules_.edgWeights( 2 ), intRules_.edgAbscissa( 2 ), false ); //ch
        case MESH_TRIANGLE_RTTI:
        case MESH_TRIANGLEFACE_RTTI:
        case MESH_TRIANGLE6_RTTI:
        case MESH_TRIANGLEFACE6_RTTI:
            return u( ent, intRules_.triWeights( 2 ), intRules_.triAbscissa( 2 ), false ); //ch
        case MESH_QUADRANGLE_RTTI:
        case MESH_QUADRANGLE8_RTTI:
            return u( ent, intRules_.quaWeights( 2 ), intRules_.quaAbscissa( 2 ), false ); //ch
        case MESH_TETRAHEDRON_RTTI:
        case MESH_TETRAHEDRON10_RTTI:
            return u( ent, intRules_.tetWeights( 2 ), intRules_.tetAbscissa( 2 ), false );

        default: std::cerr << WHERE_AM_I << " celltype not spezified " << ent.rtti() << std::endl;
    }
    return *this;
}

template < > ElementMatrix < double > & ElementMatrix < double >::u2( const MeshEntity & ent ){
    uint nVerts = ent.nodeCount();
    if ( size() != nVerts ) resize( nVerts );
    for ( uint i = 0; i < nVerts; i ++ ) idx_[ i ] = ent.node( i ).id();

    switch( ent.rtti() ){
    case MESH_BOUNDARY_NODE_RTTI:
        mat_[ 0 ][ 0 ] = 1.0;
    break;
    case MESH_EDGE_CELL_RTTI:
    case MESH_EDGE_RTTI:
    //{
//         double J = ent.shape().jacobianDeterminant();
//         if ( J < 0 ) {
//             std::cerr << WHERE_AM_I << "JacobianDeterminant < 0 (" << J << ")" << std::endl;
//             std::cerr << ent.shape() << std::endl;
//             J = std::fabs( J );
//         }
//         mat_[ 0 ][ 0 ] = J / 3.0;
//         mat_[ 1 ][ 0 ] = J / 6.0;
//         mat_[ 0 ][ 1 ] = J / 6.0;
//         mat_[ 1 ][ 1 ] = J / 3.0;
//         std::cout << "2 " << *this << std::endl;
/*}
    break;*/
        return u2( ent, intRules_.edgWeights( 2 ), intRules_.edgAbscissa( 2 ), false ); //ch
    case MESH_EDGE3_CELL_RTTI:
    case MESH_EDGE3_RTTI:
    //{
//         double J = ent.shape().jacobianDeterminant();
//         if ( J < 0 ) {
//             std::cerr << WHERE_AM_I << "JacobianDeterminant < 0 (" << J << ")" << std::endl;
//             std::cerr << ent.shape() << std::endl;
//             J = std::fabs( J );
//         }
//         mat_[ 0 ][ 0 ] =   J / 30.0 * 4.0;
//         mat_[ 1 ][ 0 ] = - J / 30.0;
//         mat_[ 2 ][ 0 ] =   J / 15.0;
//         mat_[ 0 ][ 1 ] = - J / 30.0;
//         mat_[ 1 ][ 1 ] =   J / 30.0 * 4.0;
//         mat_[ 2 ][ 1 ] =   J / 15.0;
//         mat_[ 0 ][ 2 ] =   J / 15.0;
//         mat_[ 1 ][ 2 ] =   J / 15.0;
//         mat_[ 2 ][ 2 ] =   J / 30.0 * 16.0;
//          std::cout << "2 " << *this << std::endl;
    //} break;
        return u2( ent, intRules_.edgWeights( 3 ), intRules_.edgAbscissa( 3 ), false ); //ch
    case MESH_TRIANGLE_RTTI:
    case MESH_TRIANGLEFACE_RTTI:
    //{
//         double J = ent.shape().jacobianDeterminant();
//         if ( J < 0 ) {
//             std::cerr << WHERE_AM_I << "JacobianDeterminant < 0 (" << J << ")" << std::endl;
//             std::cerr << ent.shape() << std::endl;
//             J = std::fabs( J );
//         }
//         double Jl = J / 24.0;
//         mat_[ 0 ][ 0 ] = Jl * 2.0;
//         mat_[ 1 ][ 0 ] = Jl;
//         mat_[ 2 ][ 0 ] = Jl;
//         mat_[ 0 ][ 1 ] = Jl;
//         mat_[ 1 ][ 1 ] = Jl * 2.0;
//         mat_[ 2 ][ 1 ] = Jl;
//         mat_[ 0 ][ 2 ] = Jl;
//         mat_[ 1 ][ 2 ] = Jl;
//         mat_[ 2 ][ 2 ] = Jl * 2.0;
//         std::cout << "2 " << *this << std::endl;
//}break;
        return u2( ent, intRules_.triWeights( 2 ), intRules_.triAbscissa( 2 ), false ); //ch
    case MESH_QUADRANGLE_RTTI:
    case MESH_QUADRANGLEFACE_RTTI:
        return u2( ent, intRules_.quaWeights( 2 ), intRules_.quaAbscissa( 2 ), false );
    case MESH_QUADRANGLE8_RTTI:
    case MESH_QUADRANGLEFACE8_RTTI:
        return u2( ent, intRules_.quaWeights( 3 ), intRules_.quaAbscissa( 3 ), false );
    case MESH_TRIANGLE6_RTTI:
    case MESH_TRIANGLEFACE6_RTTI:
    //{
//         double J = ent.shape().jacobianDeterminant();
//         if ( J < 0 ) {
//             std::cerr << WHERE_AM_I << "JacobianDeterminant < 0 (" << J << ")" << std::endl;
//             std::cerr << ent.shape() << std::endl;
//             J = std::fabs( J );
//         }
//         for ( uint i = 0; i < nVerts; i++ ){
//             for ( uint j = 0; j < nVerts; j++ ){
//                 mat_[ i ][ j ] = J * (*Tri6_u2)[i][j];//compound[ i ][ j ];
//             }
//         }
//         std::cout << "2 " << *this << std::endl;
    //} break;
        return u2( ent, intRules_.triWeights( 4 ), intRules_.triAbscissa( 4 ), false ); //ch
    case MESH_TETRAHEDRON_RTTI:
        return u2( ent, intRules_.tetWeights( 2 ), intRules_.tetAbscissa( 2 ), false );
    case MESH_TETRAHEDRON10_RTTI:
        return u2( ent, intRules_.tetWeights( 4 ), intRules_.tetAbscissa( 4 ), false );
    default:
      std::cerr << ent.rtti() << std::endl;
      THROW_TO_IMPL
    }

    return *this;
}

template < > ElementMatrix < double > & ElementMatrix < double >::ux2uy2uz2( const Cell & cell ){

    uint dim = cell.nodeCount();
    if ( size() != dim ) resize( dim );
    for ( uint i = 0; i < dim; i ++ ) idx_[ i ] = cell.node( i ).id();

//     double J = cell.jacobianDeterminant();
//     if ( J <= 0 ) std::cerr << WHERE_AM_I << " JacobianDeterminant < 0 (" << J << ") " << cell << std::endl;
//      std::cout << J << std::endl;

    switch ( cell.rtti() ) {
    case MESH_EDGE_CELL_RTTI:
    case MESH_EDGE3_CELL_RTTI:
        return ux2( cell, intRules_.edgWeights( 2 ), intRules_.edgAbscissa( 2 ), false );
    case MESH_TRIANGLE_RTTI: {
    ////////////////////////////////////////////////////////////////////
/*        double dN1dx = cell.shape().deriveCoordinates( 0, 0 );
        double dN2dx = cell.shape().deriveCoordinates( 1, 0 );
        double dN3dx = cell.shape().deriveCoordinates( 2, 0 );
        double dN1dy = cell.shape().deriveCoordinates( 0, 1 );
        double dN2dy = cell.shape().deriveCoordinates( 1, 1 );
        double dN3dy = cell.shape().deriveCoordinates( 2, 1 );
        mat_[ 0 ][ 0 ] = J / 2.0 * ( dN1dx * dN1dx + dN1dy * dN1dy );
        mat_[ 1 ][ 0 ] = J / 2.0 * ( dN2dx * dN1dx + dN2dy * dN1dy );
        mat_[ 2 ][ 0 ] = J / 2.0 * ( dN3dx * dN1dx + dN3dy * dN1dy );
        mat_[ 0 ][ 1 ] = J / 2.0 * ( dN1dx * dN2dx + dN1dy * dN2dy );
        mat_[ 1 ][ 1 ] = J / 2.0 * ( dN2dx * dN2dx + dN2dy * dN2dy );
        mat_[ 2 ][ 1 ] = J / 2.0 * ( dN3dx * dN2dx + dN3dy * dN2dy );
        mat_[ 0 ][ 2 ] = J / 2.0 * ( dN1dx * dN3dx + dN1dy * dN3dy );
        mat_[ 1 ][ 2 ] = J / 2.0 * ( dN2dx * dN3dx + dN2dy * dN3dy );
        mat_[ 2 ][ 2 ] = J / 2.0 * ( dN3dx * dN3dx + dN3dy * dN3dy );*/
////////////////////////////////////////////////////////////////////
/*        double x1 = cell.node( 0 ).x();
        double x2 = cell.node( 1 ).x();
        double x3 = cell.node( 2 ).x();
        double y1 = cell.node( 0 ).y();
        double y2 = cell.node( 1 ).y();
        double y3 = cell.node( 2 ).y();

        double a = ( (x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1) ) / J;
        double b = - ( (x3 - x1) * (x2 - x1) + (y3 - y1) * (y2 - y1) ) / J;
        double c = ( (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) ) / J;

        mat_[ 0 ][ 0 ] = a *  0.5 + b        + c *  0.5 ;
        mat_[ 1 ][ 0 ] = a * -0.5 + b * -0.5            ;
        mat_[ 2 ][ 0 ] =            b * -0.5 + c * -0.5 ;
        mat_[ 1 ][ 1 ] = a *  0.5                       ;
        mat_[ 2 ][ 1 ] =            b *  0.5            ;
        mat_[ 2 ][ 2 ] =                       c *  0.5 ;

        mat_[ 0 ][ 1 ] = mat_[ 1 ][ 0 ];
        mat_[ 0 ][ 2 ] = mat_[ 2 ][ 0 ];
        mat_[ 1 ][ 2 ] = mat_[ 2 ][ 1 ];

        std::cout << "2" << *this << std::endl;*/
        return ux2uy2( cell, intRules_.triWeights( 1 ), intRules_.triAbscissa( 1 ), false );
    } break;
    case MESH_TRIANGLE6_RTTI: {
///////////////////////////////////////////////////////////////////////////////////
//         double x1 = cell.node( 0 ).x();
//         double x2 = cell.node( 1 ).x();
//         double x3 = cell.node( 2 ).x();
//         double y1 = cell.node( 0 ).y();
//         double y2 = cell.node( 1 ).y();
//         double y3 = cell.node( 2 ).y();
//
//         double b1 = y2-y3;
//         double b2 = y3-y1;
//         double b3 = y1-y2;
//         double c1 = x3-x2;
//         double c2 = x1-x3;
//         double c3 = x2-x1;
//         double a = (c3*b2-c2*b3)/2.0;
//         double b1sqr=b1*b1;
//         double b2sqr=b2*b2;
//         double b3sqr=b3*b3;
//         double c1sqr=c1*c1;
//         double c2sqr=c2*c2;
//         double c3sqr=c3*c3;
//         double b1b2=b1*b2;
//         double b1b3=b1*b3;
//         double b2b3=b2*b3;
//         double c1c2=c1*c2;
//         double c1c3=c1*c3;
//         double c2c3=c2*c3;
//
//         mat_[0][0]=(b1sqr+c1sqr) / (4.0*a);
//         mat_[0][1]=(-b1b2-c1c2)/(12.0*a);
//         mat_[0][2]=(-b1b3-c1c3)/(12.0*a);
//         mat_[0][3]=(b1b2+c1c2)/(3.0*a);
//         mat_[0][4]=0.0;
//         mat_[0][5]=(b1b3+c1c3)/(3.0*a);
//         mat_[1][1]=(b2sqr+c2sqr)/(4.0*a);
//         mat_[1][2]=(-b2b3-c2c3)/(12.0*a);
//         mat_[1][3]=(b1b2+c1c2)/(3.0*a);
//         mat_[1][4]=(b2b3+c2c3)/(3.0*a);
//         mat_[1][5]=0.0;
//         mat_[2][2]=(b3sqr+c3sqr)/(4.0*a);
//         mat_[2][3]=0.0;
//         mat_[2][4]=(b2b3+c2c3)/(3.0*a);
//         mat_[2][5]=(b1b3+c1c3)/(3.0*a);
//         mat_[3][3]=2.0*((b1sqr+b1b2+b2sqr)+(c1sqr+c1c2+c2sqr))/(3.0*a);
//         mat_[3][4]=((b1b2+b2sqr+2.0*b1b3+b2b3)+(c1c2+c2sqr+2.0*c1c3+c2c3))/(3.0*a);
//         mat_[3][5]=((b1sqr+b1b3+b1b2+2.0*b2b3)+(c1sqr+c1c3+c1c2+2.0*c2c3))/(3.0*a);
//         mat_[4][4]=2.0*((b2sqr+b2b3+b3sqr)+(c2sqr+c2c3+c3sqr))/(3.0*a);
//         mat_[4][5]=((2.0*b1b2+b2b3+b1b3+b3sqr)+(2.0*c1c2+c2c3+c1c3+c3sqr))/(3.0*a);
//         mat_[5][5]=2.0*((b1sqr+b1b3+b3sqr)+(c1sqr+c1c3+c3sqr))/(3.0*a);
//
//         for ( int i = 1, imax = 6; i < imax; i ++ ){
//              for ( int j = 0, jmax = i; j < jmax; j ++ ){
//                 mat_[ i ][ j ]=mat_[ j ][ i ];
//              }
//         }

//working version
//         double x1 = cell.node( 0 ).x();
//         double x2 = cell.node( 1 ).x();
//         double x3 = cell.node( 2 ).x();
//         double y1 = cell.node( 0 ).y();
//         double y2 = cell.node( 1 ).y();
//         double y3 = cell.node( 2 ).y();
//
//         double a = ( (x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1) ) / J / 6.0;
//         double b = - ( (x3 - x1) * (x2 - x1) + (y3 - y1) * (y2 - y1) ) / J / 6.0;
//         double c = ( (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) ) / J / 6.0;
//         for ( int i = 0, imax = 6; i < imax; i ++ ){
//             for ( int j = 0, jmax = 6; j < jmax; j ++ ){
//                 mat_[ i ][ j ] = Triangle6_S1[ i ][ j ] * a +
//                                  Triangle6_S2[ i ][ j ] * b +
//                                  Triangle6_S3[ i ][ j ] * c;
//             }
//         }
//         std::cout << "2" << *this << std::endl;
        return ux2uy2( cell, intRules_.triWeights( 2 ), intRules_.triAbscissa( 2 ), false ); //ch
    } break;
    case MESH_QUADRANGLE_RTTI:
        return ux2uy2( cell, intRules_.quaWeights( 2 ), intRules_.quaAbscissa( 2 ), false );
    case MESH_QUADRANGLE8_RTTI:
        return ux2uy2( cell, intRules_.quaWeights( 3 ), intRules_.quaAbscissa( 3 ), false );
    case MESH_TETRAHEDRON_RTTI:
    //{
//         double x_xi = cell.shape().partDerivationRealToUnity( 0, 1 );
//         double x_eta = cell.shape().partDerivationRealToUnity( 0, 2 );
//         double x_zeta = cell.shape().partDerivationRealToUnity( 0, 3 );
//
//         double y_xi = cell.shape().partDerivationRealToUnity( 1, 1 );
//         double y_eta = cell.shape().partDerivationRealToUnity( 1, 2 );
//         double y_zeta = cell.shape().partDerivationRealToUnity( 1, 3 );
//
//         double z_xi = cell.shape().partDerivationRealToUnity( 2, 1 );
//         double z_eta = cell.shape().partDerivationRealToUnity( 2, 2 );
//         double z_zeta = cell.shape().partDerivationRealToUnity( 2, 3 );
//
//         double xi_x =  1.0 / J * det( y_eta, z_eta, y_zeta, z_zeta );
//         double xi_y = -1.0 / J * det( x_eta, z_eta, x_zeta, z_zeta );
//         double xi_z =  1.0 / J * det( x_eta, y_eta, x_zeta, y_zeta );
//
//         double eta_x = -1.0 / J * det( y_xi, z_xi, y_zeta, z_zeta );
//         double eta_y =  1.0 / J * det( x_xi, z_xi, x_zeta, z_zeta );
//         double eta_z = -1.0 / J * det( x_xi, y_xi, x_zeta, y_zeta );
//
//         double zeta_x =  1.0 / J * det( y_xi, z_xi, y_eta, z_eta );
//         double zeta_y = -1.0 / J * det( x_xi, z_xi, x_eta, z_eta );
//         double zeta_z =  1.0 / J * det( x_xi, y_xi, x_eta, y_eta );
//
//         double a = J / 6.0 * ( xi_x * xi_x + xi_y * xi_y + xi_z * xi_z );
//         double b = J / 6.0 * ( eta_x * eta_x + eta_y * eta_y + eta_z * eta_z );
//         double c = J / 6.0 * ( zeta_x * zeta_x + zeta_y * zeta_y + zeta_z * zeta_z );
//
//         double d = J / 6.0 * ( xi_x * eta_x + xi_y * eta_y + xi_z * eta_z );
//         double e = J / 6.0 * ( xi_x * zeta_x + xi_y * zeta_y + xi_z * zeta_z );
//         double f = J / 6.0 * ( eta_x * zeta_x + eta_y * zeta_y + eta_z * zeta_z );
//
//     double u_xi2[ 4 ][ 4 ] = {
//       {   a,  -a, 0.0, 0.0 },
//       {  -a,   a, 0.0, 0.0 },
//       { 0.0, 0.0, 0.0, 0.0 },
//       { 0.0, 0.0, 0.0, 0.0 }
//     };
//     double u_eta2[ 4 ][ 4 ] = {
//       {   b, 0.0,  -b, 0.0 },
//       { 0.0, 0.0, 0.0, 0.0 },
//       {  -b, 0.0,   b, 0.0 },
//       { 0.0, 0.0, 0.0, 0.0 }
//     };
//     double u_zeta2[ 4 ][ 4 ] = {
//       {   c, 0.0, 0.0,  -c },
//       { 0.0, 0.0, 0.0, 0.0 },
//       { 0.0, 0.0, 0.0, 0.0 },
//       {  -c, 0.0, 0.0,   c }
//     };
//     double u_xi__u_eta[ 4 ][ 4 ] = {
//       { 2.0 * d,  -d,  -d, 0.0 },
//       {      -d, 0.0,   d, 0.0 },
//       {      -d,   d, 0.0, 0.0 },
//       {     0.0, 0.0, 0.0, 0.0 }
//     };
//     double u_xi__u_zeta[ 4 ][ 4 ] = {
//       { 2.0 * e,  -e, 0.0,  -e },
//       {      -e, 0.0, 0.0,   e },
//       {     0.0, 0.0, 0.0, 0.0 },
//       {      -e,   e, 0.0, 0.0 }
//     };
//     double u_eta__u_zeta[ 4 ][ 4 ] = {
//       { 2.0 * f, 0.0,  -f,  -f },
//       {     0.0, 0.0, 0.0, 0.0 },
//       {      -f, 0.0, 0.0,   f },
//       {      -f, 0.0,   f, 0.0 }
//     };
//
//     for ( uint i = 0; i < dim; i++ ){
//       for ( uint j = 0; j < dim; j++ ){
//  mat_[ i ][ j ] = u_xi2[ i ][ j ] + u_eta2[ i ][ j ] +  u_zeta2[ i ][ j ] +
//    u_xi__u_eta[ i ][ j ] + u_xi__u_zeta[ i ][ j ] + u_eta__u_zeta[ i ][ j ];
//       }
//     }
//     std::cout << "0 " << *this << std::endl;
//} break;
        return ux2uy2uz2( cell, intRules_.tetWeights( 1 ), intRules_.tetAbscissa( 1 ), false ); //ch
    case MESH_TETRAHEDRON10_RTTI:
//     {
// ////////////////////////////////////////////////////////////////////
//     double x_xi = cell.shape().partDerivationRealToUnity( 0, 1 );
//     double x_eta = cell.shape().partDerivationRealToUnity( 0, 2 );
//     double x_zeta = cell.shape().partDerivationRealToUnity( 0, 3 );
//
//     double y_xi = cell.shape().partDerivationRealToUnity( 1, 1 );
//     double y_eta = cell.shape().partDerivationRealToUnity( 1, 2 );
//     double y_zeta = cell.shape().partDerivationRealToUnity( 1, 3 );
//
//     double z_xi = cell.shape().partDerivationRealToUnity( 2, 1 );
//     double z_eta = cell.shape().partDerivationRealToUnity( 2, 2 );
//     double z_zeta = cell.shape().partDerivationRealToUnity( 2, 3 );
//
//     double xi_x =  1.0 / J * det( y_eta, z_eta, y_zeta, z_zeta );
//     double xi_y = -1.0 / J * det( x_eta, z_eta, x_zeta, z_zeta );
//     double xi_z =  1.0 / J * det( x_eta, y_eta, x_zeta, y_zeta );
//
//     double eta_x = -1.0 / J * det( y_xi, z_xi, y_zeta, z_zeta );
//     double eta_y =  1.0 / J * det( x_xi, z_xi, x_zeta, z_zeta );
//     double eta_z = -1.0 / J * det( x_xi, y_xi, x_zeta, y_zeta );
//
//     double zeta_x =  1.0 / J * det( y_xi, z_xi, y_eta, z_eta );
//     double zeta_y = -1.0 / J * det( x_xi, z_xi, x_eta, z_eta );
//     double zeta_z =  1.0 / J * det( x_xi, y_xi, x_eta, y_eta );
//
//     double a = J / 6.0 * ( xi_x * xi_x + xi_y * xi_y + xi_z * xi_z );
//     double b = J / 6.0 * ( eta_x * eta_x + eta_y * eta_y + eta_z * eta_z );
//     double c = J / 6.0 * ( zeta_x * zeta_x + zeta_y * zeta_y + zeta_z * zeta_z );
//
//     double d = J / 6.0 * ( xi_x * eta_x + xi_y * eta_y + xi_z * eta_z );
//     double e = J / 6.0 * ( xi_x * zeta_x + xi_y * zeta_y + xi_z * zeta_z );
//     double f = J / 6.0 * ( eta_x * zeta_x + eta_y * zeta_y + eta_z * zeta_z );
//
//     RSTLMatrix compound( 10, 10 );
//     compound = a * (*Tet10_u_xi2) + b * (*Tet10_u_eta2) + c * (*Tet10_u_zeta2)
//         + (2.0*d) * (*Tet10_u_xi_u_eta)
//         + (2.0*e) * (*Tet10_u_xi_u_zeta)
//         + (2.0*f) * (*Tet10_u_eta_u_zeta);
//     for ( uint i = 0; i < dim; i++ ){
//       for ( uint j = 0; j < dim; j++ ){
//         //** * 6.0 weil a b c d e f / 6
//         mat_[ i ][ j ] = compound[ i ][ j ] * 6.0;
//       }
//     }
// //    std::cout << " 2 " << *this << std::endl;
//
//     break;
//   }
  return ux2uy2uz2( cell, intRules_.tetWeights( 2 ), intRules_.tetAbscissa( 2 ), false );
  default:
    std::cerr << cell.rtti() << std::endl;
    THROW_TO_IMPL
    break;
  }
  return *this;
}

} // namespace GIMLI
