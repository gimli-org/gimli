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

#ifndef _GIMLI_ELEMENTMATRIX__H
#define _GIMLI_ELEMENTMATRIX__H

#include "gimli.h"
#include "mesh.h"
#include "shape.h"
#include "node.h"
#include "vectortemplates.h"
#include "numericbase.h"
#include "stopwatch.h"
#include "matrix.h"

namespace GIMLI{

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

//   static double Triangle6_S4[ 6 ][ 6 ] = {
//     {  6.0, -1.0, -1.0,  0.0, -4.0,  0.0 },
//     { -1.0,  6.0, -1.0,  0.0,  0.0, -4.0 },
//     { -1.0, -1.0,  6.0, -4.0,  0.0,  0.0 },
//     {  0.0,  0.0, -4.0, 32.0, 16.0, 16.0 },
//     { -4.0,  0.0,  0.0, 16.0, 32.0, 16.0 },
//     {  0.0, -4.0,  0.0, 16.0, 16.0, 32.0 }
//   };
//  static double Edge3_Me[ 3 ][ 3 ] = { // Me = edge.length / 30.0 * Edge3_Me
//     {  4.0,  2.0, -1.0 },
//     {  2.0, 16.0,  2.0 },
//     { -1.0,  2.0,  4.0 }
//   };
//   static double Edge3_Se[ 3 ][ 3 ] = { // Se = 1 / ( 3.0 * edge.length ) * Edge3_Se
//     {  7.0, -8.0,  1.0 },
//     { -8.0, 16.0, -8.0 },
//     {  1.0, -8.0,  7.0 }
//   };

template < class T > class DLLEXPORT ElementMatrix {
public:
    ElementMatrix( ) { initBaseMatricies(); }

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

  ~ElementMatrix() {
/*    delete Tri6_u_xi2;
    delete Tri6_u_eta2;
    delete Tri6_u_xi_u_eta;
    delete Tri6_u2;

    delete Tet4_u_xi2;
    delete Tet4_u_eta2;
    delete Tet4_u_zeta2;
    delete Tet4_u_xi_u_eta;
    delete Tet4_u_xi_u_zeta;
    delete Tet4_u_eta_u_zeta;
    delete Tet4_u2;

    delete Tet10_u_xi2;
    delete Tet10_u_eta2;
    delete Tet10_u_zeta2;
    delete Tet10_u_xi_u_eta;
    delete Tet10_u_xi_u_zeta;
    delete Tet10_u_eta_u_zeta;
    delete Tet10_u2;*/
  }

    inline const RVector & operator[]( uint row ) const { return mat_[ row ]; }

   void initBaseMatricies(){
//     Tri6_u_xi2 = new RSTLMatrix( 6, 6 );
//     Tri6_u_eta2 = new RSTLMatrix( 6, 6 );
//     Tri6_u_xi_u_eta = new RSTLMatrix( 6, 6 );
//     Tri6_u2 = new RSTLMatrix( 6, 6 );
//
//     Tet4_u_xi2 = new RSTLMatrix( 4, 4 );
//     Tet4_u_eta2 = new RSTLMatrix( 4, 4 );
//     Tet4_u_zeta2 = new RSTLMatrix( 4, 4 );
//     Tet4_u_xi_u_eta = new RSTLMatrix( 4, 4 );
//     Tet4_u_xi_u_zeta = new RSTLMatrix( 4, 4 );
//     Tet4_u_eta_u_zeta = new RSTLMatrix( 4, 4 );
//     Tet4_u2= new RSTLMatrix( 4, 4 );
//
//     Tet10_u_xi2 = new RSTLMatrix( 10, 10 );
//     Tet10_u_eta2 = new RSTLMatrix( 10, 10 );
//     Tet10_u_zeta2 = new RSTLMatrix( 10, 10 );
//     Tet10_u_xi_u_eta = new RSTLMatrix( 10, 10 );
//     Tet10_u_xi_u_zeta = new RSTLMatrix( 10, 10 );
//     Tet10_u_eta_u_zeta = new RSTLMatrix( 10, 10 );
//     Tet10_u2= new RSTLMatrix( 10, 10 );
//
//     *Tri6_u_xi2 = transpose( TriangleQuadA() )
//           * createMatrix( 6, intQuadTriangle_u_xi2() ) * TriangleQuadA();
//     *Tri6_u_eta2 = transpose( TriangleQuadA() )
//         * createMatrix( 6, intQuadTriangle_u_eta2() ) * TriangleQuadA();
//     *Tri6_u_xi_u_eta = transpose( TriangleQuadA() )
//         * createMatrix( 6, intQuadTriangle_u_xi_u_eta() ) * TriangleQuadA();
//     *Tri6_u2 = transpose( TriangleQuadA() )
//         * createMatrix( 6, intQuadTriangle_u2() ) * TriangleQuadA();
//
//     *Tet4_u_xi2 = transpose( TetrahedronLinearA() )
//         * createMatrix( 4, intLinTetrahedron_u_xi2() ) * TetrahedronLinearA();
//     *Tet4_u_eta2 = transpose( TetrahedronLinearA() )
//         * createMatrix( 4, intLinTetrahedron_u_eta2() ) * TetrahedronLinearA();
//     *Tet4_u_zeta2 = transpose( TetrahedronLinearA() )
//         * createMatrix( 4, intLinTetrahedron_u_zeta2() ) * TetrahedronLinearA();
//     *Tet4_u_xi_u_eta = transpose( TetrahedronLinearA() )
//         * createMatrix( 4, intLinTetrahedron_u_xi_u_eta() ) * TetrahedronLinearA();
//     *Tet4_u_xi_u_zeta = transpose( TetrahedronLinearA() )
//         * createMatrix( 4, intLinTetrahedron_u_xi_u_zeta() ) * TetrahedronLinearA();
//     *Tet4_u_eta_u_zeta = transpose( TetrahedronLinearA() )
//         * createMatrix( 4, intLinTetrahedron_u_eta_u_zeta() ) * TetrahedronLinearA();
//
//     *Tet10_u_xi2 = transpose( TetrahedronQuadA() )
//         * createMatrix( 10, intQuadTetrahedron_u_xi2() ) * TetrahedronQuadA();
//     *Tet10_u_eta2 = transpose( TetrahedronQuadA() )
//         * createMatrix( 10, intQuadTetrahedron_u_eta2() ) * TetrahedronQuadA();
//     *Tet10_u_zeta2 = transpose( TetrahedronQuadA() )
//         * createMatrix( 10, intQuadTetrahedron_u_zeta2() ) * TetrahedronQuadA();
//     *Tet10_u_xi_u_eta = transpose( TetrahedronQuadA() )
//         * createMatrix( 10, intQuadTetrahedron_u_xi_u_eta() ) * TetrahedronQuadA();
//     *Tet10_u_xi_u_zeta = transpose( TetrahedronQuadA() )
//         * createMatrix( 10, intQuadTetrahedron_u_xi_u_zeta() ) * TetrahedronQuadA();
//     *Tet10_u_eta_u_zeta = transpose( TetrahedronQuadA() )
//         * createMatrix( 10, intQuadTetrahedron_u_eta_u_zeta() ) * TetrahedronQuadA();
   }

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
    const T & getVal( uint i, uint j ) const { return mat_[ i ][ j ]; }

    const RVector & row( uint i ) const { return mat_[ i ]; }
    const RMatrix & mat() const { return mat_; }
    const std::vector < uint > & idx() const { return idx_; }

    ElementMatrix < T> & u( const MeshEntity & ent );

    ElementMatrix < T> & u2( const MeshEntity & ent );

    ElementMatrix < T > & ux2uy2uz2( const Cell & cell );

    ElementMatrix < T > & u( const MeshEntity & ent, const RVector & w, const std::vector < RVector3 > & x, bool verbose = false );
    ElementMatrix < T > & u2( const MeshEntity & ent, const RVector & w, const std::vector < RVector3 > & x, bool verbose = false );
    ElementMatrix < T > & ux2( const MeshEntity & ent, const RVector & w, const std::vector < RVector3 > & x, bool verbose = false );
    ElementMatrix < T > & ux2uy2( const MeshEntity & ent, const RVector & w, const std::vector < RVector3 > & x, bool verbose = false );
    ElementMatrix < T > & ux2uy2uz2( const MeshEntity & ent, const RVector & w, const std::vector < RVector3 > & x, bool verbose = false );

    friend std::ostream & operator << ( std::ostream & str, const ElementMatrix< double > & pos );

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


//     RSTLMatrix *Tri3_u_xi2;
//   RSTLMatrix *Tri3_u_eta2;
//   RSTLMatrix *Tri3_u_xi_u_eta;
//   RSTLMatrix *Tri3_u2;
//
//   RSTLMatrix *Tri6_u_xi2;
//   RSTLMatrix *Tri6_u_eta2;
//   RSTLMatrix *Tri6_u_xi_u_eta;
//   RSTLMatrix *Tri6_u2;
//
//   RSTLMatrix *Tet4_u_xi2;
//   RSTLMatrix *Tet4_u_eta2;
//   RSTLMatrix *Tet4_u_zeta2;
//   RSTLMatrix *Tet4_u_xi_u_eta;
//   RSTLMatrix *Tet4_u_xi_u_zeta;
//   RSTLMatrix *Tet4_u_eta_u_zeta;
//   RSTLMatrix *Tet4_u2;
//
//   RSTLMatrix *Tet10_u_xi2;
//   RSTLMatrix *Tet10_u_eta2;
//   RSTLMatrix *Tet10_u_zeta2;
//   RSTLMatrix *Tet10_u_xi_u_eta;
//   RSTLMatrix *Tet10_u_xi_u_zeta;
//   RSTLMatrix *Tet10_u_eta_u_zeta;
//   RSTLMatrix *Tet10_u2;
};

} // namespace GIMLI{

#endif // _GIMLI_ELEMENTMATRIX__H
