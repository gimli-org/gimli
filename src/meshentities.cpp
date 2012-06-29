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

#include "meshentities.h"
#include "node.h"
#include "shape.h"

#include <map>
#include <algorithm>

namespace GIMLI{

void Edge_shapeFunctionsL( const RVector3 & L, RVector & funct ) {
    funct.resize( 2 );
    funct[ 0 ] = 1.0 - L[ 0 ];
    funct[ 1 ] = L[ 0 ];
}

void Edge3_shapeFunctionsL( const RVector3 & L, RVector & funct ) {
    funct.resize( 3 );
    double L1 = 1.0 - L[ 0 ];
    double L2 = L[ 0 ];
    funct[ 0 ] = L1 * ( L1 - L2 );
    funct[ 1 ] = L2 * ( L2 - L1 );
    funct[ 2 ] = 4.0 * L1 * L2;
}

void Triangle_shapeFunctionsL( const RVector3 & L, RVector & funct ) {
    funct.resize( 3 );
    funct[ 0 ] = 1.0 - L[ 0 ] - L[ 1 ];
    funct[ 1 ] = L[ 0 ];
    funct[ 2 ] = L[ 1 ];
}

void Triangle6_shapeFunctionsL( const RVector3 & L, RVector & funct ) {
    funct.resize( 6 );
    double L1 = 1.0 - L[ 0 ] - L[ 1 ];
    double L2 = L[ 0 ];
    double L3 = L[ 1 ];
    funct[ 0 ] = L1 * ( 2.0 * L1 - 1.0 );
    funct[ 1 ] = L2 * ( 2.0 * L2 - 1.0 );
    funct[ 2 ] = L3 * ( 2.0 * L3 - 1.0 );
    funct[ 3 ] = 4.0 * L1 * L2;
    funct[ 4 ] = 4.0 * L2 * L3;
    funct[ 5 ] = 4.0 * L3 * L1;
}

void Tetrangle_shapeFunctionsL( const RVector3 & L, RVector & funct ) {
    funct.resize( 4 );
    funct[ 0 ] = 1.0 - L[ 0 ] - L[ 1 ] - L[ 2 ];
    funct[ 1 ] = L[ 0 ];
    funct[ 2 ] = L[ 1 ];
    funct[ 3 ] = L[ 2 ];
}

Boundary * findBoundary_( const std::set < Boundary *> & common ){
    if ( common.size() == 1 ) {
        return *common.begin();
    } else {
        if ( common.size() > 1 ){
            std::cerr << WHERE_AM_I << " pls. check, this should not happen. "
                    " There is more then one boundary defined." <<
                    common.size() << std::endl;
            std::for_each( common.begin(), common.end(), cerrPtrObject() );
        }
    }
    return NULL;
}

Boundary * findBoundary( const Node & n1 ) { return findBoundary_( n1.boundSet() ); }

Boundary * findBoundary( const Node & n1, const Node & n2 ) {
    std::set < Boundary * > common;
    intersectionSet( common, n1.boundSet(), n2.boundSet() );
    return findBoundary_( common );
}

Boundary * findBoundary( const Node & n1, const Node & n2, const Node & n3 ) {
    std::set < Boundary * > common;
    intersectionSet( common, n1.boundSet(), n2.boundSet(), n3.boundSet() );
    return findBoundary_( common );
}

Boundary * findBoundary( const std::vector < Node * > & n ) {
    std::vector < std::set< Boundary * > > bs;
    for ( uint i = 0; i < n.size(); i ++ ) bs.push_back( n[ i ]->boundSet() );

    std::set < Boundary * > common;
    intersectionSet( common, bs );

    if ( common.size() == 1 ) {
        return *common.begin();
    } else {
        if ( common.size() > 1 ){
            std::cerr << WHERE_AM_I << " pls. check, this should not happen. there is zero ore more then one cell defined." <<
            common.size() << std::endl;
            return *common.begin();
        }
    }
    return NULL;
}

Cell * findCommonCell( const std::vector < Node * > & n, bool warn ) {
  //** search the cell[s] which is[are] the intersection of all cells in n
    std::vector < std::set< Cell * > > bs;
    for ( uint i = 0; i < n.size(); i ++ ) bs.push_back( n[ i ]->cellSet() );

    std::set < Cell * > common;
    intersectionSet( common, bs );

    if ( common.size() == 1 ) {
        return *common.begin();
    } else {
        if ( common.size() > 1 ){
            if ( warn ){
                std::cerr << WHERE_AM_I << " pls. check, this should not happen. there is more then one cell defined." <<
                common.size() << std::endl;
            }
            return *common.begin();
        }
    }
    return NULL;
}

std::ostream & operator << ( std::ostream & str, const MeshEntity & e ){
    str << "MeshEntity " << &e << " rtti: " << e.rtti() << " id: " << e.id() << " rtti: " << e.rtti() << "\tN: " ;
    for ( uint i = 0; i < e.nodeCount(); i ++ ) str << e.node( i ).id() << " " ;
    return str;
}

std::ostream & operator << ( std::ostream & str, const Boundary & e ){
    str << "Boundary " << &e << " rtti: " << e.rtti() << " id: " << e.id() << "\tN: " ;
    for ( uint i = 0; i < e.nodeCount(); i ++ ) str << e.node( i ).id() << " " ;
    str << " marker: " << e.marker();
    return str;
}

std::ostream & operator << ( std::ostream & str, const Edge & e ){
    str << "Edge " << &e << " id: " << e.id() << "\t"
      << e.node(0).id() << " " << e.node(1).id()
      << " marker: " << e.marker();
  return str;
}

std::ostream & operator << ( std::ostream & str, const TriangleFace & t ){
  str << "TriangleFace " << &t << " id: " << t.id() << "\t"
      << t.node(0).id() << " " << t.node(1).id() << " " << t.node(2).id()
      << " attribute: " << t.marker();
  return str;
}

std::ostream & operator << ( std::ostream & str, const EdgeCell & c ){
  str << "EdgeCell " << &c << " id: " << c.id() << "\tN: ";
  for ( uint i = 0; i < c.nodeCount(); i ++ ) str << c.node( i ).id() << " " ;
  str << " attribute: " << c.attribute();
  return str;
}

std::ostream & operator << ( std::ostream & str, const Cell & c ){
  str << "Cell " << &c << " id: " << c.id() << "\tN: ";
  for ( uint i = 0; i < c.nodeCount(); i ++ ) str << c.node( i ).id() << " " ;
  str << " attribute: " << c.attribute();
  return str;
}

std::ostream & operator << ( std::ostream & str, const Triangle & t ){
  str << "Triangle " << &t << " id: " << t.id() << "\t"
      << t.node(0).id() << " " << t.node(1).id() << " " << t.node(2).id()
      << " attribute: " << t.attribute();
  return str;
}

std::ostream & operator << ( std::ostream & str, const Quadrangle & t ){
  str << "Quadrangle " << &t << " id: " << t.id() << "\t"
      << t.node(0).id() << " " << t.node(1).id() << " " << t.node(2).id()  << " " << t.node(3).id()
      << " attribute: " << t.attribute();
  return str;
}

std::ostream & operator << ( std::ostream & str, const Tetrahedron & t ){
    str << "Tetrahedron " << &t << " id: " << t.id() << "\t"
        << t.node(0).id() << " " << t.node(1).id() << " " << t.node(2).id()  << " " << t.node(3).id()
        << " attribute: " << t.attribute();
  return str;
}

std::ostream & operator << ( std::ostream & str, const Hexahedron & t ){
    str << "Hexahedron " << &t << " id: " << t.id() << "\t"
        << t.node( 0 ).id() << " " << t.node( 1 ).id() << " " << t.node( 2 ).id() << " "
        << t.node( 3 ).id() << " " << t.node( 4 ).id() << " " << t.node( 5 ).id() << " "
        << t.node( 6 ).id() << " " << t.node( 7 ).id()
        << " attribute: " << t.attribute();
    return str;
}

std::ostream & operator << ( std::ostream & str, const TriPrism & t ){
    str << "TrianglePrism" << &t << " id: " << t.id() << "\t"
        << t.node( 0 ).id() << " " << t.node( 1 ).id() << " " << t.node( 2 ).id() << " "
        << t.node( 3 ).id() << " " << t.node( 4 ).id() << " " << t.node( 5 ).id() << " "
        << " attribute: " << t.attribute();
    return str;
}

RVector3 MeshEntity::center() const {
    return shape_->center();
}

void MeshEntity::fillShape_( ){
    if ( shape_ ){
        if ( shape_->nodeCount() > this->nodeCount() ){
            std::cerr << WHERE_AM_I << " not enough nodes to fill shape: " << shape_->rtti()
                        << " " << this->rtti() << " " << shape_->nodeCount()
                        << " " << this->nodeCount() << std::endl;
        } else {
            for ( uint i = 0; i < shape_->nodeCount(); i ++ ) shape_->setNode( i, node( i ) );
        }
    }
}

void MeshEntity::setNodes_( std::vector < Node * > & nodes ){
    if ( nodes.size() > 0 ){
        if ( nodeVector_.size() != nodes.size() ) nodeVector_.resize( nodes.size() );
        std::copy( nodes.begin(), nodes.end(), &nodeVector_[ 0 ] );
    } else {
        std::cerr << WHERE_AM_I << " not enough nodes to fill meshEntity " << std::endl;
    }
}

void MeshEntity::deRegisterNodes_(){
}

std::pair< IndexArray, RVector > MeshEntity::shapeFunctions( const RVector3 & pos ) const {
    RVector f;
    shapeFunctionsL( shape().coordinates( pos ), f );
    IndexArray idVec( nodeCount() );
    RVector sF( nodeCount() );

    for ( uint i = 0; i < nodeCount(); i ++ ) {
        idVec[ i ] = node( i ).id();
        sF[ i ] = f[ i ];
    }
    return std::pair< IndexArray, RVector >( idVec, sF );
}

void MeshEntity::shapeFunctionsDeriveL( const RVector3 & coord, uint dim, RVector & funct ) const {
    funct.resize( nodeCount() );
    switch ( dimension() ){
        case 1:
        {
            double drd = shape_->deriveCoordinates( 0, dim );
            funct = drd * deriveNdL( coord, 0 );
        } break;
        case 2:
        {
            double drd = shape_->deriveCoordinates( 0, dim );
            double dsd = shape_->deriveCoordinates( 1, dim );
            funct = drd * deriveNdL( coord, 0 ) + dsd * deriveNdL( coord, 1 );
        } break;
        case 3:
        {
            double drd = shape_->deriveCoordinates( 0, dim );
            double dsd = shape_->deriveCoordinates( 1, dim );
            double dtd = shape_->deriveCoordinates( 2, dim );
            funct = drd * deriveNdL( coord, 0 ) + dsd * deriveNdL( coord, 1 ) + dtd * deriveNdL( coord, 2 );
        } break;
    };

}


std::pair< IndexArray, std::vector < RVector3 > > MeshEntity::shapeFunctionsDerive( const RVector3 & pos ) const {
    RVector3 coords = shape().coordinates( pos );

    IndexArray idVec( nodeCount() );
    for ( uint i = 0; i < nodeCount(); i ++ ) idVec[ i ] = node( i ).id();

    std::vector< RVector3 > dSF( nodeCount(), RVector3( 0.0, 0.0 ) );
    RVector funct( nodeCount() );

    for ( uint dim = 0; dim < dimension(); dim ++ ){
        shapeFunctionsDeriveL( coords, dim, funct );

        for ( uint i = 0; i < nodeCount(); i ++ ) {
            dSF[ i ][ dim ] = funct[ i ];
        }
    }

    return std::pair< IndexArray, std::vector < RVector3 > >( idVec, dSF );
}

double MeshEntity::interpolate( const RVector3 & pos, const RVector & data ) const {
   return interpolate_( shapeFunctions( pos ), data );
}

double MeshEntity::interpolate_( const std::pair< IndexArray, RVector > & sF, const RVector & data ) const {
    return sum( data( sF.first ) * sF.second );
}

RVector3 MeshEntity::grad( const RVector3 & pos, const RVector & data ) const {
   return grad_( shapeFunctionsDerive( pos ), data );
}

RVector3 MeshEntity::grad_( const std::pair< IndexArray, std::vector < RVector3 > > & sF, const RVector & data ) const {
    RVector3 res(0.0, 0.0, 0.0);
    for ( uint i = 0; i < 3; i ++ ){
        //std::cout << i << std::endl;;
        //for ( std::map< uint, double >::const_iterator it = sF[ i ].begin(); it != sF[ i ].end(); it ++  ){
          //  std::cout << "   " << it->first << "=" << data[ it->first ] << " " << it->second<< std::endl;;
        for ( uint j = 0; j < nodeCount(); j ++ ){
            res[ i ] += data[ sF.first[ j ] ] * sF.second[ j ][ i ];
        }
    }
    //std::cout << res << std::endl;
    return res;
}


Node * Cell::oppositeTo( const Boundary & bound ){
    THROW_TO_IMPL
    //** maybee obsolete pls check

//     std::set < Node * > nodes;
//     for ( int i = 0; i < 3; i++ ) nodes.insert( &node( i ) );
//
//     for ( uint i = 0; i < bound.nodeCount(); i ++ ) nodes.erase( & bound.node( i ) );
//
//     if ( nodes.size() == 1 ){
//         return ( *nodes.begin() );
//     } else {
//         std::cerr << bound << std::endl << *this << std::endl;
//         std::cerr << WHERE_AM_I << " this triangle does not contain the edge" << std::endl;
//     }
     return NULL;
}

double Cell::jacobianDeterminant() const {
    return shape_->jacobianDeterminant();
}

void Cell::cleanNeighbourInfos( ){
    for ( uint i = 0; i < this->neighbourCellCount(); i ++ ){
        neighbourCells_[ i ] = NULL;
    }
}

void Cell::findNeighbourCell( uint i ){
    std::vector < Node * > n( boundaryNodes( i ) );

    std::set < Cell *> common;
    std::set < Cell *> commonTmp;

    if ( n.size() > 1 ) intersectionSet( common, n[ 0 ]->cellSet(), n[ 1 ]->cellSet() );
    
    for ( uint j = 2; j < n.size(); j ++ ){
        commonTmp = common;
        intersectionSet( common, commonTmp, n[ j ]->cellSet() );
    }

    common.erase( this );
    if ( common.size() == 1 ) neighbourCells_[ i ] = *common.begin(); else neighbourCells_[ i ] = NULL;
}

void Cell::registerNodes_( ){
    for ( uint i = 0; i < this->nodeCount(); i ++ ){
        if ( nodeVector_[ i ] ) nodeVector_[ i ]->insertCell( *this );
    }
}

void Cell::deRegisterNodes_(){
    for ( uint i = 0; i < this->nodeCount(); i ++ ){
        if ( nodeVector_[ i ] ) nodeVector_[ i ]->eraseCell( *this );
    }
}

void Boundary::registerNodes_(){
    for ( uint i = 0; i < this->nodeCount(); i ++ ){
        if ( nodeVector_[ i ] ) nodeVector_[ i ]->insertBoundary( *this );
    }
}

void Boundary::deRegisterNodes_(){
    for ( uint i = 0; i < this->nodeCount(); i ++ ){
        if ( nodeVector_[ i ] ) nodeVector_[ i ]->eraseBoundary( *this );
    }
}

RVector3 Boundary::norm( ) const {
    return shape_->norm();
}

NodeBoundary::NodeBoundary( std::vector < Node * > & nodes ) : Boundary( nodes ){
  shape_ = new NodeShape();
  fillShape_( );
}

NodeBoundary::NodeBoundary( Node & n1 ) : Boundary( ) {
  shape_ = new NodeShape();
  setNodes( n1, false );
}

NodeBoundary::~NodeBoundary(){
  delete shape_;
}

void NodeBoundary::setNodes( Node & n1, bool changed ){
  if ( changed ) deRegisterNodes_( );

  std::vector < Node * > nodes;
  nodes.push_back( &n1 );
  setNodes_( nodes );
  registerNodes_( );
  fillShape_( );
}

Edge::Edge( std::vector < Node * > & nodes ) : Boundary( nodes ){
  shape_ = new EdgeShape();
  fillShape_( );
}

Edge::Edge( Node & n1, Node & n2 ){
  shape_ = new EdgeShape();
  setNodes( n1, n2, false );
}

Edge::~Edge(){
    delete shape_;
}

void Edge::setNodes( Node & n1, Node & n2, bool changed ){
  if ( (&n1 == &n2) ){
    std::cerr << WHERE << " Edge nodes not valid " << n1 << " " << n2 << " " << std::endl;
    exit( EXIT_MESH_NO_ELEMENT );
  }
  if ( changed ) deRegisterNodes_( );

  std::vector < Node * > nodes;
  nodes.push_back( &n1 ); nodes.push_back( &n2 );
  setNodes_( nodes );
  registerNodes_( );
  fillShape_( );
}

int Edge::swap(){
  if ( (marker_ != 0 ) || ( leftCell_ == NULL ) || ( rightCell_ == NULL ) ) {
    //    cout<< WHERE_AM_I << " " << marker_ << " " << leftCell_  << " " << rightCell_ << std::endl;
    return 0;
  } else if ( ( leftCell_->rtti() == MESH_TRIANGLE_RTTI ) && ( rightCell_->rtti() == MESH_TRIANGLE_RTTI ) ) {

    Node * oA = &node( 0 );
    Node * oB = &node( 1 );
    Triangle *left = dynamic_cast< Triangle * >( leftCell_ );
    Triangle *right = dynamic_cast< Triangle * >( rightCell_ );

    Node * oL = left->oppositeTo( * this );
    Node * oR = right->oppositeTo( * this );

    if ( oL == NULL || oR == NULL ){
      std::cout << *this << std::endl
		<< left << std::endl
		<< right << std::endl;
      if ( oL != NULL ) std::cout << "oL " << oL->id() << std::endl;
      if ( oR != NULL ) std::cout << "oR " << oR->id() << std::endl;
      exit(0);
    }

    //** swap only when the resulting triangles have same sign
    //** abort swapping for concav domain, posible check by angles
    if ( sign( jacobianDetXY( oL->pos(), oR->pos(), oB->pos() ) ) !=
	       sign( jacobianDetXY( oL->pos(), oA->pos(), oR->pos() ) ) ){
      return 0;
    }

    right->setNodes( * oL, * oA, * oR );

    setNodes( * oL, * oR );
    //    cout << "\t" << *this << std::endl;

    if ( leftCell_ == rightCell_ ){
      std::cerr << WHERE << " Edge " << id() << " wrong swapped " << std::endl;
      std::cerr << "LeftElement: " << left->id()
	   << "; RigthElement: " << right->id() << std::endl;
      std::cerr << "NodeA: " << oA->id() << ", NodeB: " << oB->id()
	   << ", NodeL: " << oL->id() << ", NodeR: " << oR->id() << std::endl;
      return 0;
    }

    left->setNodes( * oL, * oR, * oB );
    right->setNodes( * oL, * oA, * oR );

    return 1;
  }
  return 0;
}

void Edge::shapeFunctionsL( const RVector3 & L, RVector & funct ) const {
    return Edge_shapeFunctionsL( L, funct );
}

// double Edge::interpolate( const RVector3 & queryPos, const RVector & data ) const {
//   double radius = node( 0 ).pos().distance( node( 1 ).pos() );
//   double rmin = node( 0 ).pos().distance( queryPos );
//
//   double v0 = data[ node( 0 ).id() ];
//   double v1 = data[ node( 1 ).id() ];
//   return v0 + ( ( v1 - v0 ) / radius ) * rmin;
// }

Edge3::Edge3( std::vector < Node * > & nodes ) : Edge( nodes ){
}

Edge3::~Edge3(){
}

void Edge3::shapeFunctionsL( const RVector3 & L, RVector & funct ) const {
    return Edge3_shapeFunctionsL( L, funct );
}

TriangleFace::TriangleFace( std::vector < Node * > & nodes ) : Boundary( nodes ){
  shape_ = new TriangleShape();
  fillShape_( );
}

TriangleFace::TriangleFace( Node & n1, Node & n2, Node & n3 ){
  shape_ = new TriangleShape();
  setNodes( n1, n2, n3, false  );
}

TriangleFace::~TriangleFace(){
  delete shape_;
}

void TriangleFace::setNodes( Node & n1, Node & n2, Node & n3, bool changed  ){
  if ( (&n1 == &n2) || (&n1 == &n3) || (&n2 == &n3)){
    std::cerr << WHERE << " TriangleFace nodes not valid " << n1 << " " << n2 << " " <<  n3 << std::endl;
    exit( EXIT_MESH_NO_ELEMENT );
  }
  if ( changed ) deRegisterNodes_( );

  std::vector < Node * > nodes;
  nodes.push_back( &n1 ); nodes.push_back( &n2 ); nodes.push_back( &n3 );
  setNodes_( nodes );
  registerNodes_( );
  fillShape_( );
}

void TriangleFace::shapeFunctionsL( const RVector3 & L, RVector & funct ) const {
    return Triangle_shapeFunctionsL( L, funct );
}

TriangleFace6::TriangleFace6( std::vector < Node * > & nodes ) : TriangleFace( nodes ){
}

TriangleFace6::~TriangleFace6(){
}

void TriangleFace6::shapeFunctionsL( const RVector3 & L, RVector & funct ) const {
    return Triangle6_shapeFunctionsL( L, funct );
}

QuadrangleFace::QuadrangleFace( std::vector < Node * > & nodes ) : Boundary( nodes ){
  shape_ = new QuadrangleShape();
  fillShape_( );
}

QuadrangleFace::QuadrangleFace( Node & n1, Node & n2, Node & n3, Node & n4 ){
  shape_ = new QuadrangleShape();
  setNodes( n1, n2, n3, n4, false  );
}

QuadrangleFace::~QuadrangleFace(){
  delete shape_;
}

void QuadrangleFace::setNodes( Node & n1, Node & n2, Node & n3, Node & n4, bool changed  ){
  if ( (&n1 == &n2) || (&n1 == &n3) || (&n2 == &n3)){
    std::cerr << WHERE << " QuadrangleFace nodes not valid " << n1 << " "
            << n2 << " " <<  n3 << " " << n4<< std::endl;
    exit( EXIT_MESH_NO_ELEMENT );
  }
  if ( changed ) deRegisterNodes_( );

  std::vector < Node * > nodes;
  nodes.push_back( &n1 ); nodes.push_back( &n2 );
  nodes.push_back( &n3 ); nodes.push_back( &n4 );
  setNodes_( nodes );
  registerNodes_( );
  fillShape_( );
}

EdgeCell::EdgeCell( std::vector < Node * > & nodes ) : Cell( nodes ){
    shape_ = new EdgeShape();
    fillShape_( );
    neighbourCells_.resize( this->neighbourCellCount(), NULL );
}

EdgeCell::EdgeCell( Node & n1, Node & n2 ) : Cell() {
    shape_ = new EdgeShape();
    setNodes( n1, n2, false  );
    neighbourCells_.resize( this->neighbourCellCount(), NULL );
}

EdgeCell::~EdgeCell(){
    delete shape_;
}

void EdgeCell::setNodes( Node & n1, Node & n2, bool changed ){
    if ( changed ) deRegisterNodes_( );

    std::vector < Node * > nodes;
    nodes.push_back( &n1 ); nodes.push_back( &n2 );
    setNodes_( nodes );
    registerNodes_( );
    fillShape_( );
}

void EdgeCell::findNeighbourCell( uint id ){
    //     case: 0 - 1
    //     case: 1 - 0

    std::set < Cell * > common( nodeVector_[ (id + 1 )%2 ]->cellSet() );
    common.erase( this );
    if ( common.size() == 1 ) neighbourCells_[ id ] = *common.begin(); else neighbourCells_[ id ] = NULL;
}

void EdgeCell::shapeFunctionsL( const RVector3 & L, RVector & funct ) const{
    return Edge_shapeFunctionsL( L, funct );
}

RVector EdgeCell::deriveNdL( const RVector3 & coord, uint dim ) const {
    RVector f( nodeCount() );
    f[ 0 ] = - 1.0;
    f[ 1 ] =   1.0;
    return f;
}

Edge3Cell::Edge3Cell( std::vector < Node * > & nodes ) : EdgeCell( nodes ){
}

Edge3Cell::~Edge3Cell(){
}

void Edge3Cell::shapeFunctionsL( const RVector3 & L, RVector & funct ) const {
    return Edge3_shapeFunctionsL( L, funct );
}

RVector Edge3Cell::deriveNdL( const RVector3 & coord, uint dim ) const {
    RVector f( nodeCount() );
    double r = coord[ 0 ];
    f[ 0 ] = 4.0 * r - 3.0;
    f[ 1 ] = 4.0 * r - 1.0;
    f[ 2 ] = 4.0 - 8.0 * r;
    return f;
}

Triangle::Triangle( std::vector < Node * > & nodes ) : Cell( nodes ){
    shape_ = new TriangleShape();
    fillShape_( );
    neighbourCells_.resize( this->neighbourCellCount(), NULL );
}

Triangle::Triangle( Node & n1, Node & n2, Node & n3 ): Cell(){
    shape_ = new TriangleShape();
    setNodes( n1, n2, n3, false );
    neighbourCells_.resize( this->neighbourCellCount(), NULL );
}

Triangle::~Triangle(){
//     std::cout << " delete Triangle shape_ " << std::endl;
    delete shape_;
}

void Triangle::setNodes( Node & n1, Node & n2, Node & n3, bool changed ){
    if ( (&n1 == &n2) || (&n1 == &n3) || (&n2 == &n3 ) ){
        std::cerr << WHERE << " Triangle nodes not valid " << n1 << " " << n2 << " " << n3 << std::endl;
        exit( EXIT_MESH_NO_ELEMENT );
    }

    if ( changed ) deRegisterNodes_( );

    std::vector < Node * > nodes;
    nodes.push_back( &n1 ); nodes.push_back( &n2 ); nodes.push_back( &n3 );
    setNodes_( nodes );
    registerNodes_( );
    fillShape_( );

//  double J = this->jacobianDeterminant();

//   //** switch direction to force positiv jacobian determinant
//   if ( J < 0 ) {
//     shape_->setNode( 1, n3 );
//     shape_->setNode( 2, n2 );
//     J *= -1.0;
//   }
//
//   for ( int i = 0; i < 3; i++ ){
//     Boundary * bound = findBoundary( node( i ), node( (i+1)%3 ) );
//
//     if ( bound != NULL ) {
//       if ( &bound->node( 0 ) == &node( i ) ) {
// 	bound->setLeftCell( *this );
//       } else {
// 	bound->setRightCell( *this );
//       }
//     }
//   }
}

void Triangle::findNeighbourCell( uint i ){
    std::set < Cell * > common;
    //** cell oposite to node(i)
    intersectionSet( common, nodeVector_[ ( i + 1 )%3 ]->cellSet(), nodeVector_[ ( i + 2 )%3 ]->cellSet() );
    common.erase( this );
    if ( common.size() == 1 ) neighbourCells_[ i ] = *common.begin(); else neighbourCells_[ i ] = NULL;
}

void Triangle::shapeFunctionsL( const RVector3 & L, RVector & funct ) const {
    return Triangle_shapeFunctionsL( L, funct );
}

// dNdL( Nr, pos ), return for each N
RVector Triangle::deriveNdL( const RVector3 & coord, uint dim ) const {
    RVector f( nodeCount() );
    switch ( dim ){
        case 0: // dNdr
            f[ 0 ] = - 1.0;
            f[ 1 ] =   1.0;
            f[ 2 ] =   0.0;
            return f;
        case 1:// dNds
            f[ 0 ] = - 1.0;
            f[ 1 ] =   0.0;
            f[ 2 ] =   1.0;
            return f;
    }
    throwLengthError( 1, WHERE_AM_I + " undefined values " + toStr( coord ) + " " + toStr( dim ) );
    return f;
}

Triangle6::Triangle6( std::vector < Node * > & nodes ) : Triangle( nodes ){
}

Triangle6::~Triangle6(){
}

void Triangle6::shapeFunctionsL( const RVector3 & coord, RVector & funct ) const {
    return Triangle6_shapeFunctionsL( coord, funct );
}

// void Triangle6::shapeFunctionsDeriveL( const RVector3 & L, uint coord, RVector & funct ) const {
//     funct.resize( nodeCount() );
//     double dL2d = shape_->deriveCoordinates( 1, coord );
//     double dL3d = shape_->deriveCoordinates( 2, coord );
//
//     double L2 = L[ 0 ];
//     double L3 = L[ 1 ];
//     funct[ 0 ] = dL3d * ( 4.0 * L3 + 4.0 * L2 - 3.0)   + dL2d * ( 4.0 * L3 + 4.0 * L2 - 3.0 );
//     funct[ 1 ] =                                         dL2d * ( 4.0 * L2 - 1.0 );
//     funct[ 2 ] = dL3d * ( 4.0 * L3 - 1.0 );
//     funct[ 3 ] = dL3d * -4.0 * L2                      + dL2d * ( -4.0 * L3 - 8.0 * L2 + 4.0 );
//     funct[ 4 ] = dL3d *  4.0 * L2                      + dL2d *  4.0 * L3;
//     funct[ 5 ] = dL3d * (- 8.0 * L3 - 4.0 * L2 + 4.0 ) + dL2d * -4.0 * L3;
//
//     return;
//     funct = dL3d * deriveNdL( L, 0 ) + dL2d * deriveNdL( L, 1 );
// }

RVector Triangle6::deriveNdL( const RVector3 & coord, uint dim ) const {
    RVector f( nodeCount() );

    double r = coord[ 0 ];
    double s = coord[ 1 ];
    switch ( dim ){
        case 0: //dNdr
            f[ 0 ] =  4.0 * r + 4.0 * s - 3.0;
            f[ 1 ] =  4.0 * r - 1.0;
            f[ 2 ] =  0.0;
            f[ 3 ] = -4.0 * s - 8.0 * r + 4.0;
            f[ 4 ] =  4.0 * s;
            f[ 5 ] = -4.0 * s;
            return f;
        break;
        case 1://dNds
            f[ 0 ] =  4.0 * s + 4.0 * r - 3.0;
            f[ 1 ] =  0.0;
            f[ 2 ] =  4.0 * s - 1.0;
            f[ 3 ] = -4.0 * r;
            f[ 4 ] =  4.0 * r;
            f[ 5 ] = -8.0 * s - 4.0 * r + 4.0;
            return f;
    }
    throwLengthError( 1, WHERE_AM_I + " undefined values " + toStr( coord ) + " " + toStr( dim ) );
    return f;
}

Quadrangle::Quadrangle( std::vector < Node * > & nodes ) : Cell( nodes ){
  shape_ = new QuadrangleShape();
  fillShape_( );
  neighbourCells_.resize( this->neighbourCellCount(), NULL );
}

Quadrangle::Quadrangle( Node & n1, Node & n2, Node & n3, Node & n4 ): Cell(){
  shape_ = new QuadrangleShape();
  setNodes( n1, n2, n3, n4, false  );
  neighbourCells_.resize( this->neighbourCellCount(), NULL );
}

Quadrangle::~Quadrangle(){
    delete shape_;
}

void Quadrangle::setNodes( Node & n1, Node & n2, Node & n3, Node & n4, bool changed ){
  if ( changed ) deRegisterNodes_( );

  std::vector < Node * > nodes;
  nodes.push_back( &n1 ); nodes.push_back( &n2 );
  nodes.push_back( &n3 ); nodes.push_back( &n4 );
  setNodes_( nodes );
  registerNodes_( );
  fillShape_( );
}

void Quadrangle::findNeighbourCell( uint id ){
    std::set < Cell * > common;

    intersectionSet( common, nodeVector_[ ( id + 1 )%4 ]->cellSet(), nodeVector_[ ( id + 2 )%4 ]->cellSet() );
    common.erase( this );
    if ( common.size() == 1 ) neighbourCells_[ id ] = *common.begin(); else neighbourCells_[ id ] = NULL;

    //** start check;
//     if ( this->id() == ){
//
//         int pIdx = 0;
//         shape_->touch1( neighbourCells_[ id ]->center(), true, pIdx );
//         if ( pIdx != id ){
//             std::cout << "something goes wrong here idx:" << id << " pIdx = " << pIdx << " "
//                     <<  id - pIdx<< std::endl;
//             std::cout << "orign:" << *this << std::endl;
//             std::cout << "test:" << *neighbourCells_[ id ] << std::endl;
//             for ( uint i = 0; i <4 ; i ++ ){
//                 std::cout << this->node( i )<< std::endl;
//             }
//             for ( uint i = 0; i <4 ; i ++ ){
//                 common.clear();
//                 intersectionSet( common, nodeVector_[ ( i ) ]->cellSet(),
//                                         nodeVector_[ ( i + 1 )%4 ]->cellSet() );
//                 common.erase( this );
//                 if ( common.size() == 1 ) std::cout << i << " "<< **common.begin() << std::endl;
//             }
//
//             exit(1);
//         }
//     }
}

void Quadrangle::shapeFunctionsL( const RVector3 & L, RVector & funct ) const{
    funct.resize( nodeCount() );
    double r = L[ 0 ];
    double s = L[ 1 ];
    funct[ 0 ] = ( 1.0 - r ) * ( 1.0 - s );
    funct[ 1 ] = r * ( 1.0 - s );
    funct[ 2 ] = r * s;
    funct[ 3 ] = ( 1.0 - r ) * s;
}

//     funct.resize( nodeCount() );
//     double drd = shape_->deriveCoordinates( 0, coord );
//     double dsd = shape_->deriveCoordinates( 1, coord );
//     double r = L[ 0 ];
//     double s = L[ 1 ];
//   //  std::cout << "dLdxy: " << r << " "<< s << " " << drd << " "<< dsd << std::endl;
// //     r = 0.5;
// //     s = 0.5;
//
//     funct[ 0 ] =  drd * ( s - 1.0 ) + dsd * ( r - 1.0 );
//     funct[ 1 ] =  drd * ( 1.0 - s ) - dsd * r;
//     funct[ 2 ] =  drd * s           + dsd * r;
//     funct[ 3 ] = -drd * s           + dsd * ( 1.0 - r );
//}

RVector Quadrangle::deriveNdL( const RVector3 & coord, uint dim ) const {
    RVector f( nodeCount() );
    double r = coord[ 0 ];
    double s = coord[ 1 ];
    switch ( dim ){
        case 0: //dNdr
            f[ 0 ] =  s - 1.0;
            f[ 1 ] =  1.0 - s;
            f[ 2 ] =  s;
            f[ 3 ] = -s;
            return f;
        break;
        case 1://dNds
            f[ 0 ] =  r - 1.0;
            f[ 1 ] = -r;
            f[ 2 ] =  r;
            f[ 3 ] =  1.0 - r;
            return f;
    }
    throwLengthError( 1, WHERE_AM_I + " undefined values " + toStr( coord ) + " " + toStr( dim ) );
    return f;
}

Quadrangle8::Quadrangle8( std::vector < Node * > & nodes ) : Quadrangle( nodes ){
}

Quadrangle8::~Quadrangle8(){
}

void Quadrangle8::shapeFunctionsL( const RVector3 & L, RVector & funct ) const {
    funct.resize( nodeCount() );
    double r = L[ 0 ];
    double s = L[ 1 ];

    funct[ 0 ] = ( 1.0 - r ) * ( 1.0 - s ) * ( 1.0 - 2.0 * r - 2.0 * s );
    funct[ 1 ] = -r * ( 1.0 - s )* ( 1.0 - 2.0 * r + 2.0 * s );
    funct[ 2 ] = -r * s * ( 3.0 - 2.0 * r - 2.0 * s );
    funct[ 3 ] = -s * ( 1.0 - r ) * ( 1.0 + 2.0 * r - 2.0 * s );
    funct[ 4 ] = 4.0 * r * ( 1.0 - r ) * ( 1.0 - s );
    funct[ 5 ] = 4.0 * r * s * ( 1.0 - s );
    funct[ 6 ] = 4.0 * r * s * ( 1.0 - r );
    funct[ 7 ] = 4.0 * s * ( 1.0 - r ) * ( 1.0 - s );
}

// void Quadrangle8::shapeFunctionsDeriveL( const RVector3 & L, uint coord, RVector & funct ) const{
//     double drd = shape_->deriveCoordinates( 0, coord );
//     double dsd = shape_->deriveCoordinates( 1, coord );
//     double r = L[ 0 ];
//     double s = L[ 1 ];
//
//     funct[ 0 ] = drd * ( - 2.0 * s * s  + ( 5.0 - 4.0 * r ) * s + 4.0 * r - 3.0 ) +
//                  dsd * ( ( 4.0 - 4.0 * r ) * s - 2.0 * r * r  + 5.0 * r - 3.0 );
//     funct[ 1 ] = drd * ( 2.0 * s * s + (- 4.0 * r - 1.0 ) * s + 4.0 * r - 1.0 ) +
//                  dsd * ( 4.0 * r * s - 2.0 * r * r - r);
//     funct[ 2 ] = drd * ( 2.0 * s * s + ( 4.0 * r - 3.0 ) * s ) +
//                  dsd * ( 4.0 * r * s + 2.0 * r * r - 3.0 * r );
//     funct[ 3 ] = drd * ( ( 4.0 * r - 1.0 ) * s - 2.0 * s * s ) +
//                  dsd * ( ( 4.0 - 4.0 * r ) * s + 2.0 * r * r - r - 1.0 );
//     funct[ 4 ] = drd * ( ( 8.0 * r - 4.0 ) * s - 8.0 * r + 4.0 ) +
//                  dsd * ( 4.0 * r * r - 4.0 * r );
//     funct[ 5 ] = drd * ( 4.0 * s - 4.0 * s * s ) +
//                  dsd * ( 4.0 * r - 8.0 * r * s );
//     funct[ 6 ] = drd * ( 4.0 - 8.0 * r ) * s +
//                  dsd * ( 4.0 * r - 4.0 * r * r );
//     funct[ 7 ] = drd * ( 4.0 * s * s - 4.0 * s ) +
//                  dsd * ( ( 8.0 * r - 8.0 ) * s - 4.0 * r + 4.0 );
// }

RVector Quadrangle8::deriveNdL( const RVector3 & coord, uint dim ) const {
    RVector f( nodeCount() );
    double r = coord[ 0 ];
    double s = coord[ 1 ];
    switch ( dim ){
        case 0: //dNdr
            f[ 0 ] = ( - 2.0 * s * s  + ( 5.0 - 4.0 * r ) * s + 4.0 * r - 3.0 );
            f[ 1 ] = ( 2.0 * s * s + (- 4.0 * r - 1.0 ) * s + 4.0 * r - 1.0 );
            f[ 2 ] = 2.0 * s * s + ( 4.0 * r - 3.0 ) * s;
            f[ 3 ] = ( 4.0 * r - 1.0 ) * s - 2.0 * s * s;
            f[ 4 ] = ( 8.0 * r - 4.0 ) * s - 8.0 * r + 4.0;
            f[ 5 ] = 4.0 * s - 4.0 * s * s;
            f[ 6 ] = ( 4.0 - 8.0 * r ) * s;
            f[ 7 ] = 4.0 * s * s - 4.0 * s;
            return f;
        case 1://dNds
            f[ 0 ] = ( 4.0 - 4.0 * r ) * s - 2.0 * r * r  + 5.0 * r - 3.0 ;
            f[ 1 ] = 4.0 * r * s - 2.0 * r * r - r;
            f[ 2 ] = 4.0 * r * s + 2.0 * r * r - 3.0 * r;
            f[ 3 ] = ( 4.0 - 4.0 * r ) * s + 2.0 * r * r - r - 1.0;
            f[ 4 ] = 4.0 * r * r - 4.0 * r;
            f[ 5 ] = 4.0 * r - 8.0 * r * s;
            f[ 6 ] = 4.0 * r - 4.0 * r * r;
            f[ 7 ] = ( 8.0 * r - 8.0 ) * s - 4.0 * r + 4.0 ;
            return f;
    }
    throwLengthError( 1, WHERE_AM_I + " undefined values " + toStr( coord ) + " " + toStr( dim ) );
    return f;
}

Tetrahedron::Tetrahedron( std::vector < Node * > & nodes ) : Cell( nodes ){
  shape_ = new TetrahedronShape();
  fillShape_( );
  neighbourCells_.resize( this->neighbourCellCount(), NULL );
}

Tetrahedron::Tetrahedron( Node & n1, Node & n2, Node & n3, Node & n4 ) : Cell() {
  shape_ = new TetrahedronShape();
  setNodes( n1, n2, n3, n4, false  );
  neighbourCells_.resize( this->neighbourCellCount(), NULL );
}

Tetrahedron::~Tetrahedron(){
  delete shape_;
}

void Tetrahedron::setNodes( Node & n1, Node & n2, Node & n3, Node & n4, bool changed ){
  if ( changed ) deRegisterNodes_( );

  std::vector < Node * > nodes;
  nodes.push_back( &n1 ); nodes.push_back( &n2 );
  nodes.push_back( &n3 ); nodes.push_back( &n4 );
  setNodes_( nodes );
  registerNodes_( );
  fillShape_( );
}

void Tetrahedron::findNeighbourCell( uint i ){
    std::set < Cell * > c1, common;
    intersectionSet( c1, nodeVector_[ ( i + 1 )%4 ]->cellSet(), nodeVector_[ ( i + 2 )%4 ]->cellSet() );
    intersectionSet( common, c1, nodeVector_[ ( i + 3 )%4 ]->cellSet() );

    common.erase( this );
    if ( common.size() == 1 ) neighbourCells_[ i ] = *common.begin(); else neighbourCells_[ i ] = NULL;
}

void Tetrahedron::shapeFunctionsL( const RVector3 & L, RVector & funct ) const{
    return Tetrangle_shapeFunctionsL( L, funct );
}

RVector Tetrahedron::deriveNdL( const RVector3 & coord, uint dim ) const {
    RVector f( nodeCount() );
    switch ( dim ){
        case 0: // dNdr
            f[ 0 ] = - 1.0;
            f[ 1 ] =   1.0;
            f[ 2 ] =   0.0;
            f[ 3 ] =   0.0;
            return f;
        case 1:// dNds
            f[ 0 ] = - 1.0;
            f[ 1 ] =   0.0;
            f[ 2 ] =   1.0;
            f[ 3 ] =   0.0;
            return f;
        case 2:// dNdt
            f[ 0 ] = - 1.0;
            f[ 1 ] =   0.0;
            f[ 2 ] =   0.0;
            f[ 3 ] =   1.0;
            return f;
    }
    throwLengthError( 1, WHERE_AM_I + " undefined values " + toStr( coord ) + " " + toStr( dim ) );
    return f;
}

Tetrahedron10::Tetrahedron10( std::vector < Node * > & nodes ) : Tetrahedron( nodes ){
}

Tetrahedron10::~Tetrahedron10(){
}

void Tetrahedron10::shapeFunctionsL( const RVector3 & L, RVector & funct ) const {
    funct.resize( nodeCount() );
    double L1 = 1.0 - L[ 0 ] - L[ 1 ] - L[ 2 ];
    double L2 = L[ 0 ];
    double L3 = L[ 1 ];
    double L4 = L[ 2 ];
    funct[ 0 ] = L1 * ( 2.0 * L1 - 1.0 );
    funct[ 1 ] = L2 * ( 2.0 * L2 - 1.0 );
    funct[ 2 ] = L3 * ( 2.0 * L3 - 1.0 );
    funct[ 3 ] = L4 * ( 2.0 * L4 - 1.0 );

    //*zienkiewitz numbering dcfemlib(meshconvert -p)
    funct[ 4 ] = 4.0 * L1 * L2;
    funct[ 5 ] = 4.0 * L1 * L3;
    funct[ 6 ] = 4.0 * L1 * L4;
    funct[ 7 ] = 4.0 * L2 * L3;
    funct[ 8 ] = 4.0 * L3 * L4;
    funct[ 9 ] = 4.0 * L4 * L2;

    //*lokal numbering
//     sF[ node( 4 ).id() ] = 4.0 * z1 * z2;
//     sF[ node( 5 ).id() ] = 4.0 * z2 * z3;
//     sF[ node( 6 ).id() ] = 4.0 * z3 * z1;
//     sF[ node( 7 ).id() ] = 4.0 * z1 * z4;
//     sF[ node( 8 ).id() ] = 4.0 * z2 * z4;
//     sF[ node( 9 ).id() ] = 4.0 * z3 * z4;
}

// void Tetrahedron10::shapeFunctionsDeriveL( const RVector3 & L, uint coord, RVector & funct ) const {
//     funct.resize( nodeCount() );
//
//     double dL2 = shape_->deriveCoordinates( 1, coord );
//     double dL3 = shape_->deriveCoordinates( 2, coord );
//     double dL4 = shape_->deriveCoordinates( 3, coord );
//     double dL1 = 1- dL2-dL3-dL4;
//     double L1 = 1.0 - L[ 0 ] - L[ 1 ] - L[ 2 ];
//     double L2 = L[ 0 ];
//     double L3 = L[ 1 ];
//     double L4 = L[ 2 ];
//
//     funct[ 0 ] = dL1 * ( 4.0 * L1 -1.0 );
//     funct[ 1 ] = dL2 * ( 4.0 * L2 -1.0 );
//     funct[ 2 ] = dL3 * ( 4.0 * L3 -1.0 );
//     funct[ 3 ] = dL4 * ( 4.0 * L4 -1.0 );
//
//     //*zienkiewitz numbering dcfemlib(meshconvert -p)
//     funct[ 4 ] = ( dL1 * L2 + dL2 * L1 ) * 4.0;
//     funct[ 5 ] = ( dL1 * L3 + dL3 * L1 ) * 4.0;
//     funct[ 6 ] = ( dL1 * L4 + dL4 * L1 ) * 4.0;
//     funct[ 7 ] = ( dL2 * L3 + dL3 * L2 ) * 4.0;
//     funct[ 8 ] = ( dL3 * L4 + dL4 * L3 ) * 4.0;
//     funct[ 9 ] = ( dL4 * L2 + dL2 * L4 ) * 4.0;
// //         ret.back()[ node( 4 ).id() ] = ( dL1 * L2 + dL2 * L1 ) * 4.0;
// //         ret.back()[ node( 5 ).id() ] = ( dL2 * L3 + dL3 * L2 ) * 4.0;
// //         ret.back()[ node( 6 ).id() ] = ( dL3 * L1 + dL1 * L3 ) * 4.0;
// //         ret.back()[ node( 7 ).id() ] = ( dL1 * L4 + dL4 * L1 ) * 4.0;
// //         ret.back()[ node( 8 ).id() ] = ( dL2 * L4 + dL4 * L2 ) * 4.0;
// //         ret.back()[ node( 9 ).id() ] = ( dL3 * L4 + dL4 * L3 ) * 4.0;
//
// }

RVector Tetrahedron10::deriveNdL( const RVector3 & coord, uint dim ) const {
    RVector f( nodeCount() );
    double r = coord[ 0 ];
    double s = coord[ 1 ];
    double t = coord[ 2 ];

    switch ( dim ){
        case 0: // dNdr
            f[ 0 ] =  4.0 * r + 4.0 * s + 4.0 * t - 3.0;
            f[ 1 ] =  4.0 * r - 1.0;
            f[ 2 ] =  0.0;
            f[ 3 ] =  0.0;
            f[ 4 ] = -8.0 * r - 4.0 * s - 4.0 * t + 4.0;
            f[ 5 ] = -4.0 * s;
            f[ 6 ] = -4.0 * t;
            f[ 7 ] =  4.0 * s;
            f[ 8 ] =  0.0;
            f[ 9 ] =  4.0 * t;
            return f;
        case 1:// dNds
            f[ 0 ] =  4.0 * r + 4.0 * s + 4.0 * t - 3.0;
            f[ 1 ] =  0.0;
            f[ 2 ] =  4.0 * s - 1.0;
            f[ 3 ] =  0.0;
            f[ 4 ] = -4.0 * r;
            f[ 5 ] = -4.0 * r - 8.0 * s - 4.0 * t + 4.0;
            f[ 6 ] = -4.0 * t;
            f[ 7 ] =  4.0 * r;
            f[ 8 ] =  4.0 * t;
            f[ 9 ] =  0.0;
            return f;
        case 2:// dNdt
            f[ 0 ] =  4.0 * r + 4.0 * s + 4.0 * t - 3.0;
            f[ 1 ] =  0.0;
            f[ 2 ] =  0.0;
            f[ 3 ] =  4.0 * t - 1.0;
            f[ 4 ] = -4.0 * r;
            f[ 5 ] = -4.0 * s;
            f[ 6 ] = -4.0 * r - 4.0 * s - 8.0 * t + 4.0;
            f[ 7 ] =  0.0;
            f[ 8 ] =  4.0 * s;
            f[ 9 ] =  4.0 * r;
            return f;
    }
    throwLengthError( 1, WHERE_AM_I + " undefined values " + toStr( coord ) + " " + toStr( dim ) );
    return f;
}

Hexahedron::Hexahedron( std::vector < Node * > & nodes ) : Cell( nodes ) {
    shape_ = new HexahedronShape();
    fillShape_( );
    neighbourCells_.resize( this->neighbourCellCount(), NULL );
}

Hexahedron::~Hexahedron(){
    delete shape_;
}

// void Hexahedron::findNeighbourCell( uint i ){
//     std::vector < Node * > n( boundaryNodes( i ) );
// 
//     std::set < Cell *> common;
//     std::set < Cell *> commonTmp;
// 
//     if ( n.size() > 1 ) intersectionSet( common, n[ 0 ]->cellSet(), n[ 1 ]->cellSet() );
//     
//     for ( uint j = 2; j < n.size(); j ++ ){
//         commonTmp = common;
//         intersectionSet( common, commonTmp, n[ j ]->cellSet() );
//     }
// 
//     common.erase( this );
//     if ( common.size() == 1 ) neighbourCells_[ i ] = *common.begin(); else neighbourCells_[ i ] = NULL;
// }

void Hexahedron::shapeFunctionsL( const RVector3 & L, RVector & funct ) const{
    funct.resize( nodeCount() );
    THROW_TO_IMPL
}

RVector Hexahedron::deriveNdL( const RVector3 & coord, uint dim ) const {
    RVector f(nodeCount() );
    THROW_TO_IMPL
    return f;
}

std::vector < Node * > Hexahedron::boundaryNodes( uint i ){
    std::vector < Node * > nodes( 4 );
    nodes[ 0 ] = nodeVector_[ HexahedronFacesID[ i ][ 0 ] ];
    nodes[ 1 ] = nodeVector_[ HexahedronFacesID[ i ][ 1 ] ];
    nodes[ 2 ] = nodeVector_[ HexahedronFacesID[ i ][ 2 ] ];
    nodes[ 3 ] = nodeVector_[ HexahedronFacesID[ i ][ 3 ] ];
    return nodes;
}


TriPrism::TriPrism( std::vector < Node * > & nodes ) : Cell( nodes ){
    shape_ = new TriPrismShape();
    fillShape_( );
    neighbourCells_.resize( this->neighbourCellCount(), NULL );
}

TriPrism::~TriPrism(){
    delete shape_;
}

void TriPrism::shapeFunctionsL( const RVector3 & L, RVector & funct ) const {
    funct.resize( nodeCount() );
    THROW_TO_IMPL
}

RVector TriPrism::deriveNdL( const RVector3 & coord, uint dim ) const {
    RVector f(nodeCount() );
    THROW_TO_IMPL
    return f;
}

std::vector < Node * > TriPrism::boundaryNodes( uint i ){
    std::vector < Node * > nodes;
    for ( uint j = 0; j < 3; j ++ ){
        nodes.push_back( nodeVector_[ TriPrismFacesID[ i ][ j ] ] );
    }    
    if ( TriPrismFacesID[ i ][ 3 ] != -1 ) nodes.push_back( nodeVector_[ TriPrismFacesID[ i ][ 3 ] ] );
    
    return nodes;
}

} // namespace GIMLI{
