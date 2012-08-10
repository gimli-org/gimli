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
    
MeshEntity::MeshEntity( ){
    marker_ = 0;
}

MeshEntity::MeshEntity( std::vector < Node * > & nodes ){
    marker_ = 0;
    setNodes_( nodes );
}

MeshEntity::~MeshEntity(){
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
            for ( Index i = 0; i < shape_->nodeCount(); i ++ ) shape_->setNode( i, node( i ) );
            
            //* create Shape function and cache them to avoid multithreading problems ... 
            ShapeFunctionCache::instance().shapeFunctions( *shape_ );
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

IndexArray MeshEntity::ids() const {
    IndexArray idVec( nodeCount() );
    
    for ( uint i = 0; i < nodeCount(); i ++ ) {
        idVec[ i ] = node( i ).id();
    }
    return idVec;
}

RVector3 MeshEntity::rst( uint i ) const {
    return shape_->rst( node( i ).pos() ).round( TOLERANCE );
}

std::vector < PolynomialFunction < double > > MeshEntity::createShapeFunctions( ) const{
    std::cerr << "need shape function implementation for meshEntity " << rtti() << std::endl;
    THROW_TO_IMPL
    return std::vector < PolynomialFunction < double > > ();
}
 
RVector MeshEntity::N( const RVector3 & rst ) const {
    RVector n( nodeCount(), 0.0 );
    this->N( rst, n );
    return n;
}
    
void MeshEntity::N( const RVector3 & rst, RVector & n ) const {
    const std::vector< PolynomialFunction < double > > &N = ShapeFunctionCache::instance().shapeFunctions( *this );
    
    for ( Index i = 0; i < N.size(); i ++ ) {
        n[ i ] = N[ i ]( rst );
    }
}

RVector MeshEntity::dNdL( const RVector3 & rst, uint i ) const {
        
    const std::vector< PolynomialFunction < double > > &dNL = ShapeFunctionCache::instance().deriveShapeFunctions( *this, i );
    
    RVector ret( dNL.size() );
    
    for ( Index i = 0; i < dNL.size(); i ++ ) {
         ret[ i ] = dNL[ i ]( rst );
    }
    
    return ret;
}
    
    
double MeshEntity::pot( const RVector3 & xyz, const RVector & u ) const {
    return sum( u( this->ids() ) * this->N( shape().rst( xyz ) ) );
}

RVector3 MeshEntity::grad( const RVector3 & xyz, const RVector & u ) const {
  
    RVector3 rst( shape_->rst( xyz ) );
        
    RMatrix MdNdL;
    MdNdL.push_back( dNdL( rst, 0 ) );
    MdNdL.push_back( dNdL( rst, 1 ) );
    MdNdL.push_back( dNdL( rst, 2 ) );
        
    RVector up( u( this->ids() ) );
    RVector3 gr( 3 );
    gr[0] = sum( up * MdNdL.transMult( shape_->invJacobian().col( 0 ) ) );
    gr[1] = sum( up * MdNdL.transMult( shape_->invJacobian().col( 1 ) ) );
    gr[2] = sum( up * MdNdL.transMult( shape_->invJacobian().col( 2 ) ) );
    
    return gr;
    //return grad_( shapeFunctionsDerive( xyz ), u );
}

// RVector3 MeshEntity::grad_( const std::pair< IndexArray, std::vector < RVector3 > > & dN, const RVector & u ) const {
//     
//     RVector3 res(0.0, 0.0, 0.0);
//     for ( uint i = 0; i < 3; i ++ ){
//         //std::cout << i << std::endl;;
//         //for ( std::map< uint, double >::const_iterator it = sF[ i ].begin(); it != sF[ i ].end(); it ++  ){
//           //  std::cout << "   " << it->first << "=" << data[ it->first ] << " " << it->second<< std::endl;;
// 
//         for ( uint j = 0; j < nodeCount(); j ++ ){
//             res[ i ] += u[ dN.first[ j ] ] * dN.second[ j ][ i ];
//         }
//     }
//     //std::cout << res << std::endl;
//     return res;
// }

// void MeshEntity::shapeFunctionsDeriveL( const RVector3 & coord, uint d, RVector & funct ) const {
//     funct.resize( nodeCount() );
//     switch ( this->dim() ){
//         case 1:
//         {
//             double drd = shape_->deriveCoordinates( 0, d );
//             funct = drd * dNdL( coord, 0 );
//         } break;
//         case 2:
//         {
//             double drd = shape_->deriveCoordinates( 0, d );
//             double dsd = shape_->deriveCoordinates( 1, d );
//             funct = drd * dNdL( coord, 0 ) + dsd * dNdL( coord, 1 );
//         } break;
//         case 3:
//         {
//             double drd = shape_->deriveCoordinates( 0, d );
//             double dsd = shape_->deriveCoordinates( 1, d );
//             double dtd = shape_->deriveCoordinates( 2, d );
//             funct = drd * dNdL( coord, 0 ) + dsd * dNdL( coord, 1 ) + dtd * dNdL( coord, 2 );
//         } break;
//     };
// 
// }

// // std::pair< IndexArray, std::vector < RVector3 > > MeshEntity::shapeFunctionsDerive( const RVector3 & xyz ) const {
// //     RVector3 coords = shape().rst( xyz );
// // 
// //     std::vector< RVector3 > dSF( nodeCount(), RVector3( 0.0, 0.0 ) );
// //     RVector funct( nodeCount() );
// // 
// //     for ( uint d = 0; d < dim(); d ++ ){
// //         shapeFunctionsDeriveL( coords, d, funct );
// // 
// //         for ( uint i = 0; i < nodeCount(); i ++ ) {
// //             dSF[ i ][ d ] = funct[ i ];
// //         }
// //     }
// // 
// //     return std::pair< IndexArray, std::vector < RVector3 > >( this->ids(), dSF );
// // }

//############### Cell ##################

Cell::Cell(  ) : MeshEntity( ), attribute_( ), tagged_( false ) {
}
    
Cell::Cell( std::vector < Node * > & nodes ) : MeshEntity( nodes ), attribute_( 0.0 ), tagged_( false ) {
    registerNodes_( );
}

Cell::~Cell(){
    deRegisterNodes_();
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

RVector3 Boundary::rst( uint i ) const {
    if ( this->nodeCount() == shape_->nodeCount() ) return shape_->rst( i );
    std::cerr << "need local coordinate function implementation for meshEntity " << rtti() << std::endl;
    THROW_TO_IMPL
    return RVector3( 0.0, 0.0, 0.0 );
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


std::vector < PolynomialFunction < double > > Edge::createShapeFunctions( ) const {
     return createPolynomialShapeFunctions( *this, 2, true, false );
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

// RVector Edge::N( const RVector3 & L ) const {
//     return shape_->N( L );
// }

Edge3::Edge3( std::vector < Node * > & nodes ) : Edge( nodes ){
}

Edge3::~Edge3(){
}

RVector3 Edge3::rst( uint i ) const {
    if ( i == 2 ) return RVector3( 0.5, 0.0, 0.0 );
    return shape_->rst( i );
}

std::vector < PolynomialFunction < double > > Edge3::createShapeFunctions( ) const {
    return createPolynomialShapeFunctions( *this, 3, true, false );
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

std::vector < PolynomialFunction < double > > TriangleFace::createShapeFunctions( ) const {
    return createPolynomialShapeFunctions( *this, 2, true, false );
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

Triangle6Face::Triangle6Face( std::vector < Node * > & nodes ) : TriangleFace( nodes ){
}

Triangle6Face::~Triangle6Face(){
}

RVector3 Triangle6Face::rst( uint i ) const {
    if ( i == 3 ) return RVector3( ( shape_->rst( 0 ) + shape_->rst( 1 ) ) / 2 );
    if ( i == 4 ) return RVector3( ( shape_->rst( 1 ) + shape_->rst( 2 ) ) / 2 );
    if ( i == 5 ) return RVector3( ( shape_->rst( 2 ) + shape_->rst( 0 ) ) / 2 );
    return shape_->rst( i );
}

std::vector < PolynomialFunction < double > > Triangle6Face::createShapeFunctions( ) const {
    return createPolynomialShapeFunctions( *this, 3, true, false );
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

std::vector < PolynomialFunction < double > > QuadrangleFace::createShapeFunctions( ) const {
    return createPolynomialShapeFunctions( *this, 2, true, true );
}

Quadrangle8Face::Quadrangle8Face( std::vector < Node * > & nodes ) : QuadrangleFace( nodes ){
}

Quadrangle8Face::~Quadrangle8Face(){
}

RVector3 Quadrangle8Face::rst( uint i ) const {
    if ( i == 4 ) return RVector3( ( shape_->rst( 0 ) + shape_->rst( 1 ) ) / 2 );
    if ( i == 5 ) return RVector3( ( shape_->rst( 1 ) + shape_->rst( 2 ) ) / 2 );
    if ( i == 6 ) return RVector3( ( shape_->rst( 2 ) + shape_->rst( 3 ) ) / 2 );
    if ( i == 7 ) return RVector3( ( shape_->rst( 3 ) + shape_->rst( 0 ) ) / 2 );
    return shape_->rst( i );
}

std::vector < PolynomialFunction < double > > Quadrangle8Face::createShapeFunctions( ) const {
    return createPolynomialShapeFunctions( *this, 3, true, true );
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

std::vector < PolynomialFunction < double > > EdgeCell::createShapeFunctions( ) const{
    return createPolynomialShapeFunctions( *this, 2, true, false );
}

// RVector EdgeCell::N( const RVector3 & L ) const{
//     return shape_->N( L );
// }
// 
// RVector EdgeCell::dNdL( const RVector3 & coord, uint dim ) const {
//     RVector dN( nodeCount() );
//     dN[ 0 ] = - 1.0;
//     dN[ 1 ] =   1.0;
//     return dN;
// }

Edge3Cell::Edge3Cell( std::vector < Node * > & nodes ) : EdgeCell( nodes ){
}

Edge3Cell::~Edge3Cell(){
}

std::vector < PolynomialFunction < double > > Edge3Cell::createShapeFunctions( ) const{
    return createPolynomialShapeFunctions( *this, 3, true, false );
}

// RVector Edge3Cell::N( const RVector3 & rst ) const{
//     RVector n( 3 );
//     Edge3_shapeFunctions( rst[ 0 ], n );
//     return n;
// }
// 
// RVector Edge3Cell::dNdL( const RVector3 & coord, uint dim ) const {
//     RVector f( nodeCount() );
//     double r = coord[ 0 ];
//     f[ 0 ] = 4.0 * r - 3.0;
//     f[ 1 ] = 4.0 * r - 1.0;
//     f[ 2 ] = 4.0 - 8.0 * r;
//     return f;
// }

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

std::vector < PolynomialFunction < double > > Triangle::createShapeFunctions( ) const{
    return createPolynomialShapeFunctions( *this, 2, true, false );
}


// RVector Triangle::N( const RVector3 & L ) const {
//     return shape_->N( L );
// }

// // dNdL( Nr, pos ), return for each N
// RVector Triangle::dNdL( const RVector3 & coord, uint dim ) const {
//     RVector f( nodeCount(), 0.0 );
//     switch ( dim ){
//         case 0: // dNdr
//             f[ 0 ] = - 1.0;
//             f[ 1 ] =   1.0;
//             f[ 2 ] =   0.0;
//             return f;
//         case 1:// dNds
//             f[ 0 ] = - 1.0;
//             f[ 1 ] =   0.0;
//             f[ 2 ] =   1.0;
//             return f;
//     }
//     return f;
// }

Triangle6::Triangle6( std::vector < Node * > & nodes ) : Triangle( nodes ){
}

Triangle6::~Triangle6(){
}

std::vector < PolynomialFunction < double > > Triangle6::createShapeFunctions( ) const{
    return createPolynomialShapeFunctions( *this, 3, true, false );
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

std::vector < PolynomialFunction < double > > Quadrangle::createShapeFunctions( ) const{
    return createPolynomialShapeFunctions( *this, 2, true, true );
}

Quadrangle8::Quadrangle8( std::vector < Node * > & nodes ) : Quadrangle( nodes ){
}

Quadrangle8::~Quadrangle8(){
}

std::vector < PolynomialFunction < double > > Quadrangle8::createShapeFunctions( ) const{
//     return createPolynomialShapeFunctions( *this, 3, false, false );
    return createPolynomialShapeFunctions( *this, 3, true, true );
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

std::vector < PolynomialFunction < double > > Tetrahedron::createShapeFunctions( ) const{
    return createPolynomialShapeFunctions( *this, 2, true, false );
}

Tetrahedron10::Tetrahedron10( std::vector < Node * > & nodes ) : Tetrahedron( nodes ){
}

Tetrahedron10::~Tetrahedron10(){
}

std::vector < PolynomialFunction < double > > Tetrahedron10::createShapeFunctions( ) const{
    return createPolynomialShapeFunctions( *this, 3, true, false );
}

Hexahedron::Hexahedron( std::vector < Node * > & nodes ) : Cell( nodes ) {
    shape_ = new HexahedronShape();
    fillShape_( );
    neighbourCells_.resize( this->neighbourCellCount(), NULL );
}

Hexahedron::~Hexahedron(){
    delete shape_;
}

std::vector < PolynomialFunction < double > > Hexahedron::createShapeFunctions( ) const{
//     std::cout << "return createPolynomialShapeFunctions( *this, 2, false, false );" << std::endl;
//     return createPolynomialShapeFunctions( *this, 2, false, false );
    return createPolynomialShapeFunctions( *this, 2, true, true );
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

std::vector < Node * > Hexahedron::boundaryNodes( uint i ){
    std::vector < Node * > nodes( 4 );
    nodes[ 0 ] = nodeVector_[ HexahedronFacesID[ i ][ 0 ] ];
    nodes[ 1 ] = nodeVector_[ HexahedronFacesID[ i ][ 1 ] ];
    nodes[ 2 ] = nodeVector_[ HexahedronFacesID[ i ][ 2 ] ];
    nodes[ 3 ] = nodeVector_[ HexahedronFacesID[ i ][ 3 ] ];
    return nodes;
}

Hexahedron20::Hexahedron20( std::vector < Node * > & nodes ) : Hexahedron( nodes ) {
}

Hexahedron20::~Hexahedron20(){
    
}

std::vector < PolynomialFunction < double > > Hexahedron20::createShapeFunctions( ) const{
//     return createPolynomialShapeFunctions( *this, 3, false, false);
    return createPolynomialShapeFunctions( *this, 3, true, true );
}


TriPrism::TriPrism( std::vector < Node * > & nodes ) : Cell( nodes ){
    shape_ = new TriPrismShape();
    fillShape_( );
    neighbourCells_.resize( this->neighbourCellCount(), NULL );
}

TriPrism::~TriPrism(){
    delete shape_;
}

std::vector < PolynomialFunction < double > > TriPrism::createShapeFunctions( ) const{
    //return createPolynomialShapeFunctions( *this, 2, true, true);
    RVector e2( 2 ); e2[ 0 ] = 0; e2[ 1 ] =  1; // x

    RPolynomialFunction T3_2( e2, RVector( 0 ) );
    RPolynomialFunction T3_3( RVector( 0 ), e2 );
    RPolynomialFunction T3_1 = -( -1.0 + T3_2 + T3_3 );
        
    RPolynomialFunction E2_2T( RVector( 0 ), RVector( 0 ), e2 );
    RPolynomialFunction E2_1T = -( -1.0 + E2_2T );
    std::vector < PolynomialFunction < double > > ret;
    
    ret.push_back( T3_1 * E2_1T );
    ret.push_back( T3_2 * E2_1T );
    ret.push_back( T3_3 * E2_1T );
    ret.push_back( T3_1 * E2_2T );
    ret.push_back( T3_2 * E2_2T );
    ret.push_back( T3_3 * E2_2T );
    return ret;
    
}

std::vector < Node * > TriPrism::boundaryNodes( uint i ){
    std::vector < Node * > nodes;
    for ( uint j = 0; j < 3; j ++ ){
        nodes.push_back( nodeVector_[ TriPrismFacesID[ i ][ j ] ] );
    }    
    if ( TriPrismFacesID[ i ][ 3 ] != 255 ) nodes.push_back( nodeVector_[ TriPrismFacesID[ i ][ 3 ] ] );
    
    return nodes;
}

TriPrism15::TriPrism15( std::vector < Node * > & nodes ) : TriPrism( nodes ){
}

TriPrism15::~TriPrism15(){
}

std::vector < PolynomialFunction < double > > TriPrism15::createShapeFunctions( ) const{
    //return createPolynomialShapeFunctions( *this, 3, false, false);
    
    RVector xy( 9, 0.0 );
    xy[ 0 ] = 1.;  // 1
    xy[ 1 ] = 1.;  // x
    xy[ 2 ] = 1.;  // x^2
    xy[ 3 ] = 1.;  // y
    xy[ 4 ] = 1.;  // yx
    xy[ 5 ] = 0.;  // yx^2
    xy[ 6 ] = 1.;  // y^2
    xy[ 7 ] = 0.;  // y^2x
    xy[ 8 ] = 0.;  // y^2x^2
    
    RVector start( 27 );
    start.setVal( xy, 0, 9 );
    start.setVal( xy, 9, 18 );
    start.setVal( xy, 18, 27 );
    start[ 18 + 2 ] = 0.; // z^2x^2
    start[ 18 + 4 ] = 0.; // z^2xy
    start[ 18 + 6 ] = 0.; // z^2y^2
    
    return createPolynomialShapeFunctions( *this, 3, false, false, start);
    
    
    //#return createPolynomialShapeFunctions( *this, 3, true, true);
}

Pyramid::Pyramid( std::vector < Node * > & nodes ): Cell( nodes ){
    shape_ = new PyramidShape();
    fillShape_( );
    neighbourCells_.resize( this->neighbourCellCount(), NULL );
}

Pyramid::~Pyramid(){
}

std::vector < PolynomialFunction < double > > Pyramid::createShapeFunctions( ) const{
    THROW_TO_IMPL
    return std::vector < PolynomialFunction < double > >();
}
 
std::vector < Node * > Pyramid::boundaryNodes( uint i ){
THROW_TO_IMPL
    return std::vector < Node * > ();
}

Pyramid13::Pyramid13( std::vector < Node * > & nodes ): Pyramid( nodes ){
}

Pyramid13::~Pyramid13(){
}

std::vector < PolynomialFunction < double > > Pyramid13::createShapeFunctions( ) const{
THROW_TO_IMPL
    return std::vector < PolynomialFunction < double > >();
}


    
} // namespace GIMLI{
