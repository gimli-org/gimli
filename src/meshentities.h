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

#ifndef _GIMLI_MESHENTITIES__H
#define _GIMLI_MESHENTITIES__H

#include "gimli.h"
#include "vector.h"
#include "matrix.h"
#include "elementmatrix.h"
#include "baseentity.h"

#include <vector>
#include <set>

namespace GIMLI{

template < class Typname > bool lesserId( const Typname * a, const Typname * b){
    return ( a->id() < b->id() );
}

DLLEXPORT Boundary * findBoundary( const Node & n1 );
DLLEXPORT Boundary * findBoundary( const Node & n1, const Node & n2 );
DLLEXPORT Boundary * findBoundary( const Node & n1, const Node & n2, const Node & n3 );

DLLEXPORT Boundary * findBoundary( const std::vector < Node * > & n );

DLLEXPORT Cell * findCommonCell( const std::vector < Node * > & n, bool warn = true );

/*! Functor to collect all nodes of a MeshEntity into a given std::set. */
class CollectNodeFunctor{
public:
    CollectNodeFunctor( std::set< Node * > & c ):c_( & c ){}

    template < class E > void operator( )( const E * e ){
        for ( uint i = 0; i < e->nodeCount( ); i++ ) c_->insert( &e->node( i ) );
    }
    std::set< Node * > * c_;
};

/*! Collect all common \ref Node's for a container of MeshEntities, e.g., commonNodes( std::list< Cell > & cellList ), returns all common nodes that are associated to each cell at the cellList. */
template < class ContainerOfMeshEntities > std::set< Node * > commonNodes( const ContainerOfMeshEntities & c ){
    std::set< Node * > commons;

    CollectNodeFunctor gc( commons );
    for_each( c.begin(), c.end(), gc );
    return commons;
}

class DLLEXPORT MeshEntity : public BaseEntity {
public:
    MeshEntity( ){
        marker_ = 0;
    }

    MeshEntity( std::vector < Node * > & nodes ){
        marker_ = 0;
        setNodes_( nodes );
    }

    virtual ~MeshEntity(){};

    inline Node & node( uint i ) { return *nodeVector_[ i ]; }

    inline Node & node( uint i ) const { return *nodeVector_[ i ]; }

    inline uint nodeCount() const { return nodeVector_.size(); }

    virtual uint dimension() const { return 0; }

    virtual uint rtti() const { return MESH_MESHENTITY_RTTI; }

    virtual uint parentType() const { return MESH_MESHENTITY_RTTI; }

    inline Shape & shape() { return *shape_; }

    inline Shape & shape() const { return *shape_; }

    inline Shape * pShape() { return shape_; }

    RVector3 center() const;

    virtual double attribute() const { return -1.0; }

    double interpolate( const RVector3 & queryPos, const RVector & sol ) const;

    RVector3 grad( const RVector3 & pos, const RVector & data ) const;

    /*! Return shapefunctions N_i < i, N > i=nodeId() for cartesion coordinates at entity pos(x,y,z) */
    std::pair< IndexArray, RVector > shapeFunctions( const RVector3 & pos ) const;

    /*! CHANGE dNdCartesian( pos )
    Return derived shapefunctions for cartesion coordinates dN(L_i)/dxyz at pos
    */
    std::pair< IndexArray, std::vector < RVector3 > > shapeFunctionsDerive( const RVector3 & pos ) const;

    /*!
        Interface for shapefunctions in natural coordinates N(L_i).
        L need to contain L( L_i+1, ... ) cause the relation 1 = \sum L_i*/
    virtual void shapeFunctionsL( const RVector3 & L, RVector & funct ) const{
        std::cout << "mesh entity shapeFunctions " << rtti() << std::endl;
        THROW_TO_IMPL
    }

    /*! CHANGE dNdC( Pos L, coord c, RVec( nN)  )
        Interface for derived shapefunctions dN(L_i)/dL_coord.
        L need to contain L(L_i+1,...) cause the relation 1 = \sum L_i*/
    virtual void shapeFunctionsDeriveL( const RVector3 & L, uint coord, RVector & funct ) const;


    /*! CHANGE to dNdL( No, L ) */
    virtual RVector deriveNdL( const RVector3 & L, uint No ) const{
        std::cout << "mesh entity deriveNdL" << rtti() << std::endl;
        THROW_TO_IMPL
        return RVector(0);
    }

    friend std::ostream & operator << ( std::ostream & str, const MeshEntity & c );

    const std::vector< Node * > & nodes() const { return nodeVector_; }
//   inline std::iterator < std::vector < Node * > > begin() { return nodeVector_.begin(); }
//
//   inline std::vector < Node * >::iterator * end() { return nodeVector_.end(); }

//    void setUxCache( const ElementMatrix < double >  & mat ) const { uxCache_ = mat; }
 
 //   const ElementMatrix < double > & uxCache() const { return uxCache_; }
     
    void setUxCache( const RMatrix & mat ) const { uxCache_ = mat; }
   
    const RMatrix & uxCache() const { return uxCache_; }

protected:
    void fillShape_();

    void setNodes_( std::vector < Node * > & nodes );

    void deRegisterNodes_();

    double interpolate_( const std::pair< IndexArray, RVector > & sF, const RVector & data ) const;

    RVector3 grad_( const std::pair< IndexArray, std::vector < RVector3 > > & sF, const RVector & data ) const;

    Shape * shape_;

    std::vector < Node * > nodeVector_;

    /*! Cache for derivation matrixes */
    //mutable ElementMatrix < double > uxCache_; // to expensive maybe, lightweight baseclass here
    mutable RMatrix uxCache_;
    
private:
    /*! do not copy a mesh entity at all */
    MeshEntity( const MeshEntity & ent ){
        THROW_TO_IMPL
    }
    
    /*! do not assign a mesh entity at all */
    MeshEntity & operator = ( const MeshEntity & ent ){
         if ( this != & ent ){
             THROW_TO_IMPL
         } return *this;
    }
};

class DLLEXPORT Cell : public MeshEntity {
public:
    Cell() : MeshEntity( ), attribute_( 0.0 ), tagged_( false ) { }

    Cell( std::vector < Node * > & nodes ) : MeshEntity( nodes ), attribute_( 0.0 ), tagged_( false ) {
        registerNodes_( );
    }

    ~Cell(){
        deRegisterNodes_();
    }

    inline virtual uint rtti() const { return MESH_CELL_RTTI; }
    inline virtual uint parentType() const { return MESH_CELL_RTTI; }
    inline virtual uint neighbourCellCount() const { return 0; }

    void cleanNeighbourInfos( );

    /*! Return the direct neighbour cell corresponding to local node i. The cell will be searched and stored by the virtual method \ref findNeighbourCell. All neighbouring relationships have to be initialized ones by calling \ref Mesh.createNeighborInfos( ).
     If no cell can be found NULL is returned. */
    inline Cell * neighbourCell( uint i ){ return neighbourCells_[ i ]; }

    /*! Find neighbour cell regarding to the i-th Boundary and store them in neighbourCells_. */
    virtual void findNeighbourCell( uint i );
    
    inline double attribute() const { return attribute_; }

    inline void setAttribute( double attr ) { attribute_ = attr; }

    virtual double jacobianDeterminant() const;

    /*! DEPRECATED????
      Find the node of this cell which is in opposite position to the given boundary. Returns a pointer to the node. The boundary must be part of the cell otherwise, a NULL pointer returns. Works for triangle/edge and tetrahedron/triangleFace*/
    Node * oppositeTo( const Boundary & bound );

//     /*! Find the Boundary of this cell which is in opposite position to the given node. Returns a pointer to the boundary.
//       or a NULL pointer. Works for triangle/edge and tetrahedron/triangleFace*/
//     Boundary * oppositeTo( const Node & node );

    DLLEXPORT friend std::ostream & operator << ( std::ostream & str, const Cell & c );

    /*! Mark the cell. Don't use the tag when you use some cell search. */
    inline void setTagged( bool tagged ){ tagged_ = tagged; }
    /*! Untag the cell */
    inline void untag() { setTagged( false ); }
    /*! Tag the cell */
    inline void tag() { setTagged( true ); }
    /*! Return true if the cell is tagged */
    inline bool tagged() const { return tagged_; }

    /*! Experimental */
    virtual std::vector < Node * > boundaryNodes( uint i ){
        CERR_TO_IMPL
        std::cout << rtti() << std::endl;
        std::vector < Node * > n;
        return n;
    }

protected:
    void registerNodes_( );

    void deRegisterNodes_( );

    
    
    std::vector < Cell * > neighbourCells_;

    double attribute_;

    bool tagged_;

    //std::vector < Cell * > neighbourCells_;

private:
    /*! Don't call this class directly */
    Cell( const Cell & cell ){}
};

class DLLEXPORT Boundary : public MeshEntity{
public:
    Boundary( ) : MeshEntity(){
    }

    Boundary( std::vector < Node * > & nodes )
        : MeshEntity( nodes ), leftCell_( NULL ), rightCell_( NULL ) {
        registerNodes_();
    }

    ~Boundary(){
        deRegisterNodes_();
    }

    inline virtual uint rtti() const { return MESH_BOUNDARY_RTTI; }
    inline virtual uint parentType() const { return MESH_BOUNDARY_RTTI; }

    inline const Cell & leftCell() const { return *leftCell_; }
    inline Cell * leftCell() { return leftCell_; }

    inline const Cell & rightCell() const { return *rightCell_; }
    inline Cell * rightCell() { return rightCell_; }

    inline void setLeftCell(  Cell * cell ) { leftCell_  = cell; }
    inline void setRightCell( Cell * cell ) { rightCell_ = cell; }

    friend std::ostream & operator << ( std::ostream & str, const Boundary & e );

    virtual RVector3 norm( ) const;

protected:
    void registerNodes_();

    void deRegisterNodes_();

    Cell *leftCell_;
    Cell *rightCell_;

private:
    /*! Don't call this class directly */
    Boundary( const Boundary & bound ){
        std::cerr << "Boundary( const Boundary & bound )" << std::endl;
    }
};

class DLLEXPORT NodeBoundary : public Boundary{
public:
    NodeBoundary( Node & n1 );

    NodeBoundary( std::vector < Node * > & nodes );

    virtual ~NodeBoundary();

    inline virtual double jacobianDeterminant() const { return 0.0; }

    inline virtual uint dimension() const { return 1; }

    inline virtual uint rtti() const { return MESH_BOUNDARY_NODE_RTTI; }

    void setNodes( Node & n1, bool changed = true );

    friend std::ostream & operator << ( std::ostream & str, const NodeBoundary & e );

protected:
};

class DLLEXPORT Edge : public Boundary{
public:
    Edge( Node & n1, Node & n2 );

    Edge( std::vector < Node * > & nodes );

    virtual ~Edge();

    inline virtual uint dimension() const { return 2; }

    inline virtual uint rtti() const { return MESH_EDGE_RTTI; }

    void setNodes( Node & n1, Node & n2, bool changed = true  );

    /*! Swap edge between two triangular neighbor cells. Only defined if both neighbors are triangles. */
    int swap();

    void shapeFunctionsL( const RVector3 & L, RVector & funct ) const;

    friend std::ostream & operator << ( std::ostream & str, const Edge & e );

protected:
};

class DLLEXPORT Edge3 : public Edge{
public:

    Edge3( std::vector < Node * > & nodes );

    virtual ~Edge3();

    inline virtual uint rtti() const { return MESH_EDGE3_RTTI; }

    void shapeFunctionsL( const RVector3 & L, RVector & funct ) const;

// friend std::ostream & operator << ( std::ostream & str, const Edge & e );

protected:
};


class DLLEXPORT TriangleFace : public Boundary{
public:
    TriangleFace( Node & n1, Node & n2, Node & n3 );

    TriangleFace( std::vector < Node * > & nodes );

    virtual ~TriangleFace();

    inline virtual uint dimension() const { return 3; }

    inline virtual uint rtti() const { return MESH_TRIANGLEFACE_RTTI; }

    void setNodes( Node & n1, Node & n2, Node & n3, bool changed = true  );

    friend std::ostream & operator << ( std::ostream & str, const TriangleFace & e );

    void shapeFunctionsL( const RVector3 & L, RVector & funct ) const;

//   virtual double interpolate( const RVector3 & queryPos, const RVector & sol ) const;

protected:
};

class DLLEXPORT TriangleFace6 : public TriangleFace{
public:
    TriangleFace6( std::vector < Node * > & nodes );

    ~TriangleFace6();

    virtual uint rtti() const { return MESH_TRIANGLEFACE6_RTTI; }

    void shapeFunctionsL( const RVector3 & L, RVector & funct ) const;

protected:
};

class DLLEXPORT QuadrangleFace : public Boundary{
public:
  QuadrangleFace( Node & n1, Node & n2, Node & n3, Node & n4 );

  QuadrangleFace( std::vector < Node * > & nodes );

  virtual ~QuadrangleFace();

  inline virtual uint dimension() const { return 3; }

  inline virtual uint rtti() const { return MESH_QUADRANGLEFACE_RTTI; }

  void setNodes( Node & n1, Node & n2, Node & n3, Node & n4, bool changed = true  );

  friend std::ostream & operator << ( std::ostream & str, const TriangleFace & e );

//   virtual double interpolate( const RVector3 & queryPos, const RVector & sol ) const;

protected:
private:
    QuadrangleFace( const QuadrangleFace & quad ){
        std::cerr << "QuadrangleFace( const QuadrangleFace & quad )" << std::endl;
    }

};

class DLLEXPORT EdgeCell : public Cell {
public:
    EdgeCell( Node & n1, Node & n2 );

    EdgeCell( std::vector < Node * > & nodes );

    virtual ~EdgeCell();

    inline virtual uint dimension() const { return 1; }

    inline virtual uint rtti() const { return MESH_EDGE_CELL_RTTI; }

    inline virtual uint neighbourCellCount() const { return 2; }

    virtual void findNeighbourCell( uint id );

    void setNodes( Node & n1, Node & n2, bool changed = true  );

    virtual void shapeFunctionsL( const RVector3 & L, RVector & funct ) const;

    /* Return dN d zeta_No( coord ) */
    virtual RVector deriveNdL( const RVector3 & coord, uint No ) const;

    friend std::ostream & operator << ( std::ostream & str, const EdgeCell & t );

protected:
};


/*! count: 0-1, 2(0-1) */
class DLLEXPORT Edge3Cell : public EdgeCell {
public:
    Edge3Cell( std::vector < Node * > & nodes );

    virtual ~Edge3Cell();

    inline virtual uint rtti() const { return MESH_EDGE3_CELL_RTTI; }

    virtual void shapeFunctionsL( const RVector3 & L, RVector & funct ) const;

    virtual RVector deriveNdL( const RVector3 & coord, uint No ) const;

protected:
};

class DLLEXPORT Triangle : public Cell {
public:
    Triangle( Node & n1, Node & n2, Node & n3 );

    Triangle( std::vector < Node * > & nodes );

    virtual ~Triangle();

    inline virtual uint dimension() const { return 2; }

    inline virtual uint rtti() const { return MESH_TRIANGLE_RTTI; }

    inline virtual uint neighbourCellCount() const { return 3; }

    virtual void findNeighbourCell( uint i );

    void setNodes( Node & n1, Node & n2, Node & n3, bool changed = true  );

    virtual void shapeFunctionsL( const RVector3 & L, RVector & funct ) const;

    virtual RVector deriveNdL( const RVector3 & coord, uint No ) const;

    friend std::ostream & operator << ( std::ostream & str, const Triangle & t );

protected:
};

class DLLEXPORT Triangle6 : public Triangle {
public:
    Triangle6( std::vector < Node * > & nodes );

    virtual ~Triangle6();

    inline virtual uint rtti() const { return MESH_TRIANGLE6_RTTI; }

    virtual void shapeFunctionsL( const RVector3 & coord, RVector & funct ) const;

    virtual RVector deriveNdL( const RVector3 & coord, uint No ) const;

protected:
};

//! Quadrangle
/*! Quadrangle

Node direction:
   3-----2
  /     /
 /     /
0-----1

Neighbourship relations:

Neighbour Nr, on Boundary a->b
    0           2->3
    1           3->0
    2           0->1
    3           1->2
*/
class DLLEXPORT Quadrangle : public Cell {
public:
    Quadrangle( Node & n1, Node & n2, Node & n3, Node & n4 );

    Quadrangle( std::vector < Node * > & nodes );

    virtual ~Quadrangle();

    inline virtual uint dimension() const { return 2; }

    virtual uint rtti() const { return MESH_QUADRANGLE_RTTI; }

    void setNodes( Node & n1, Node & n2, Node & n3, Node & n4, bool changed = true  );

    inline virtual uint neighbourCellCount() const { return 4; }

    virtual void findNeighbourCell( uint i );

    virtual void shapeFunctionsL( const RVector3 & L, RVector & funct ) const;

    virtual RVector deriveNdL( const RVector3 & coord, uint dim ) const;

    friend std::ostream & operator << ( std::ostream & str, const Quadrangle & t );

protected:
};

//! Quadrangle8 serendipity class
/*! Quadrangle with 8 nodes for quadratic shapefunctions

Node direction:
   3-----2
  /     /
 /     /
0-----1

count: 0-1-2-3, 4(0-1), 5(1-2), 6(2-3), 7(3-0)
*/
class DLLEXPORT Quadrangle8 : public Quadrangle {
public:
    Quadrangle8( std::vector < Node * > & nodes );

    virtual ~Quadrangle8();

    virtual uint rtti() const { return MESH_QUADRANGLE8_RTTI; }

    virtual void shapeFunctionsL( const RVector3 & L, RVector & funct ) const;

    virtual RVector deriveNdL( const RVector3 & coord, uint dim ) const;

protected:
};

//! A Tetrahedron
/*! A Tet

Node direction:

Neighbourship relations:

Neighbour Nr, on Boundary a-b-c
  intersectionSet( c1, nodeVector_[ ( i + 1 )%4 ]->cellSet(), nodeVector_[ ( i + 2 )%4 ]->cellSet() );
    intersectionSet( common, c1, nodeVector_[ ( i + 3 )%4 ]->cellSet() );
    0           1-2-3
    1           2-3-0
    2           3-2-1
    3           0-1-2
*/
class DLLEXPORT Tetrahedron : public Cell {
public:
    Tetrahedron( Node & n1, Node & n2, Node & n3, Node & n4 );

    Tetrahedron( std::vector < Node * > & nodes );

    virtual ~Tetrahedron();

    inline virtual uint dimension() const { return 3; }

    virtual uint rtti() const { return MESH_TETRAHEDRON_RTTI; }

    void setNodes( Node & n1, Node & n2, Node & n3, Node & n4, bool changed = true  );

    inline virtual uint neighbourCellCount() const { return 4; }

    virtual void findNeighbourCell( uint i );

    virtual void shapeFunctionsL( const RVector3 & L, RVector & funct ) const;

    virtual RVector deriveNdL( const RVector3 & coord, uint dim ) const;

    friend std::ostream & operator << ( std::ostream & str, const Tetrahedron & t );

protected:
};

//*! Zienkiewicz count: 1-2-3-4, 5(1-2), 6(1-3), 7(1-4), 8(2-3), 9(3-4), 10(4-2)* //
//*! VTK,Flaherty,Gimli count: 1-2-3-4, 5(1-2), 6(2-3), 7(3-1), 8(1-4), 9(2-4), 10(3-4)* //
class DLLEXPORT Tetrahedron10 : public Tetrahedron {
public:
    Tetrahedron10( std::vector < Node * > & nodes );

    virtual ~Tetrahedron10();

    virtual uint rtti() const { return MESH_TETRAHEDRON10_RTTI; }

//    void shapeFunctionsDeriveL( const RVector3 & L, uint coord, RVector & funct ) const;
    virtual void shapeFunctionsL( const RVector3 & L, RVector & funct ) const;

    virtual RVector deriveNdL( const RVector3 & coord, uint dim ) const;

protected:
};

//! A Hexahedron
/*! A Hexahedron

Node direction:

  7------6  \n
 /|     /|  \n
4------5 |  \n
| 3----|-2  \n
|/     |/   \n
0------1    \n

Neighbourship relations:

Neighbour Nr, on Boundary a-b-c-d
    0           1-2-6-5
    1           2-3-7-6
    2           0-4-7-3
    3           0-1-5-4
    4           4-5-6-7
    5           0-3-2-1

T.~Apel and N.~Düvelmeyer, Transformation of Hexaedral Finite Element Meshes into Tetrahedral Meshes According to Quality Criteria,
Computing Volume 71, Number 4 / November, 2003, DOI 10.1007/s00607-003-0031-5, Pages   293-304
5-Tet-split: type 6(2) 1-4-5-6, 3-7-6-4, 1-4-0-3, 1-2-3-6, 1-6-4-3
6-Tet-split: type 1    0-1-2-6, 0-2-3-6, 0-1-6-5, 0-4-5-6, 0-3-7-6, 0-4-6-7

*/

static const uint8 HexahedronFacesID[ 6 ][ 4 ] = {
    {1, 2, 6, 5},
    {2, 3, 7, 6},
    {0, 4, 7, 3},
    {0, 1, 5, 4},
    {4, 5, 6, 7},
    {0, 3, 2, 1}
};

static const uint8 HexahedronSplit5TetID[5][4] = {
    {1, 4, 5, 6},
    {3, 6, 7, 4},
    {1, 0, 4, 3},
    {1, 2, 3, 6},
    {1, 4, 6, 3}
};

static const uint8 HexahedronSplit6TetID[6][4] = {
    {0, 1, 2, 6},
    {0, 2, 3, 6},
    {0, 1, 6, 5},
    {0, 4, 5, 6},
    {0, 3, 7, 6},
    {0, 4, 6, 7}
};

class DLLEXPORT Hexahedron: public Cell {
public:
    Hexahedron( std::vector < Node * > & nodes );

    virtual ~Hexahedron();

    inline virtual uint dimension() const { return 3; }

    virtual uint rtti() const { return MESH_HEXAHEDRON_RTTI; }

    inline virtual uint neighbourCellCount() const { return 6; }

    friend std::ostream & operator << ( std::ostream & str, const Hexahedron & t );

    virtual void shapeFunctionsL( const RVector3 & L, RVector & funct ) const;

    virtual RVector deriveNdL( const RVector3 & coord, uint dim ) const;

    /*! Experimental */
    virtual std::vector < Node * > boundaryNodes( uint i );

protected:
};

//! Triangular prism
/*! 
 * A Triangular prism is a three-sided prism. Equivalently, it is a pentahedron of which two faces are parallel.
 * Node direction:

 5    \n
 /|\   \n
3---4  \n
| 2 |  \n
|/ \|  \n
0---1  \n

*/

static const uint8 TriPrismFacesID[ 5 ][ 4 ] = {
    {1, 2, 5, 4},
    {0, 2, 5, 3},
    {0, 1, 4, 3},
    {3, 4, 5, 255},
    {0, 2, 1, 255},
};

class DLLEXPORT TriPrism : public Cell {
public:
    TriPrism( std::vector < Node * > & nodes );

    virtual ~TriPrism();

    inline virtual uint dimension() const { return 3; }

    virtual uint rtti() const { return MESH_TRIPRISM_RTTI; }

    inline virtual uint neighbourCellCount() const { return 5; }

    friend std::ostream & operator << ( std::ostream & str, const Hexahedron & t );

    virtual void shapeFunctionsL( const RVector3 & L, RVector & funct ) const;

    virtual RVector deriveNdL( const RVector3 & coord, uint dim ) const;

    /*! Experimental */
    virtual std::vector < Node * > boundaryNodes( uint i );

protected:
};

} // namespace GIMLI

#endif // MESHENTITIES__H
