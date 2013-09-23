/***************************************************************************
 *   Copyright (C) 2006-2013 by the resistivity.net development team       *
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

#include "polynomial.h"

#include <vector>
#include <set>

namespace GIMLI{

template < class Typname > bool lesserId(const Typname * a, const Typname * b){
    return (a->id() < b->id());
}

DLLEXPORT Boundary * findBoundary(const Node & n1);
DLLEXPORT Boundary * findBoundary(const Node & n1, const Node & n2);
DLLEXPORT Boundary * findBoundary(const Node & n1, const Node & n2, const Node & n3);

DLLEXPORT Boundary * findBoundary(const std::vector < Node * > & n);

DLLEXPORT Cell * findCommonCell(const std::vector < Node * > & n, bool warn = true);

/*! Functor to collect all nodes of a MeshEntity into a given std::set. */
class CollectNodeFunctor{
public:
    CollectNodeFunctor(std::set< Node * > & c):c_(& c){}

    template < class E > void operator()(E * e){
        for (Index i = 0; i < e->nodeCount(); i++) c_->insert(& e->node(i));
    }
    std::set< Node * > * c_;
};

/*! Collect all common \ref Node's for a container of MeshEntities, e.g.,
 * commonNodes(std::list< Cell > & cellList), returns all common nodes that are
 * associated to each cell at the cellList. */
template < class ContainerOfMeshEntities > 
std::set< Node * > commonNodes(const ContainerOfMeshEntities & c){
    std::set< Node * > commons;

    CollectNodeFunctor gc(commons);
    for_each(c.begin(), c.end(), gc);
    return commons;
}

class DLLEXPORT MeshEntity : public BaseEntity {
public:

    /*! Default constructor.*/
    MeshEntity();
    
    /*! Construct the entity from the given nodes.*/
    MeshEntity(std::vector < Node * > & nodes);
        
    /*! Default destructor.*/
    virtual ~MeshEntity();

    /*! Return the dimension for this MeshEntity. */
    virtual uint dim() const { return 0; }

    /*! Return the runtime identification for this MeshEntity. */
    virtual uint rtti() const { return MESH_MESHENTITY_RTTI; }
    
    /*! To separate between major MeshEntity families e.g. Cell and Boundary. */
    virtual uint parentType() const { return MESH_MESHENTITY_RTTI; }
    
    inline Node & node(uint i) { return *nodeVector_[ i ]; }

    inline Node & node(uint i) const { return *nodeVector_[ i ]; }

    inline uint nodeCount() const { return nodeVector_.size(); }

    inline Shape & shape() { return *shape_; }

    inline Shape & shape() const { return *shape_; }

    inline Shape * pShape() { return shape_; }

    /*! See Shape::rst */
    RVector3 rst(uint i) const;
    
    /*! Return the center coordinates of this MeshEntity. */
    RVector3 center() const;

    virtual double attribute() const { return -1.0; }

    /*! Return IndexArray of all node ids. */
    IndexArray ids() const ;

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;
    
    /*! Return a \ref RVector for the \f$ n=[0,\mathrm{nodeCount()}] \f$ shape functions \f$ N_n(L_1,L_2,L_3)\f$ for the local coordinate \f$ (L_1,L_2,L_3)\f$ */     
    virtual RVector N(const RVector3 & rst) const;
    
    virtual void N(const RVector3 & rst, RVector & n) const;
    
    /*! Return a \ref RVector of the derivation for the \f$ n=[0,\mathrm{nodeCount()}] \f$ shape functions \f$ N_n(L_1,L_2,L_3)\f$ for the local coordinate \f$ (L_1,L_2,L_3)\f$
     *   regarding to the local coordinates \f$ L_i \f$ \n
     * \f$ \frac{\partial N_n(L_1,L_2,L_3)}{\partial L_i} \f$ with may be \f$ i = 0,1,2 \f$*/
    virtual RVector dNdL(const RVector3 & rst, uint i) const;
           
    /*! Interpolate a scalar field at position p for the scalar field u regarding to the shape functions of the entity. 
     * \param p Cartesian coordinates (x,y,z) need to be inside, or on the boundary, of the entity.
     * \param u The field vector u need to be of size mesh.nodeCount() for the corresponding mesh.  */
    double pot(const RVector3 & p, const RVector & u) const;
         
    /*! Return gradient at position pos for field u regarding to the shape functions of the entity. 
    * The field vector u need to be of size mesh.nodeCount() for the corresponding mesh.
    * The position pos, in Cartesian coordinates (x,y,z), need to be inside, or on the boundary, of the entity.  */
    RVector3 grad(const RVector3 & p, const RVector & u) const;

    friend std::ostream & operator << (std::ostream & str, const MeshEntity & c);

    const std::vector< Node * > & nodes() const { return nodeVector_; }
//   inline std::iterator < std::vector < Node * > > begin() { return nodeVector_.begin(); }
//
//   inline std::vector < Node * >::iterator * end() { return nodeVector_.end(); }

//    void setUxCache(const ElementMatrix < double >  & mat) const { uxCache_ = mat; }
 
 //   const ElementMatrix < double > & uxCache() const { return uxCache_; }
     
    void setUxCache(const RMatrix & mat) const { uxCache_ = mat; }
   
    const RMatrix & uxCache() const { return uxCache_; }

protected:
    void fillShape_();

    void setNodes_(std::vector < Node * > & nodes);

    void deRegisterNodes_();

    Shape * shape_;

    std::vector < Node * > nodeVector_;

    /*! Cache for derivation matrixes */
    //mutable ElementMatrix < double > uxCache_; // to expensive maybe, lightweight baseclass here
    mutable RMatrix uxCache_;
    
private:
    /*! do not copy a mesh entity at all */
    MeshEntity(const MeshEntity & ent){
        THROW_TO_IMPL
    }
    
    /*! do not assign a mesh entity at all */
    MeshEntity & operator = (const MeshEntity & ent){
         if (this != & ent){
             THROW_TO_IMPL
         } return *this;
    }
};

//! A abstract cell
/*! Interface class for all cells. */
class DLLEXPORT Cell : public MeshEntity {
public:
    DLLEXPORT friend std::ostream & operator << (std::ostream & str, const Cell & c);
    
    /*! Default constructor. */
    Cell();
    
    /*! Construct cell from vector of nodes. */
    Cell(std::vector < Node * > & nodes);

    /*! Default destructor. */
    ~Cell();

    virtual uint rtti() const { return MESH_CELL_RTTI; }
    virtual uint parentType() const { return MESH_CELL_RTTI; }
    virtual uint neighbourCellCount() const { return 0; }
    inline uint boundaryCount() const { return neighbourCellCount(); }

    void cleanNeighbourInfos();

    Cell * neighbourCell(const RVector & sf);
    
    /*! Return the direct neighbor cell corresponding to local node i.
     * The cell will be searched and stored by the virtual method 
     * \ref findNeighbourCell. 
     * All neighboring relationships have to be initialized ones by calling 
     * \ref Mesh::createNeighborInfos(). 
     * If no cell can be found NULL is returned. */
    inline Cell * neighbourCell(uint i){ return neighbourCells_[ i ]; }

    /*! Find neighbor cell regarding to the i-th Boundary and store them 
     * in neighbourCells_. */
    virtual void findNeighbourCell(uint i);
    
    inline double attribute() const { return attribute_; }

    inline void setAttribute(double attr) { attribute_ = attr; }

    /*! DEPRECATED????
      Find the node of this cell which is in opposite position to the given boundary. Returns a pointer to the node. The boundary must be part of the cell otherwise, a NULL pointer returns. Works for triangle/edge and tetrahedron/triangleFace*/
    Node * oppositeTo(const Boundary & bound);

    /*! Find the nearest boundary to be crossed in direction to the point 
     * responsible for the shape function. */
    Boundary * boundaryTo(const RVector & sf);
    
    /*! Mark the cell. Don't use the tag when you use some cell search. */
    inline void setTagged(bool tagged){ tagged_ = tagged; }
    
    /*! Untag the cell */
    inline void untag() { setTagged(false); }
    
    /*! Tag the cell */
    inline void tag() { setTagged(true); }
    
    /*! Return true if the cell is tagged */
    inline bool tagged() const { return tagged_; }

    /*! Experimental */
    virtual std::vector < Node * > boundaryNodes(Index i){
        CERR_TO_IMPL
        std::cout << rtti() << std::endl;
        std::vector < Node * > n;
        return n;
    }
    
protected:
    void registerNodes_();

    void deRegisterNodes_();
    
    std::vector < Cell * > neighbourCells_;

    double attribute_;

    bool tagged_;

private:
    /*! Don't call this class directly */
    Cell(const Cell & cell){}
};

class DLLEXPORT Boundary : public MeshEntity{
public:
    Boundary() : MeshEntity(){
    }

    Boundary(std::vector < Node * > & nodes)
        : MeshEntity(nodes), leftCell_(NULL), rightCell_(NULL) {
        registerNodes_();
    }

    ~Boundary(){
        deRegisterNodes_();
    }

    virtual uint rtti() const { return MESH_BOUNDARY_RTTI; }
    virtual uint parentType() const { return MESH_BOUNDARY_RTTI; }

    /*! return these coordinates manual until boundary coordinate transformation is done. */
    virtual RVector3 rst(uint i) const;
    
    inline Cell & leftCell() const { return *leftCell_; }
    inline Cell * leftCell() { return leftCell_; }

    inline Cell & rightCell() const { return *rightCell_; }
    inline Cell * rightCell() { return rightCell_; }

    inline void setLeftCell(Cell * cell) { leftCell_  = cell; }
    inline void setRightCell(Cell * cell) { rightCell_ = cell; }

    friend std::ostream & operator << (std::ostream & str, const Boundary & e);

    virtual RVector3 norm() const;

protected:
    void registerNodes_();

    void deRegisterNodes_();

    Cell *leftCell_;
    Cell *rightCell_;

private:
    /*! Don't call this class directly */
    Boundary(const Boundary & bound){
        std::cerr << "Boundary(const Boundary & bound)" << std::endl;
    }
};

class DLLEXPORT NodeBoundary : public Boundary{
public:
    NodeBoundary(Node & n1);

    NodeBoundary(std::vector < Node * > & nodes);

    virtual ~NodeBoundary();

    virtual uint dim() const { return 1; }

    virtual uint rtti() const { return MESH_BOUNDARY_NODE_RTTI; }

    void setNodes(Node & n1, bool changed = true);

    friend std::ostream & operator << (std::ostream & str, const NodeBoundary & e);

protected:
};

class DLLEXPORT Edge : public Boundary{
public:
    Edge(Node & n1, Node & n2);

    Edge(std::vector < Node * > & nodes);

    virtual ~Edge();

    virtual uint dim() const { return 1; }

    virtual uint rtti() const { return MESH_EDGE_RTTI; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    void setNodes(Node & n1, Node & n2, bool changed = true);

    /*! Swap edge between two triangular neighbor cells. Only defined if both neighbors are triangles. */
    int swap();
    
//     /*! See ref MeshEntity::N() */
//     virtual RVector N(const RVector3 & L) const;
    
//     void shapeFunctionsL(const RVector3 & L, RVector & funct) const;

    friend std::ostream & operator << (std::ostream & str, const Edge & e);

protected:
};

class DLLEXPORT Edge3 : public Edge{
public:

    Edge3(std::vector < Node * > & nodes);

    virtual ~Edge3();

    virtual uint rtti() const { return MESH_EDGE3_RTTI; }

    /*! return these coordinates manual until boundary coordinate transformation is done. */
    virtual RVector3 rst(uint i) const;
    
    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;
        
//     /*! See MeshEntity::N() */
//     virtual RVector N(const RVector3 & L) const;
    /*
    void shapeFunctionsL(const RVector3 & L, RVector & funct) const;*/

// friend std::ostream & operator << (std::ostream & str, const Edge & e);

protected:
};


class DLLEXPORT TriangleFace : public Boundary{
public:
    TriangleFace(Node & n1, Node & n2, Node & n3);

    TriangleFace(std::vector < Node * > & nodes);

    virtual ~TriangleFace();

    virtual uint dim() const { return 2; }

    virtual uint rtti() const { return MESH_TRIANGLEFACE_RTTI; }

    void setNodes(Node & n1, Node & n2, Node & n3, bool changed = true);
    
    friend std::ostream & operator << (std::ostream & str, const TriangleFace & e);

    /*! hate this method need refactoring */
    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;
    
protected:
};

class DLLEXPORT Triangle6Face : public TriangleFace{
public:
    Triangle6Face(std::vector < Node * > & nodes);

    ~Triangle6Face();

    virtual uint rtti() const { return MESH_TRIANGLEFACE6_RTTI; }
    
    /*! hate this method need refactoring */
    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;
    
    /*! hate this method need refactoring */
    virtual RVector3 rst(uint i) const;

protected:
};

class DLLEXPORT QuadrangleFace : public Boundary{
public:
    QuadrangleFace(Node & n1, Node & n2, Node & n3, Node & n4);

    QuadrangleFace(std::vector < Node * > & nodes);

    virtual ~QuadrangleFace();

    virtual uint dim() const { return 2; }

    virtual uint rtti() const { return MESH_QUADRANGLEFACE_RTTI; }
    
    /*! hate this method need refactoring */
    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;
    
    void setNodes(Node & n1, Node & n2, Node & n3, Node & n4, bool changed = true);

    friend std::ostream & operator << (std::ostream & str, const TriangleFace & e);

protected:
    
private:
    QuadrangleFace(const QuadrangleFace & quad){
        std::cerr << "QuadrangleFace(const QuadrangleFace & quad)" << std::endl;
    }
};

class DLLEXPORT Quadrangle8Face : public QuadrangleFace{
public:
    
    Quadrangle8Face(std::vector < Node * > & nodes);

    virtual ~Quadrangle8Face();

    virtual uint rtti() const { return MESH_QUADRANGLEFACE8_RTTI; }
    
    /*! hate this method need refactoring */
    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;
    
    /*! hate this method need refactoring */
    virtual RVector3 rst(uint i) const;

protected:
    
private:
    
};

class DLLEXPORT EdgeCell : public Cell {
public:
    EdgeCell(Node & n1, Node & n2);

    EdgeCell(std::vector < Node * > & nodes);

    virtual ~EdgeCell();

    virtual uint dim() const { return 1; }

    virtual uint rtti() const { return MESH_EDGE_CELL_RTTI; }

    virtual uint neighbourCellCount() const { return 2; }

//     virtual void findNeighbourCell(uint id);

    void setNodes(Node & n1, Node & n2, bool changed = true);
    
    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;
     
    virtual std::vector < Node * > boundaryNodes(Index i);
    
    friend std::ostream & operator << (std::ostream & str, const EdgeCell & t);

protected:
};


/*! count: 0-1, 2(0-1) */
class DLLEXPORT Edge3Cell : public EdgeCell {
public:
    Edge3Cell(std::vector < Node * > & nodes);

    virtual ~Edge3Cell();

    virtual uint rtti() const { return MESH_EDGE3_CELL_RTTI; }

    virtual RVector3 rst(uint i) const;
    
    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;
     
protected:
};

class DLLEXPORT Triangle : public Cell {
public:
    Triangle(Node & n1, Node & n2, Node & n3);

    Triangle(std::vector < Node * > & nodes);

    virtual ~Triangle();

    virtual uint dim() const { return 2; }

    virtual uint rtti() const { return MESH_TRIANGLE_RTTI; }

    virtual uint neighbourCellCount() const { return 3; }

//     virtual void findNeighbourCell(uint i);

    void setNodes(Node & n1, Node & n2, Node & n3, bool changed = true);

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;
     
    friend std::ostream & operator << (std::ostream & str, const Triangle & t);

    virtual std::vector < Node * > boundaryNodes(Index i);
    
protected:
};

class DLLEXPORT Triangle6 : public Triangle {
public:
    Triangle6(std::vector < Node * > & nodes);

    virtual ~Triangle6();

    virtual uint rtti() const { return MESH_TRIANGLE6_RTTI; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;
    
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
    Quadrangle(Node & n1, Node & n2, Node & n3, Node & n4);

    Quadrangle(std::vector < Node * > & nodes);

    virtual ~Quadrangle();

    virtual uint dim() const { return 2; }

    virtual uint rtti() const { return MESH_QUADRANGLE_RTTI; }

    void setNodes(Node & n1, Node & n2, Node & n3, Node & n4, bool changed = true);

    virtual uint neighbourCellCount() const { return 4; }

//     virtual void findNeighbourCell(uint i);

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    friend std::ostream & operator << (std::ostream & str, const Quadrangle & t);

    virtual std::vector < Node * > boundaryNodes(Index i);
    
protected:
};

//! Quadrangle8 for serendipity type
/*! Quadrangle serendipity type with 8 nodes for quadratic shape functions
 * 
 * Node direction: \n
 *     3---6---2  \n
 *    /       /   \n
 *   7       5    \n
 *  /       /     \n
 * 0---4---1      \n
 * 
 * count: 0-1-2-3, 4(0-1), 5(1-2), 6(2-3), 7(3-0)
*/
class DLLEXPORT Quadrangle8 : public Quadrangle {
public:
    Quadrangle8(std::vector < Node * > & nodes);

    virtual ~Quadrangle8();

    virtual uint rtti() const { return MESH_QUADRANGLE8_RTTI; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;
    
protected:
};

static const uint8 TetrahedronFacesID[ 4 ][ 3 ] = {
    {1, 2, 3},
    {2, 0, 3},
    {0, 1, 3},
    {0, 2, 1}
};

//! A Tetrahedron
/*! 
 * A Tetrahedron for linear shape functions
*/
class DLLEXPORT Tetrahedron : public Cell {
public:
    Tetrahedron(Node & n1, Node & n2, Node & n3, Node & n4);

    Tetrahedron(std::vector < Node * > & nodes);

    virtual ~Tetrahedron();

    virtual uint dim() const { return 3; }

    virtual uint rtti() const { return MESH_TETRAHEDRON_RTTI; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;
    
    void setNodes(Node & n1, Node & n2, Node & n3, Node & n4, bool changed = true);

    virtual uint neighbourCellCount() const { return 4; }

//     virtual void findNeighbourCell(uint i);

    friend std::ostream & operator << (std::ostream & str, const Tetrahedron & t);

    virtual std::vector < Node * > boundaryNodes(Index i);
    
protected:
};

//*! VTK,Flaherty,Gimli count: 1-2-3-4, 5(1-2), 6(2-3), 7(3-1), 8(1-4), 9(2-4), 10(3-4)* //
static const uint8 Tet10NodeSplit[ 10 ][ 2 ] = {
    {0,0},{1,1},{2,2},{3,3},
    {0,1},{1,2},{2,0},
    {0,3},{1,3},{2,3}
};

//*! Zienkiewicz count: 1-2-3-4, 5(1-2), 6(1-3), 7(1-4), 8(2-3), 9(3-4), 10(4-2)* //
static const uint8 Tet10NodeSplitZienk[ 10 ][ 2 ] = {
    {0,0},{1,1},{2,2},{3,3},
    {0,1},{0,2},{0,3},
    {1,2},{2,3},{3,1}
};

class DLLEXPORT Tetrahedron10 : public Tetrahedron {
public:
    Tetrahedron10(std::vector < Node * > & nodes);

    virtual ~Tetrahedron10();

    virtual uint rtti() const { return MESH_TETRAHEDRON10_RTTI; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;
     
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

class DLLEXPORT Hexahedron: public Cell {
public:
    Hexahedron(std::vector < Node * > & nodes);

    virtual ~Hexahedron();

    virtual uint dim() const { return 3; }

    virtual uint rtti() const { return MESH_HEXAHEDRON_RTTI; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;
    
    virtual uint neighbourCellCount() const { return 6; }

    friend std::ostream & operator << (std::ostream & str, const Hexahedron & t);

    /*! Experimental */
    virtual std::vector < Node * > boundaryNodes(Index i);

protected:
};

static const uint8 Hexahedron20FacesID[ 6 ][ 8 ] = {
    {0,1,5,4,8,17,12,16},
    {1,2,6,5,9,18,13,17},
    {2,3,7,6,10,19,14,18},
    {3,0,4,7,11,16,15,19},
    {0,3,2,1,11,10,9,8},
    {4,5,6,7,12,13,14,15},
};

static const uint8 Hex20NodeSplit[ 20 ][ 2 ] = {
    {0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{6,6},{7,7},
    {0,1},{1,2},{2,3},{3,0},
    {4,5},{5,6},{6,7},{7,4},
    {0,4},{1,5},{2,6},{3,7}
};

//! A Hexahedron with 20 nodes
/*! A Hexahedron with 20 nodes for quadratic base functions of serendipity style

Node direction:

          7-----14------6  \n      
         /|            /|  \n
    t   / |    s      / |  \n
    | 15 19   /     13 18  \n
    | /   |  /      /   |  \n
    |/    | /      /    |  \n
    4-----12------5     |  \n
    |     3-----10|-----2  \n
    |    /        |    /   \n
   16   /         17  /    \n
    | 11          |  9     \n 
    | /           | /      \n
    |/            |/       \n
    0------8-----1-------r \n

*/
class DLLEXPORT Hexahedron20: public Hexahedron {
public:
    Hexahedron20(std::vector < Node * > & nodes);

    virtual ~Hexahedron20();

    virtual uint rtti() const { return MESH_HEXAHEDRON20_RTTI; }
    
    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    /*! Experimental */
    virtual std::vector < Node * > boundaryNodes(Index i);
    
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
    TriPrism(std::vector < Node * > & nodes);

    virtual ~TriPrism();

    virtual uint dim() const { return 3; }

    virtual uint rtti() const { return MESH_TRIPRISM_RTTI; }
    
    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    virtual uint neighbourCellCount() const { return 5; }

    friend std::ostream & operator << (std::ostream & str, const Hexahedron & t);

    /*! Experimental */
    virtual std::vector < Node * > boundaryNodes(Index i);

protected:
};

static const uint8 Prism15NodeSplit[ 15 ][ 2 ] = {
    {0,0},{1,1},{2,2},{3,3},{4,4},{5,5},
    {0,1},{1,2},{2,0},
    {3,4},{4,5},{5,3},
    {0,3},{1,4},{2,5}
};

//! Triangular15 prism
/*! 
 * A Triangular prism with 15 Nodes is a three-sided prism for quadratic base functions.
 
        5      \n
       / \     \n 
     11 14 10  \n
     /  |  \   \n
    3---9---4  \n
    |   2   |  \n
    |  / \  |  \n
   12 8   7 13 \n
    |/     \|  \n
    0---6---1  \n

*/

class DLLEXPORT TriPrism15 : public TriPrism {
public:
    TriPrism15(std::vector < Node * > & nodes);

    virtual ~TriPrism15();

    virtual uint rtti() const { return MESH_TRIPRISM15_RTTI; }
    
    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

protected:
};

//! A Pyramid
/*! A Pyramid cell

Node direction:

          
    t      3------------2  \n
    |     /            /   \n
    |    /            /    \n
    |   /            /     \n
    |  /      5     /      \n 
    | /            /       \n
    |/            /        \n
    0------------1-------r \n

*/

class DLLEXPORT Pyramid : public Cell {
public:
    Pyramid(std::vector < Node * > & nodes);

    virtual ~Pyramid();

    virtual uint dim() const { return 3; }

    virtual uint rtti() const { return MESH_PYRAMID_RTTI; }
    
    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    virtual uint neighbourCellCount() const { return 5; }

    /*! Experimental */
    virtual std::vector < Node * > boundaryNodes(Index i);

protected:
};

//*! VTK,Flaherty,Gimli count: 1-2-3-4, 5(1-2), 6(2-3), 7(3-1), 8(1-4), 9(2-4), 10(3-4)* //
static const uint8 Pyramid13NodeSplit[ 13 ][ 2 ] = {
    {0,0},{1,1},{2,2},{3,3},{4,4},
    {0,1},{1,2},{2,3},{3,0},
    {0,4},{1,4},{2,4},{3,4}
};

//! A Pyramid
/*! A Pyramid cell with 13 nodes for quadratic base functions

Node direction:

          
    t       3-----7------2  \n
    |      /            /    \n
    |     /  12    11  /     \n
    |    /            /      \n
    |   8     4      6       \n
    |  /            /        \n 
    | /  9      10 /         \n
    |/            /          \n
    0------5-----1-------r   \n

*/

class DLLEXPORT Pyramid13 : public Pyramid {
public:
    Pyramid13(std::vector < Node * > & nodes);

    virtual ~Pyramid13();

    virtual uint rtti() const { return MESH_PYRAMID13_RTTI; }
    
    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

protected:
};

} // namespace GIMLI

#endif // MESHENTITIES__H
