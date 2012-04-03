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

#ifndef _GIMLI_Shape__H
#define _GIMLI_Shape__H

#include "gimli.h"

namespace GIMLI{

class DLLEXPORT Shape {
public:
    Shape();

    virtual ~Shape();

    /*! Pure virtual methode for runtime identification */
    virtual int rtti() const = 0;

    inline virtual std::string name() const { return "Shape"; }

    /*! Get the domain size of this shape, i.e. length, area or volume */
    double domainSize() const;

    /*! Get the jacobian determinant of the shapes domain. */
    double jacobianDeterminant() const;

    /*! Return element coordinates for pos regarding this shape i.e. return (s,t,r) for (x,y,z) */
    virtual RVector3 coordinates( const RVector3 & pos ) const;

    /*! Return real coordinates for locale coordinates. */
    virtual RVector3 xyz( const RVector3 & rst ) const;
    
    /*! Return local coordinates for real coordinates. */
    virtual RVector3 rst( const RVector3 & xyz ) const;
    
    /*! Pure virtual methode for detecting collisions of a pos and the shape */
    virtual bool touch( const RVector3 & pos, bool verbose = false ) const ;

    virtual bool touch1( const RVector3 & pos, bool verbose, int & pFunIdx ) const {
        CERR_TO_IMPL; return 0; }

    /*!DEPRECATED*/
    virtual double partDerivationRealToUnity( uint koord, uint pos ) const {
        CERR_TO_IMPL; return 0;}

    /*! Return \partial naturalCoordinate_i / \partial coord. coord 0 == x; 1 == y; 2 == z. */
    virtual double deriveCoordinates( uint i, uint coord ) const {
        CERR_TO_IMPL; return 0;
    }

    void setNode( uint i, Node & n );

    const Node & node( uint i ) const;

    Node & node( uint i );

    std::vector < Node * > & nodes( ) { return nodeVector_ ; }

    inline void resizeNodeSize( uint n ) { nodeVector_.resize( n, NULL );  }

    uint nodeCount() const { return nodeVector_.size(); }

    const std::vector< Node * > & nodes() const { return nodeVector_; }

    /*! Returns the middle position of this shape */
    RVector3 center() const;

    /*! Returns the normvector if possible otherwise returns nonValid Vector3 */
    virtual RVector3 norm() const;

protected:
    /*! Virtual method to calculate the jacobian determinant of the shapes */
    virtual double jacobianDeterminant_() const { return 0.0; }
    /*! Virtual method to calculate the domain size i.e length, area, volume of the shapes */
    virtual double domainSize_() const { return 0.0; }

    mutable double jacDeterminant_;
    mutable bool hasJacDeterminant_;

    mutable double domSize_;
    mutable bool hasDomSize_;
    
    std::vector < Node * > nodeVector_;
};

DLLEXPORT std::ostream & operator << ( std::ostream & str, const Shape & c );

class DLLEXPORT NodeShape : public Shape{
public:
    NodeShape(){ resizeNodeSize( 1 ); }

    virtual ~NodeShape(){ }

    inline virtual int rtti() const { return MESH_SHAPE_NODE_RTTI; }

    inline virtual std::string name() const { return "NodeShape"; }

    virtual double domainSize_() const { return 0.0; }

    virtual bool touch1( const RVector3 & pos, bool verbose, int & pFunIdx ) const;

    virtual RVector3 norm() const;

protected:
    virtual double jacobianDeterminant_() const { return 0.0; }
};

class DLLEXPORT EdgeShape : public Shape {
public:
    EdgeShape(){ resizeNodeSize( 2 ); }

    virtual ~EdgeShape(){ }

    inline virtual int rtti() const { return MESH_SHAPE_EDGE_RTTI; }

    inline virtual std::string name() const { return "EdgeShape"; }

    virtual RVector3 coordinates( const RVector3 & pos ) const;

    /*! Return n0 + (n1-n0)*rst. */
    virtual RVector3 xyz( const RVector3 & rst ) const;
    
    /*! Returns the partial derivative of the natural coordinates \f$L_i, i=0..1\f$ with respect to the cartesian coordinates
    \ref RVector3 pos \f$ (x) \f$:
        \f$ \frac{ \partial L_i(x) }{ \partial x } \f$ \n
    */
    virtual double deriveCoordinates( uint i, uint coord ) const;

    virtual bool touch( const RVector3 & pos, bool verbose = false ) const ;

    virtual bool touch1( const RVector3 & pos, bool verbose, int & pFunIdx ) const;

    double length() const;

    virtual RVector3 norm() const;

protected:
    /*! Return jacobian determinant of this edge */
    inline virtual double jacobianDeterminant_() const { return length(); }

    inline virtual double domainSize_() const { return length(); }
};

/*! Triangle shape.

Cartesian \f$ (x,y) \f$ to natural coordinates \f$ (L1,L2,L3) \f$ relation:\n
\f{eqnarray*}{
        x-x_1 &=& x_{21}L_2 + x_{31} L_3\\
        y-y_1 &=& y_{21}L_2 + y_{31} L_3\\
        L_1   &=& 1-L_2-L_3
   \f}
whereas \f$ x_{ij} \f$ reads \f$ x_i - x_j \f$ for the \f$x\f$-coordinate of \ref Node \f$i\f$ and \f$j\f$ respectively.
*/
class DLLEXPORT TriangleShape : public Shape {
public:

    TriangleShape(){ resizeNodeSize( 3 ); }

    virtual ~TriangleShape(){ }

    inline virtual int rtti() const { return MESH_SHAPE_TRIANGLE_RTTI; }

    inline virtual std::string name() const { return "TriangleShape"; }

    /*! Returns the unique natural coordinates \ref RVector3 (\f$ L_2(x,y),L_3(x,y),0) \f$
        for a cartesian \ref RVector3 pos \f$ (x,y) \f$ inside this \ref TriangleShape . \n
        m: matrix( [x21,x31], [y21,y31], [z21,z31] )$ \n
        b: matrix( [x-x1],[y-y1] )$ \n
        ls : linsolve_by_lu( m, b )$ \n
        L : ratsimp( first(ls) );
    */
    virtual RVector3 coordinates( const RVector3 & pos ) const;

    virtual RVector3 xyz( const RVector3 & rst ) const;
    
    /*! Returns the partial derivative of the natural coordinates \f$L_i, i=0..2\f$ with respect to the cartesian coordinates
    \ref RVector3 pos \f$ (x,y) \f$:
        \f$ \frac{ \partial L_i(x,y) }{ \partial coord } \f$ \n
        m: matrix( [x21,x31], [y21,y31] )$ \n
        b: matrix( [x-x1],[y-y1] )$ \n
        ls : linsolve_by_lu( m, b )$ \n
        L : ratsimp( first(ls) ); \n
        diff( L, x ); \n
        diff( L, y ); \n
    */
    virtual double deriveCoordinates( uint i, uint coord ) const;

    virtual bool touch1( const RVector3 & pos, bool verbose, int & pFunIdx ) const;

    double area() const;

    virtual RVector3 norm() const;

protected:
    /*! Calculate jacobian determinant
        m: matrix( [x21,x31], [y21,y31], [z21,z31] )$
        ratsimp( determinant( m ) );
    */
    virtual double jacobianDeterminant_() const;

    /*! Interface to get the size i.e. area of this \ref TriangleShape */
    inline virtual double domainSize_() const { return area(); }
};

class DLLEXPORT QuadrangleShape : public Shape {
public:

    QuadrangleShape(){ resizeNodeSize( 4 ); }

    virtual ~QuadrangleShape(){ }

    inline virtual int rtti() const { return MESH_SHAPE_QUADRANGLE_RTTI; }

    inline virtual std::string name() const { return "QuadrangleShape"; }

    double area() const;

    virtual RVector3 coordinates( const RVector3 & pos ) const;

    virtual RVector3 xyz( const RVector3 & rst ) const;
        
    virtual double deriveCoordinates( uint i, uint coord ) const;

    virtual bool touch1( const RVector3 & pos, bool verbose, int & pFunIdx  ) const;

    virtual RVector3 norm() const;

protected:
    virtual double jacobianDeterminant_() const;

    inline virtual double domainSize_() const { return area(); }
};

/*! Tetrahedral shape.

Cartesian \f$ (x,y,z) \f$ to natural coordinates \f$ (L1,L2,L3,L4) \f$ relation:\n
\f{eqnarray*}{
        x-x_1 &=& x_{21} L_2 + x_{31} L_3 + x_{41} L_4\\
        y-y_1 &=& y_{21} L_2 + y_{31} L_3 + y_{41} L_4\\
        z-z_1 &=& z_{21} L_2 + z_{31} L_3 + z_{41} L_4\\
        L_1   &=& 1-L_2-L_3
   \f}
whereas \f$ x_{ij} \f$ reads \f$ x_i - x_j \f$ for the \f$ x \f$-coordinate of \ref Node \f$ i \f$ and \f$ j \f$ respectively.
*/
class DLLEXPORT TetrahedronShape : public Shape {
public:
    TetrahedronShape(){ resizeNodeSize( 4 ); }

    virtual ~TetrahedronShape(){ }

    inline virtual int rtti() const { return MESH_SHAPE_TETRAHEDRON_RTTI; }

    inline virtual std::string name() const { return "TetrahedronShape"; }

    /*! Returns the unique natural coordinates \ref RVector3 (\f$ L_2(x,y,z),L_3(x,y,z),L_4(x,y,z) \f$ for a cartesian
    \ref RVector3 pos \f$ (x,y,z) \f$ inside this \ref TetrahedronShape . \n
        m: matrix( [x21,x31,x41], [y21,y31,y41], [z21,z31,z41] )$ \n
        b: matrix( [x-x1],[y-y1],[z-z1] )$ \n
        ls : linsolve_by_lu( m, b )$ \n
        L : ratsimp( first(ls) );
    */
    virtual RVector3 coordinates( const RVector3 & pos ) const;

    /*! Returns the partial derivative of the natural coordinates \f$ L_coord, i=0..3 \f$ with respect to the cartesian coordinates
    \ref RVector3 pos \f$ (x,y,z) \f$ inside this \ref TetrahedronShape \n
        \f$ \frac{ \partial L_{coord(x,y,z)} }{ \partial dim } \f$ \n
        m: matrix( [x21,x31,x41], [y21,y31,y41], [z21,z31,z41] )$ \n
        b: matrix( [x-x1],[y-y1],[z-z1] )$ \n
        ls : linsolve_by_lu( m, b )$ \n
        L : ratsimp( first(ls) ); \n
        diff( L, x ); \n
        diff( L, y ); \n
        diff( L, z ); \n
    */
    virtual double deriveCoordinates( uint coord, uint dim ) const;

    virtual bool touch1( const RVector3 & pos, bool verbose, int & pFunIdx ) const ;

    /*! Return koord_i - koord_0, koord are 0,1,2 for x,y,z and i = 1,..,3 representing the nodes of hte tetrahedron. No checks are performend. So ensure that the koordinates are valid. */
    double partDerivationRealToUnity( uint koord, uint i ) const ;

    void setNodes( Node * n0, Node * n1, Node * n2, Node * n3 );

    double volume() const;

protected:
    /*! Calculate jacobian determinant \n
        m: matrix( [x21,x31,x41], [y21,y31,y41], [z21,z31,z41] )$ \n
        ratsimp( determinant( m ) );
    */
    virtual double jacobianDeterminant_() const;

    inline virtual double domainSize_() const { return volume(); }
};

class DLLEXPORT HexahedronShape : public Shape {
public:
    HexahedronShape(){ resizeNodeSize( 8 ); }

    virtual ~HexahedronShape(){ }

    inline virtual int rtti() const { return MESH_SHAPE_HEXAHEDRON_RTTI; }

    inline virtual std::string name() const { return "HexahedronShape"; }

    double volume() const;

    virtual RVector3 coordinates( const RVector3 & pos ) const;

    virtual bool touch1( const RVector3 & pos, bool verbose, int & pFunIdx ) const ;

protected:
    virtual double jacobianDeterminant_() const;

    inline virtual double domainSize_() const { return volume(); }
};

} // namespace GIMLI

#endif // _GIMLI_SHAPE__H

