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

#include "shape.h"
#include "node.h"
#include "line.h"
#include "plane.h"

#include "meshentities.h"

namespace GIMLI{

std::ostream & operator << ( std::ostream & str, const Shape & c ){
    str << c.name() << " " << std::endl;
    for ( uint i = 0; i < c.nodeVector_.size(); i ++ ) {
        str << c.nodeVector_[ i ]->pos() << " ";
    }
    return str;
}

Shape::Shape(){
    jacDeterminant_ = 0.0;
    hasJacDeterminant_ = false;
    domSize_ = 0.0;
    hasDomSize_ = false;
}

Shape::~Shape(){
}

Node & Shape::node( uint i ) {
    if ( i < 0 || i > nodeCount() - 1 ){
        std::cerr << WHERE_AM_I << " requested shape node: " << i << " does not exist." << std::endl;
        exit( EXIT_MESH_NO_NODE );
    }
    return *nodeVector_[ i ];
}

const Node & Shape::node( uint i ) const {
  if ( i < 0 || i > nodeCount() - 1 ){
    std::cerr << WHERE_AM_I << " requested shape node: " << i << " does not exist." << std::endl;
    exit( EXIT_MESH_NO_NODE );
  }
  return *nodeVector_[ i ];
}

void Shape::setNode( uint i, Node & n ) {
  if ( i < 0 || i > nodeCount() - 1 ){
    std::cerr << WHERE_AM_I << " requested shape node: " << i << " does not exist." << std::endl;
    exit( EXIT_MESH_NO_NODE );
  }
  nodeVector_[ i ] = &n;
  hasJacDeterminant_ = false;
}

RVector3 Shape::center() const {
  RVector3 center( 0.0, 0.0, 0.0 );
  for ( uint i = 0; i < nodeVector_.size(); i ++ ) {
    center += nodeVector_[ i ]->pos();
  }
  center /= nodeVector_.size();
  return center;
}

RVector3 Shape::norm() const {
    return RVector3();
}

RVector3 Shape::coordinates( const RVector3 & pos ) const {
    std::cout << "shape " << rtti() << std::endl;
    THROW_TO_IMPL
    return RVector3( 0.0, 0.0, 0.0 );
}

bool Shape::touch( const RVector3 & pos, bool verbose ) const {
    int tmp; return touch1( pos, verbose, tmp );
}

double Shape::jacobianDeterminant() const {
    if ( !hasJacDeterminant_ ) {
        jacDeterminant_ = jacobianDeterminant_();
        hasJacDeterminant_ = true;
    }
    return jacDeterminant_;
}

double Shape::domainSize() const {
    if ( !hasDomSize_ ) {
        domSize_ = domainSize_();
        //** gefaehrlich nach nem node::trans, rot, scale stimmt das nichtmehr (siehe unittest:testFEM)
        hasDomSize_ = true;
    }
    return domSize_;
}

bool NodeShape::touch1( const RVector3 & pos, bool verbose, int & tmp ) const {
    tmp = 1;
    if ( pos == nodeVector_[ 0 ]->pos() ) return true; else return false;
}

RVector3 NodeShape::norm() const {
    return RVector3( 0.0, 0.0, 0.0 );
}

//** Start EDGE specific implementation
double EdgeShape::length() const {
  return nodeVector_[ 0 ]->dist( *nodeVector_[ 1 ] );
}

RVector3 EdgeShape::norm() const {
  return nodeVector_[ 0 ]->pos().normXY( nodeVector_[ 1 ]->pos() );
}

RVector3 EdgeShape::coordinates( const RVector3 & pos ) const {
 //** Coordinate transformation
//**     xp = x1 + (x2 - x1) * r
//**     xp1 = x21 * r
//**     r = xp1/x21

    double x21 = nodeVector_[ 1 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];
    double xp1 = pos[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];

    //** find pos in isoparametric koordinates
    //** shape functions N1 = 1-r-s; N2 = r; N3 = s;
    //** if max( N1, N2, N3 ) == 1 pos match the respective node
    RVector3 coords;
    coords[ 0 ] = xp1 / x21;
    return coords;
}

double EdgeShape::deriveCoordinates( uint i, uint coord ) const {
    return 1.0 / ( nodeVector_[ 1 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ] );
//     if ( i == 0 ) return -1.0 / ( nodeVector_[ 1 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ] );
//     else if ( i == 1 ) return 1.0 / ( nodeVector_[ 1 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ] );
    return 0.0;
}

bool EdgeShape::touch( const RVector3 & pos, bool verbose ) const {
  if ( pos == nodeVector_[ 0 ]->pos() || pos == nodeVector_[ 1 ]->pos() ) return true;

  Line line( nodeVector_[ 0 ]->pos(), nodeVector_[ 1 ]->pos() );
  int touch = line.touch( pos );
  if ( touch == 2 || touch == 3 || touch == 4 ) return true;

  //CERR_TO_IMPL; return false ;
  return false;
}

bool EdgeShape::touch1( const RVector3 & pos, bool verbose, int & pFunIdx ) const {
    //** Assuming 1-D mesh here,
    RVector3 coords = coordinates( pos );
    double r = coords[ 0 ];

    double N1 = 1.0 - r;
    double N2 = r;

    //** fun < 0 outside; fun == 0 obBound; fun > 0 inside; max(fun ) = barycenter (1/5 for Edge)
    //** sum of N1,N2,N3 is probably 1.0 if its inside
    double fun = std::min( N1, N2 );
    if ( ( N1 - fun ) < TOUCH_TOLERANCE ) pFunIdx = 0;
    else if ( ( N2 - fun ) < TOUCH_TOLERANCE ) pFunIdx = 1;

    if ( verbose ){
        std::cout << "r: " << r << " Edg: pFunIdx: " << pFunIdx<< std::endl;
        std::cout << " fun: " << fun << ": " << N1 << " " << N2 << std::endl;
    }

    if ( std::fabs( fun ) < TOUCH_TOLERANCE  ) return true; //** on boundary
    if ( fun > 0.0 ) return true; //** inside
    //** outside
    return false;
}

//** Start TRIANGLE specific implementation
double TriangleShape::area() const {
//   RVector3 p1 = nodeVector_[ 0 ]->pos();
//   RVector3 p2 = nodeVector_[ 1 ]->pos();
//   RVector3 p3 = nodeVector_[ 2 ]->pos();

  RVector3 a( nodeVector_[ 1 ]->pos() - nodeVector_[ 0 ]->pos() );
  RVector3 b( nodeVector_[ 2 ]->pos() - nodeVector_[ 0 ]->pos() );

  return ( ( a ).cross( b ) ).abs() * 0.5;

//   return ( nodeVector_[ 1 ]->pos() - nodeVector_[ 0 ]->pos() )
//     .cross( nodeVector_[ 2 ]->pos() - nodeVector_[ 0 ]->pos() ).abs() * 0.5;
}

RVector3 TriangleShape::norm() const{
  RVector3 a( nodeVector_[ 1 ]->pos() - nodeVector_[ 0 ]->pos() );
  RVector3 b( nodeVector_[ 2 ]->pos() - nodeVector_[ 0 ]->pos() );
  RVector3 n( ( a ).cross( b ) );
  return n / n.abs();
}

double TriangleShape::jacobianDeterminant_( ) const {
//   double x1 = nodeVector_[ 0 ]->pos()[ 0 ];
//   double x2 = nodeVector_[ 1 ]->pos()[ 0 ];
//   double x3 = nodeVector_[ 2 ]->pos()[ 0 ];
//   double y1 = nodeVector_[ 0 ]->pos()[ 1 ];
//   double y2 = nodeVector_[ 1 ]->pos()[ 1 ];
//   double y3 = nodeVector_[ 2 ]->pos()[ 1 ];
//   return ( x2 - x1 ) * ( y3 - y1 ) - ( x3 - x1 ) * ( y2 - y1 );

  return area() * 2.0;
  //** this is 2 * area, for testing;
}

RVector3 TriangleShape::coordinates( const RVector3 & pos ) const {
 //** Coordinate transformation
//**     xp = x1 + (x2 - x1) * r + (x3 - x1) * s
//**     yp = y1 + (y2 - y1) * r + (y3 - y1) * s
//**     xp1 = x21 * r + x31 * s
//**     yp1 = y21 * r + y31 * s
//** maxima: eqns: [xp1 = x21 * r + x31 * s, yp1 = y21 * r + y31 * s ]
//** maxima: linsolve(eqns, [r,s])
//**    [r=-(xp1*y31-x31*yp1)/(x31*y21-x21*y31),s=(xp1*y21-x21*yp1)/(x31*y21-x21*y31)]
//**     J = x21 * y31 - x31 * y21

    double x21 = nodeVector_[ 1 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];
    double x31 = nodeVector_[ 2 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];
    double y21 = nodeVector_[ 1 ]->pos()[ 1 ] - nodeVector_[ 0 ]->pos()[ 1 ];
    double y31 = nodeVector_[ 2 ]->pos()[ 1 ] - nodeVector_[ 0 ]->pos()[ 1 ];
    double xp1 = pos[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];
    double yp1 = pos[ 1 ] - nodeVector_[ 0 ]->pos()[ 1 ];

    //** use here the local J instead of jacobianDeterminant_( ), while they use the area-hack
    double J = x21 * y31 - x31 * y21;

    RVector3 coords;
    coords[ 0 ] = ( y31 * xp1 - x31 * yp1 ) / J; // r
    coords[ 1 ] = ( x21 * yp1 - y21 * xp1 ) / J; // s
    return coords;
}

double TriangleShape::deriveCoordinates( uint coord, uint dim ) const{
// return d coord / d dim
    if ( dim == 0 ){
        // drdx == y31
        if ( coord == 0 ) return ( nodeVector_[ 2 ]->pos()[ 1 ] - nodeVector_[ 0 ]->pos()[ 1 ] ) / jacobianDeterminant();
        // dsdx == -y21 == y12
        if ( coord == 1 ) return ( nodeVector_[ 0 ]->pos()[ 1 ] - nodeVector_[ 1 ]->pos()[ 1 ] ) / jacobianDeterminant();

    } else {
        // drdy == -x31 == x13
        if ( coord == 0 ) return ( nodeVector_[ 0 ]->pos()[ 0 ] - nodeVector_[ 2 ]->pos()[ 0 ] ) / jacobianDeterminant();
        // dsdy == x21
        if ( coord == 1 ) return ( nodeVector_[ 1 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ] ) / jacobianDeterminant();
    }
    throwLengthError( 1, WHERE_AM_I + " undefined values " + toStr( coord ) + " " + toStr( dim ) );
    return 0;
}

bool TriangleShape::touch1( const RVector3 & pos, bool verbose, int & pFunIdx ) const {
    RVector3 coords = coordinates( pos );
    double r = coords[ 0 ];
    double s = coords[ 1 ];

    //** find pos in isoparametric koordinates
    //** shape functions N1 = 1-r-s; N2 = r; N3 = s;
    //** if max( N1, N2, N3 ) == 1 pos match the respective node
    double N1 = 1.0 - r - s;
    double N2 = r;
    double N3 = s;

    //** fun < 0 outside; fun == 0 obBound; fun > 0 inside; max(fun ) =barycenter (1/3 for triangle)
    //** sum of N1,N2,N3 is probably 1.0 if its inside
    double fun = std::min( std::min( N1, N2 ), N3 );
    if ( N1 == fun ) pFunIdx = 0;
    else if ( N2 == fun ) pFunIdx = 1;
    else if ( N3 == fun ) pFunIdx = 2;

    if ( verbose ){
        std::cout << "J: jacobianDeterminant() " << std::endl;
        std::cout << "Tri: pFunIdx: " << pFunIdx<< std::endl;
        std::cout << "fun: " << fun << ": " << N1 << " " << N2 << " "<< N3 << std::endl;
    }

    if ( std::fabs( fun ) < TOUCH_TOLERANCE  ) return true; //** on boundary
    if ( fun > 0.0 ) return true; //** inside
    //** outside

    return false;
}

RVector3 QuadrangleShape::coordinates( const RVector3 & pos ) const {
    //** Coordinate transformation
//**     xp = x1 + (x2 - x1 ) * r + (x4 - x1 ) * s - [(x2 - x1 ) + (x4 - x3 )] * r * s
//**     yp = y1 + (y2 - y1 ) * r + (y4 - y1 ) * s - [(y2 - y1 ) + (y4 - y3 )] * r * s
//**     a = xp1
//**     b = x21
//**     c = x41
//**     d = ( x21 + x43 )
//**     e = yp1
//**     f = y21
//**     g = y41
//**     h = ( y21 + y43 )
//**     a = b * r + c * s - d * r * s
//**     e = f * r + g * s - h * r * s
//** maxima: eqns : [a = b * r + c * s - d * r * s, e = f * r + g * s - h * r * s]
//** maxima: solve(eqns, [r,s])
//** mathematica: Solve[eqns, {r,s}]
    double b = nodeVector_[ 1 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];
    double c = nodeVector_[ 3 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];

    double d = b + nodeVector_[ 3 ]->pos()[ 0 ]- nodeVector_[ 2 ]->pos()[ 0 ];
    double f = nodeVector_[ 1 ]->pos()[ 1 ]- nodeVector_[ 0 ]->pos()[ 1 ];
    double g = nodeVector_[ 3 ]->pos()[ 1 ]- nodeVector_[ 0 ]->pos()[ 1 ];
    double h = f + nodeVector_[ 3 ]->pos()[ 1 ]- nodeVector_[ 2 ]->pos()[ 1 ];

    double a = pos[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];
    double e = pos[ 1 ] - nodeVector_[ 0 ]->pos()[ 1 ];

    double r = 0.0;
    double s = 0.0;
    //** assuming the quadrange is a parallelogram
    if ( std::fabs( d ) < TOLERANCE && std::fabs( h ) < TOLERANCE ){
      //** maxima: eqns : [a = b*r + c*s, e = f*r + g*s]
      //** solve(eqns, [r,s])
        r =  ( c*e - a*g ) / ( c*f - b*g );
        s = -( b*e - a*f ) / ( c*f - b*g );
//      std::cout << "0.0 == " << b * r + c * s - d * r * s - a  << std::endl;
//      std::cout << "0.0 == " << f * r + g * s - h * r * s - e  << std::endl;
    } else { // no parallelogramm
        double tmpSqrd = std::pow( ( d*e - c*f + b*g - a*h ), 2.0 ) - 4.0 * ( b*e - a*f )*( d*g - c*h );
        if ( tmpSqrd < TOLERANCE ) {
            std::cout << WHERE_AM_I << " quad: to small " << std::endl;
            return RVector3( 0.0, 0.0, 0.0 );
        }
        double r1 = 1.0 / ( 2.0 * ( d*f - b*h ) ) * ( d*e + c*f - b*g - a*h + std::sqrt( tmpSqrd ) );
        double s1 = 1.0 / ( 2.0 * ( d*g - c*h ) ) * ( d*e - c*f + b*g - a*h - std::sqrt( tmpSqrd ) );
        double r2 = 1.0 / ( 2.0 * ( d*f - b*h ) ) * ( d*e + c*f - b*g - a*h - std::sqrt( tmpSqrd ) );
        double s2 = 1.0 / ( 2.0 * ( d*g - c*h ) ) * ( d*e - c*f + b*g - a*h + std::sqrt( tmpSqrd ) );

        double scale = std::max( std::fabs( r1 ),
                               std::max( std::fabs( s1 ),
                                        std::max( std::fabs( r2 ), std::fabs( s2 ) ) ) );

        double ctr1_1 = b * r1 + c * s1 - d * r1 * s1 - a;
        double ctr1_2 = f * r1 + g * s1 - h * r1 * s1 - e;
        double ctr2_1 = b * r2 + c * s2 - d * r2 * s2 - a;
        double ctr2_2 = f * r2 + g * s2 - h * r2 * s2 - e;
        if ( std::fabs( ctr1_1 ) < TOLERANCE * scale &&  std::fabs( ctr1_2 ) < TOLERANCE * scale ){
        //std::cout << "arb. quad: "  << r1 << " " << s1 << " " << r2 << " " << s2 << std::endl;
            r = r1;
            s = s1;
        } else {
            std::cout << " tmpSqrd " << tmpSqrd << std::endl;
            std::cout << d << " " << h << std::endl;
            std::cout << d*f - b*h << " " << d*g - c*h << std::endl;
            std::cout << r1 << " " << s1 << " " << r2 << " " << s2 << std::endl;

            std::cout << "0 == " << ctr1_1  << std::endl;
            std::cout << "0 == " << ctr1_2  << std::endl;
            std::cout << "0 == " << ctr2_1  << std::endl;
            std::cout << "0 == " << ctr2_2  << std::endl;
            return RVector3( 0.0, 0.0, 0.0 );
        }
    }
    RVector3 coords;
    coords[ 0 ] = r;
    coords[ 1 ] = s;
    return coords;
}

double QuadrangleShape::deriveCoordinates( uint coord, uint dim ) const{
// return d L_coord/d dim
    double b = nodeVector_[ 1 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];
    double c = nodeVector_[ 3 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];

    double d = b + nodeVector_[ 3 ]->pos()[ 0 ]- nodeVector_[ 2 ]->pos()[ 0 ];
    double f = nodeVector_[ 1 ]->pos()[ 1 ]- nodeVector_[ 0 ]->pos()[ 1 ];
    double g = nodeVector_[ 3 ]->pos()[ 1 ]- nodeVector_[ 0 ]->pos()[ 1 ];
    double h = f + nodeVector_[ 3 ]->pos()[ 1 ]- nodeVector_[ 2 ]->pos()[ 1 ];

    //** assuming the quadrangle is a parallelogram
    if ( std::fabs( d ) < TOLERANCE && std::fabs( h ) < TOLERANCE ){

//**         r =  ( c*e - a*g ) / ( c*f - b*g );
//**         s = -( b*e - a*f ) / ( c*f - b*g );
//**         dr/dx; x = a, y = e

        if ( dim == 0 ){
            // dr / dx == -g
            if ( coord == 0 ) return -g / ( c * f - b * g );
            // ds / dx == f
            if ( coord == 1 ) return  f / ( c * f - b * g );

        } else {
            // dr / dy == c
            if ( coord == 0 ) return  c / ( c * f - b * g );
            // ds / dy == -b
            if ( coord == 1 ) return -b / ( c * f - b * g );
        }

    } else {
        THROW_TO_IMPL
        return 0.0;
    }
    throwLengthError( 1, WHERE_AM_I + " undefined values " + toStr( coord ) + " " + toStr( dim ) );
    return 0.0;
}

double QuadrangleShape::area() const {
  //** Gaußsche Trapezformel
  double x13 = nodeVector_[ 0 ]->pos()[ 0 ]- nodeVector_[ 2 ]->pos()[ 0 ];
  double x42 = nodeVector_[ 3 ]->pos()[ 0 ]- nodeVector_[ 1 ]->pos()[ 0 ];

  double y13 = nodeVector_[ 0 ]->pos()[ 1 ]- nodeVector_[ 2 ]->pos()[ 1 ];
  double y24 = nodeVector_[ 1 ]->pos()[ 1 ]- nodeVector_[ 3 ]->pos()[ 1 ];
  return 0.5 * std::fabs( y13 * x42 + y24 * x13 );
    //return std::fabs( this->jacobianDeterminant() );
}

double QuadrangleShape::jacobianDeterminant_() const {
    double x21 = nodeVector_[ 1 ]->pos()[ 0 ]- nodeVector_[ 0 ]->pos()[ 0 ];
    double x32 = nodeVector_[ 2 ]->pos()[ 0 ]- nodeVector_[ 1 ]->pos()[ 0 ];
    double x41 = nodeVector_[ 3 ]->pos()[ 0 ]- nodeVector_[ 0 ]->pos()[ 0 ];
    double x43 = nodeVector_[ 3 ]->pos()[ 0 ]- nodeVector_[ 2 ]->pos()[ 0 ];
    double y21 = nodeVector_[ 1 ]->pos()[ 1 ]- nodeVector_[ 0 ]->pos()[ 1 ];
    double y32 = nodeVector_[ 2 ]->pos()[ 1 ]- nodeVector_[ 1 ]->pos()[ 1 ];
    double y41 = nodeVector_[ 3 ]->pos()[ 1 ]- nodeVector_[ 0 ]->pos()[ 1 ];
    double y43 = nodeVector_[ 3 ]->pos()[ 1 ]- nodeVector_[ 2 ]->pos()[ 1 ];

    return x21 * y41 - y21 * x41 + y21 * x43 - y43 * x21 + y41 * x32 - y32 * x41;
}

RVector3 QuadrangleShape::norm() const {
    //  CERR_TO_IMPL; return RVector3(0.0, 0.0, 0.0);//::ZERO;
    RVector3 a( nodeVector_[ 1 ]->pos() - nodeVector_[ 0 ]->pos() );
    RVector3 b( nodeVector_[ 2 ]->pos() - nodeVector_[ 0 ]->pos() );
    RVector3 n( ( a ).cross( b ) );
    return n / n.abs();
}

bool QuadrangleShape::touch1( const RVector3 & pos, bool verbose, int & pFunIdx  ) const {
    RVector3 coords = coordinates( pos );
    double r = coords[ 0 ];
    double s = coords[ 1 ];

//** shape functions
//     double N1 = ( 1+r )*( 1+s ) * 0.25;
//     double N2 = ( 1-r )*( 1+s ) * 0.25;
//     double N3 = ( 1-r )*( 1-s ) * 0.25;
//     double N4 = ( 1+r )*( 1-s ) * 0.25;

    double N1 = ( 1.0 - r ) * ( 1.0 - s );
    double N2 = r * ( 1.0 - s );
    double N3 = r * s;
    double N4 = s * ( 1.0 - r );

    double fun = std::min( std::min( std::min( N1, N2 ), N3 ), N4 );
    if ( N1 == fun ) pFunIdx = 0;
    else if ( N2 == fun ) pFunIdx = 1;
    else if ( N3 == fun ) pFunIdx = 2;
    else if ( N4 == fun ) pFunIdx = 3;

    if ( std::fabs( fun ) < TOUCH_TOLERANCE  ) return true; //** on boundary
    if ( fun > 0.0 ) return true; //** inside
    //** outside

//     if ( verbose  ) {
//         std::cout << "fun: " << fun << ": " << N1 << " " << N2 << " " << N3 << " " << N4 << std::endl;
//         std::cout << "Quat: pFunIdx " << pFunIdx << std::endl;
//     }

    TriangleShape t;
    //** split quad into 3angle, and test separate
    //** test first
    t.setNode( 0, *nodeVector_[0] );
    t.setNode( 1, *nodeVector_[1] );
    t.setNode( 2, *nodeVector_[3] );
    t.touch1( pos, verbose, pFunIdx );

    //** pFunId == 0, no match -> test 2nd.
    if ( pFunIdx != 0 ){
        pFunIdx += 1;
        if ( verbose  ) {
            std::cout << "  Quat: pFunIdx " << pFunIdx << std::endl;
        }
        return false;
    }
    t.setNode( 0, *nodeVector_[1] );
    t.setNode( 1, *nodeVector_[2] );
    t.setNode( 2, *nodeVector_[3] );
    t.touch1( pos, verbose, pFunIdx );

    //** pFunId == 1, old match -> something wrong
    if ( pFunIdx == 1 ){
        throwLengthError( 1, WHERE_AM_I + "somehting goes wrong here" );
    }
    pFunIdx = (pFunIdx+1)%3;
    if ( verbose  ) {
        std::cout << "  Quat: pFunIdx " << pFunIdx << std::endl;
    }
    return false;
}

//** Start TETRAHEDRON specific implementation
double TetrahedronShape::volume() const {
    RVector3 a( nodeVector_[ 1 ]->pos() - nodeVector_[ 0 ]->pos() );
    RVector3 b( nodeVector_[ 2 ]->pos() - nodeVector_[ 0 ]->pos() );
    RVector3 c( nodeVector_[ 3 ]->pos() - nodeVector_[ 0 ]->pos() );
//     std::cout << a << " " << b << " " << c << " " << a.cross( b ) << std::endl;
//     std::cout << std::fabs( ( a.cross( b ).dot( c ) ) ) << std::endl;
    return 1.0 / 6.0 * std::fabs( ( a.cross( b ).dot( c ) ) );
    //** pls check whish way is faster, profile, with valgrind and count ticks
    return fabs( jacobianDeterminant() / 6.0 );
}

double TetrahedronShape::partDerivationRealToUnity( uint koord, uint i ) const {
    return nodeVector_[ i ]->pos()[ koord ] - nodeVector_[ 0 ]->pos()[ koord ];
}

// double TetrahedronShape::partDerivationRealToUnityOld( char absPos, int relPos ) const {
//
//     double result = 0.0;
//   switch( absPos ){
//   case 'x':
//     result = node( relPos ).x() - node( 0 ).x();
//     break;
//   case 'y':
//     result = node( relPos ).y() - node( 0 ).y();
//     break;
//   case 'z':
//     result = node( relPos ).z() - node( 0 ).z();
//     break;
//   default : std::cerr << WHERE_AM_I << " Dimension invalid. " << std::endl;
//   }
//   return result;
// }

double TetrahedronShape::jacobianDeterminant_() const {
    double x21 = nodeVector_[ 1 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];
    double x31 = nodeVector_[ 2 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];
    double x41 = nodeVector_[ 3 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];
    double y21 = nodeVector_[ 1 ]->pos()[ 1 ] - nodeVector_[ 0 ]->pos()[ 1 ];
    double y31 = nodeVector_[ 2 ]->pos()[ 1 ] - nodeVector_[ 0 ]->pos()[ 1 ];
    double y41 = nodeVector_[ 3 ]->pos()[ 1 ] - nodeVector_[ 0 ]->pos()[ 1 ];
    double z21 = nodeVector_[ 1 ]->pos()[ 2 ] - nodeVector_[ 0 ]->pos()[ 2 ];
    double z31 = nodeVector_[ 2 ]->pos()[ 2 ] - nodeVector_[ 0 ]->pos()[ 2 ];
    double z41 = nodeVector_[ 3 ]->pos()[ 2 ] - nodeVector_[ 0 ]->pos()[ 2 ];

    return x21 * ( y31 * z41 - y41 * z31 ) +
           x31 * ( y41 * z21 - y21 * z41 ) +
           x41 * ( y21 * z31 - y31 * z21 );
}

void TetrahedronShape::setNodes( Node * n0, Node * n1, Node * n2, Node * n3 ){
  setNode( 0, *n0 ); setNode( 1, *n1 ); setNode( 2, *n2 ); setNode( 3, *n3 );
}

RVector3 TetrahedronShape::coordinates( const RVector3 & pos ) const {
      //** Coordinate transformation
//**     xp = x1 + (x2 - x1) * r + (x3 - x1) * s + (x4 - x1) * t
//**     yp = y1 + (y2 - y1) * r + (y3 - y1) * s + (y4 - y1) * t
//**     zp = z1 + (z2 - z1) * r + (z3 - z1) * s + (z4 - z1) * t
//**     xp1 = x21 * r + x31 * s + x41 * t
//**     yp1 = y21 * r + y31 * s + y41 * t
//**     zp1 = z21 * r + z31 * s + z41 * t
//** maxima: eqns: [xp1 = x21 * r + x31 * s + x41 * t, yp1 = y21 * r + y31 * s + y41 * t, zp1 = z21 * r + z31 * s + z41 * t]
//** maxima: linsolve(eqns, [r,s,t])
//     r=(x31*(yp1*z41-y41*zp1)+x41*(y31*zp1-yp1*z31)+xp1*(y41*z31-y31*z41))/
//             (x21*(y41*z31-y31*z41)+x31*(y21*z41-y41*z21)+x41*(y31*z21-y21*z31)),
//     s=-(x21*(yp1*z41-y41*zp1)+x41*(y21*zp1-yp1*z21)+xp1*(y41*z21-y21*z41))/
//             (x21*(y41*z31-y31*z41)+x31*(y21*z41-y41*z21)+x41*(y31*z21-y21*z31)),
//     t=(x21*(yp1*z31-y31*zp1)+x31*(y21*zp1-yp1*z21)+xp1*(y31*z21-y21*z31))/
//             (x21*(y41*z31-y31*z41)+x31*(y21*z41-y41*z21)+x41*(y31*z21-y21*z31))]
//               Jac = -(x21*(y41*z31-y31*z41)+x31*(y21*z41-y41*z21)+x41*(y31*z21-y21*z31))
    double x21 = nodeVector_[ 1 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];
    double x31 = nodeVector_[ 2 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];
    double x41 = nodeVector_[ 3 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];
    double y21 = nodeVector_[ 1 ]->pos()[ 1 ] - nodeVector_[ 0 ]->pos()[ 1 ];
    double y31 = nodeVector_[ 2 ]->pos()[ 1 ] - nodeVector_[ 0 ]->pos()[ 1 ];
    double y41 = nodeVector_[ 3 ]->pos()[ 1 ] - nodeVector_[ 0 ]->pos()[ 1 ];
    double z21 = nodeVector_[ 1 ]->pos()[ 2 ] - nodeVector_[ 0 ]->pos()[ 2 ];
    double z31 = nodeVector_[ 2 ]->pos()[ 2 ] - nodeVector_[ 0 ]->pos()[ 2 ];
    double z41 = nodeVector_[ 3 ]->pos()[ 2 ] - nodeVector_[ 0 ]->pos()[ 2 ];
    double xp1 = pos[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];
    double yp1 = pos[ 1 ] - nodeVector_[ 0 ]->pos()[ 1 ];
    double zp1 = pos[ 2 ] - nodeVector_[ 0 ]->pos()[ 2 ];

    //** find pos in isoparametric koordinates
    //** shape functions N1 = 1-r-s-t; N2 = r; N3 = s; N4 = t
    //** if max( N1, N2, N3 ) == 1 pos match the respective node
//     det( n1, n2, n3, n4 ) =  x21 * ( y31 * z41 - y41 * z31 ) +
//                              x31 * ( y41 * z21 - y21 * z41 ) +
//                              x41 * ( y21 * z31 - y31 * z21 );
//     N2 = det( n1, p, n3, n4 ) / det; // replace every 2 by p
//     N3 = det( n1, n2, p, n4 ) / det; // replace every 3 by p
//     N4 = det( n1, n2, n3, p ) / det; // replace every 4 by p

    RVector3 coords;
    coords[ 0 ] = ( xp1 * ( y31 * z41 - y41 * z31 ) +
                    x31 * ( y41 * zp1 - yp1 * z41 ) +
                    x41 * ( yp1 * z31 - y31 * zp1 ) ) / jacobianDeterminant();

    coords[ 1 ] = ( x21 * ( yp1 * z41 - y41 * zp1 ) +
                    xp1 * ( y41 * z21 - y21 * z41 ) +
                    x41 * ( y21 * zp1 - yp1 * z21 ) ) / jacobianDeterminant();

    coords[ 2 ] = ( x21 * ( y31 * zp1 - yp1 * z31 ) +
                    x31 * ( yp1 * z21 - y21 * zp1 ) +
                    xp1 * ( y21 * z31 - y31 * z21 ) ) / jacobianDeterminant();

    return coords;
}

double TetrahedronShape::deriveCoordinates( uint coord, uint dim) const{
    double x21 = nodeVector_[ 1 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];
    double x31 = nodeVector_[ 2 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];
    double x41 = nodeVector_[ 3 ]->pos()[ 0 ] - nodeVector_[ 0 ]->pos()[ 0 ];
    double y21 = nodeVector_[ 1 ]->pos()[ 1 ] - nodeVector_[ 0 ]->pos()[ 1 ];
    double y31 = nodeVector_[ 2 ]->pos()[ 1 ] - nodeVector_[ 0 ]->pos()[ 1 ];
    double y41 = nodeVector_[ 3 ]->pos()[ 1 ] - nodeVector_[ 0 ]->pos()[ 1 ];
    double z21 = nodeVector_[ 1 ]->pos()[ 2 ] - nodeVector_[ 0 ]->pos()[ 2 ];
    double z31 = nodeVector_[ 2 ]->pos()[ 2 ] - nodeVector_[ 0 ]->pos()[ 2 ];
    double z41 = nodeVector_[ 3 ]->pos()[ 2 ] - nodeVector_[ 0 ]->pos()[ 2 ];

    double ret = 0;

    if ( dim == 0 ){ // dx
//     (y31 z41 - y41 z31) /((x21 y31 - x31 y21) z41 + (x41 y21 - x21 y41) z31 + (x31 y41 - x41 y31) z21)
//     (y41 z21 - y21 z41) /((x21 y31 - x31 y21) z41 + (x41 y21 - x21 y41) z31 + (x31 y41 - x41 y31) z21)
//     (y21 z31 - y31 z21) /((x21 y31 - x31 y21) z41 + (x41 y21 - x21 y41) z31 + (x31 y41 - x41 y31) z21)
        switch ( coord ){
        case 0: // dr/dx
            ret = y31 * z41 - y41 * z31; break;
        case 1: // ds/dx
            ret = y41 * z21 - y21 * z41; break;
        case 2: // dt/dx
            ret = y21 * z31 - y31 * z21; break;
        }
    } else if ( dim == 1 ){ // dy
//         (x41 z31 - x31 z41) / ((x21 y31 - x31 y21) z41 + (x41 y21 - x21 y41) z31 + (x31 y41 - x41 y31) z21)
//         (x21 z41 - x41 z21) / ((x21 y31 - x31 y21) z41 + (x41 y21 - x21 y41) z31 + (x31 y41 - x41 y31) z21)
//         (x31 z21 - x21 z31) / ((x21 y31 - x31 y21) z41 + (x41 y21 - x21 y41) z31 + (x31 y41 - x41 y31) z21)
        switch ( coord ){
        case 0: // dr/dy
            ret = x41 * z31 - x31 * z41; break;
        case 1: // ds/dy
            ret = x21 * z41 - x41 * z21; break;
        case 2: // dt/dy
            ret = x31 * z21 - x21 * z31; break;
        }

    } else { // dz
//         (x31 y41 - x41 y31) / ((x21 y31 - x31 y21) z41 + (x41 y21 - x21 y41) z31 + (x31 y41 - x41 y31) z21)
//         (x41 y21 - x21 y41) / ((x21 y31 - x31 y21) z41 + (x41 y21 - x21 y41) z31 + (x31 y41 - x41 y31) z21)
//         (x21 y31 - x31 y21) / ((x21 y31 - x31 y21) z41 + (x41 y21 - x21 y41) z31 + (x31 y41 - x41 y31) z21)
    switch ( coord ){
        case 0: // dr/dz
            ret = x31 * y41 - x41 * y31; break;
        case 1: // ds/dz
            ret = x41 * y21 - x21 * y41; break;
        case 2: // dt/dz
            ret = x21 * y31 - x31 * y21; break;
        }
    }
   // std::cout << i << " " << coord << " " << ret << std::endl;
    return ret / jacobianDeterminant();
}

bool TetrahedronShape::touch1( const RVector3 & pos, bool verbose, int & pFunIdx ) const {
    RVector3 coords = coordinates( pos );
    double r = coords[ 0 ];
    double s = coords[ 1 ];
    double t = coords[ 2 ];

    double N1 = 1.0 - r - s - t;
    double N2 = r;
    double N3 = s;
    double N4 = t;

    //** fun < 0 outside; fun == 0 bound; fun > 0 inside; max(fun ) =barycenter (1/4 for tetrahedron)
    double fun = std::min( std::min( std::min( N1, N2 ), N3 ), N4 );
    if ( N1 == fun ) pFunIdx = 0;
    else if ( N2 == fun ) pFunIdx = 1;
    else if ( N3 == fun ) pFunIdx = 2;
    else if ( N4 == fun ) pFunIdx = 3;

    if ( verbose ){
        std::cout << "Jac: " << jacobianDeterminant() << " ";
        std::cout << "Tet: pFunIdx: " << pFunIdx<< std::endl;
        std::cout << " fun: " << fun << ": " << N1 << " " << N2 << " "<< N3 << " " << N4 << std::endl;
    }

    if ( std::fabs( fun ) < max( TOUCH_TOLERANCE, TOUCH_TOLERANCE * pos.abs() ) ) return true; //** on boundary
    if ( fun > 0.0 ) return true; //** inside
    //** outside
    return false;
}

double HexahedronShape::volume() const {
    double sum = 0.0;
    TetrahedronShape tet;
    for ( uint i = 0; i < 5; i ++ ){
        tet.setNodes( nodeVector_[ HexahedronSplit5TetID[ i ][ 0 ] ], nodeVector_[ HexahedronSplit5TetID[ i ][ 1 ] ],
                      nodeVector_[ HexahedronSplit5TetID[ i ][ 2 ] ], nodeVector_[ HexahedronSplit5TetID[ i ][ 3 ] ] );
        sum += tet.volume();
    }

    return sum;
}

RVector3 HexahedronShape::coordinates( const RVector3 & pos ) const {
    THROW_TO_IMPL
    return RVector3();
}

bool HexahedronShape::touch1( const RVector3 & pos, bool verbose, int & pFunIdx ) const {
    THROW_TO_IMPL
    return false;
}

double HexahedronShape::jacobianDeterminant_() const {
    THROW_TO_IMPL
    return -1;
}


} // namespace GIMLI

