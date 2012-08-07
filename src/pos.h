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

#ifndef _GIMLI_POS__H
#define _GIMLI_POS__H

#include "gimli.h"
#include "vector.h"

namespace GIMLI{

template < class ValueType > std::ostream & operator << ( std::ostream & str, const Pos< ValueType > & pos );

DLLEXPORT RVector3 center( const std::vector < RVector3 > & vPos );
DLLEXPORT double jacobianDetXY( const RVector3 & p1, const RVector3 & p2, const RVector3 & p3 );
DLLEXPORT double angle( const RVector3 & p1, const RVector3 & p2, const RVector3 & p3 );

DLLEXPORT bool xVari( const std::vector < RVector3 > & electrodeList );
DLLEXPORT bool yVari( const std::vector < RVector3 > & electrodeList );
DLLEXPORT bool zVari( const std::vector < RVector3 > & electrodeList );

DLLEXPORT RVector x( const std::vector < RVector3 > & rv );
DLLEXPORT RVector y( const std::vector < RVector3 > & rv );
DLLEXPORT RVector z( const std::vector < RVector3 > & rv );

DLLEXPORT void swapXY( std::vector < RVector3 > & rv );
DLLEXPORT void swapXZ( std::vector < RVector3 > & rv );
DLLEXPORT void swapYZ( std::vector < RVector3 > & rv );


DLLEXPORT std::vector < RVector3 > loadRVector3( const std::string & fileName );
DLLEXPORT void saveRVector3( const std::vector < RVector3 > l, const std::string & fileName );

//! 3 dimensional vector
/*! 3 dimensional vector */

template < class ValueType >
class DLLEXPORT Pos {
//template < class ValueType > class DLLEXPORT Pos : public TinyVector< ValueType >{
public:

  //  static const RVector3 ZERO(0.0, 0.0, 0.0);
//   static const RVector3 UNIT_X;
//   static const RVector3 UNIT_Y;
//   static const RVector3 UNIT_Z;
//   static const RVector3 NEGATIVE_UNIT_X;
//   static const RVector3 NEGATIVE_UNIT_Y;
//   static const RVector3 NEGATIVE_UNIT_Z;
//   static const RVector3 UNIT_SCALE;

    /*! Construct an empty in 3d at ( 0, 0, 0 ) */
    Pos( ) : valid_( true ) { assign( 0.0, 0.0, 0.0 ); }

    /*! Construct an empty in 3d at ( 0, 0, 0 ). Optional set valid flag. */
    Pos( bool valid ) : valid_( valid ) { assign( 0.0, 0.0, 0.0 ); }

    //Pos( ValueType x ) : valid_( true ) { assign( x, 0.0, 0.0 ); }
    Pos( ValueType x, ValueType y ) : valid_( true ) { assign( x, y, 0.0 ); }
    Pos( ValueType x, ValueType y, ValueType z ) : valid_( true ) { assign( x, y, z ); }

    Pos( const Pos< ValueType > & pos ) { copy_( pos ); }

    /*! Assignment operator */
    Pos < ValueType > & operator = ( const Pos < ValueType > & pos ){
        if ( this != & pos ){ copy_( pos ); } return *this; }
        
    /*! Assignment operator */
    Pos < ValueType > & operator = ( const Vector < ValueType > & v ){
        if ( v.size() > 2 ) {
            mat_[ 0 ] = v[ 0 ]; 
            mat_[ 1 ] = v[ 1 ]; 
            mat_[ 2 ] = v[ 2 ];
        } else {
            throwLengthError( 1, WHERE_AM_I + " v.size() < 2 " + toStr( v.size() ) );
        }
        return *this; 
        
    }

    inline ValueType & operator [] ( uint i ) { return mat_[ i ]; }

    inline const ValueType & operator [] ( uint i ) const { return mat_[ i ]; }

    #define DEFINE_UNARY_MOD_OPERATOR__( OP ) \
    inline Pos < double > & operator OP##= ( const double & b ){ mat_[ 0 ] OP##= b; mat_[ 1 ] OP##= b; mat_[ 2 ] OP##= b; return *this; }\
    inline Pos < double > & operator OP##= ( const Pos < double > & b ){ mat_[ 0 ] OP##= b[ 0 ]; mat_[ 1 ] OP##= b[ 1 ]; mat_[ 2 ] OP##= b[ 2 ]; return *this; }\

    DEFINE_UNARY_MOD_OPERATOR__( + )
    DEFINE_UNARY_MOD_OPERATOR__( - )
    DEFINE_UNARY_MOD_OPERATOR__( / )
    DEFINE_UNARY_MOD_OPERATOR__( * )

    #undef DEFINE_UNARY_MOD_OPERATOR__

    Pos < ValueType > operator - () const { return *this * -1.0; }

    inline void setValid( bool valid ) { valid_ = valid; }
    inline bool valid() const { return valid_; }

    inline void assign( const ValueType & x, const ValueType & y, const ValueType & z ) {
        mat_[ 0 ] = x; mat_[ 1 ] = y; mat_[ 2 ] = z;
      //  x_ = x; y_ = y; z_ = z;
    }

    inline const ValueType & x() const { return mat_[ 0 ]; }
    inline const ValueType & y() const { return mat_[ 1 ]; }
    inline const ValueType & z() const { return mat_[ 2 ]; }
    inline void setX( ValueType x ) { mat_[ 0 ] = x; }
    inline void setY( ValueType y ) { mat_[ 1 ] = y; }
    inline void setZ( ValueType z ) { mat_[ 2 ] = z; }

    /*! Set a value. Throws out of range exception if index check fails. */
    inline void setVal( const ValueType & val, size_t i ) {
        if ( i >= 0 && i < 3 ) {
            mat_[ i ] = val;
        } else {
            throwRangeError( 1, WHERE_AM_I, i, 0, 3 );
        }
    }

    /*! Get a value. Throws out of range exception if index check fails. */
    inline const ValueType & getVal( size_t i ) const {
        if ( i >= 0 && i < 3 ) {
            return mat_[ i ];
        } else {
            throwRangeError( 1, WHERE_AM_I, i, 0, 3 );
        }
        return mat_[ 0 ];
    }

    Pos< ValueType > & round( double tol ){
        mat_[ 0 ] = rint( mat_[ 0 ] / tol ) * tol;
        mat_[ 1 ] = rint( mat_[ 1 ] / tol ) * tol;
        mat_[ 2 ] = rint( mat_[ 2 ] / tol ) * tol;
        return *this;
    }

    inline double distSquared( const Pos < ValueType > & p ) const {
        return  (  ( mat_[ 0 ] - p[ 0 ] ) * ( mat_[ 0 ] - p[ 0 ] ) +
                   ( mat_[ 1 ] - p[ 1 ] ) * ( mat_[ 1 ] - p[ 1 ] ) +
                   ( mat_[ 2 ] - p[ 2 ] ) * ( mat_[ 2 ] - p[ 2 ] ) );
    }

    inline double dist( const Pos < ValueType > & p ) const { return std::sqrt( distSquared( p ) ); }

    inline double distance( const Pos < ValueType > & p ) const { return dist( p ); }

//     /*! Test if this is faster than dist() */
//     inline double distT( const Pos < ValueType > & p ) const {
//         return std::sqrt(   ( x_ - p.x_ ) * ( x_ - p.x_ ) +
//                             ( y_ - p.y_ ) * ( y_ - p.y_ ) +
//                             ( z_ - p.z_ ) * ( z_ - p.z_ ) );
//     }

    inline double abs( ) const { return length(); }

    inline double distSquared( ) const {
        return  mat_[ 0 ] * mat_[ 0 ] + mat_[ 1 ] * mat_[ 1 ] + mat_[ 2 ] * mat_[ 2 ];
    }

    inline double length( ) const { return std::sqrt( distSquared() ); }

    double angle( const Pos < ValueType > & p ) const {
        double result = acos( this->dot( p ) / ( this->abs() * p.abs() ) );
        if ( std::isnan( result ) || std::isinf( result ) ) {
            result = 0.0;
        }
        return result;
    }

    double angle( const RVector3 & p1, const RVector3 & p3 ) const {
        RVector3 a( p1 - (*this) );
        RVector3 b( p3 - (*this) );
        return ( a ).angle( b );
    }

    inline ValueType dot( const Pos < ValueType > & p ) const {
        return mat_[ 0 ] * p[ 0 ] + mat_[ 1 ] * p[ 1 ] + mat_[ 2 ] * p[ 2 ];
    }

    inline ValueType sum( ) const {
        return mat_[ 0 ] + mat_[ 1 ] + mat_[ 2 ];
    }
    
    Pos< ValueType > norm( const Pos< ValueType > & p1, const Pos< ValueType > & p2 ) const;

    /*! Return normalised copy of this Pos. */
    Pos< ValueType > norm( ) const {
        double t = this->abs();
        if ( t > TOLERANCE ) return *this / t;
        return RVector3( 0.0, 0.0, 0.0 );
    }

    /*! Normalise this Pos and return itself. */
    Pos< ValueType > & normalise( ){
        double t = this->abs();
        if ( t > TOLERANCE ) *this /= t;
        return *this;
    }

    Pos< ValueType > cross( const Pos< ValueType > & p ) const;

    Pos< ValueType > normXY( const Pos< ValueType > & p ) const;

    template < class Matrix > Pos < ValueType > & transform( const Matrix & wm ){
        double x = mat_[ 0 ], y = mat_[ 1 ], z = mat_[ 2 ];

        mat_[ 0 ] = x * wm[ 0 ][ 0 ] + y * wm[ 0 ][ 1 ] + z * wm[ 0 ][ 2 ];
        mat_[ 1 ] = x * wm[ 1 ][ 0 ] + y * wm[ 1 ][ 1 ] + z * wm[ 1 ][ 2 ];
        mat_[ 2 ] = x * wm[ 2 ][ 0 ] + y * wm[ 2 ][ 1 ] + z * wm[ 2 ][ 2 ];
        return *this;
    }

    inline Pos < ValueType > & rotateX( double phi ){
        double mat[ 3 ][ 3 ] ={     { 1.0, 0.0,              0.0},
                                    { 0.0, std::cos( phi ), -std::sin( phi ) },
                                    { 0.0, std::sin( phi ),  std::cos( phi ) } };
        return this->transform( mat );
    }
    inline Pos < ValueType > & rotateY( double phi ){
        double mat[ 3 ][ 3 ] =  {   { std::cos( phi ),  0.0, std::sin( phi ) },
                                    { 0.0,              1.0, 0.0 },
                                    { -std::sin( phi ), 0.0, std::cos( phi ) } };

        return this->transform( mat );
    }
    inline Pos < ValueType > & rotateZ( double phi ){
        double mat[ 3 ][ 3 ] = {    { std::cos( phi ), -std::sin( phi ), 0.0 },
                                    { std::sin( phi ),  std::cos( phi ), 0.0 },
                                    { 0.0,              0.0,             1.0 } };
        return this->transform( mat );
    }

    inline Pos < ValueType > & rotate( const RVector3 & r ){
        return this->rotateX( r[ 0 ] ).rotateY( r[ 1 ] ).rotateZ( r[ 2 ] );
    }
    inline Pos < ValueType > & rotate( double phiX, double phiY, double phiZ ){
        return this->rotateX( phiX ).rotateY( phiY ).rotateZ( phiZ );
    }

    inline Pos < ValueType > & scale( const RVector3 & s ){ return (*this) *= s;}

    inline Pos < ValueType > & translate( const RVector3 & t ){ return (*this) += t;}

    RVector vec() const {
        RVector tmp( 3 );
        tmp[ 0 ] = mat_[ 0 ];
        tmp[ 1 ] = mat_[ 1 ];
        tmp[ 2 ] = mat_[ 2 ];
        return tmp;
    }

//     ValueType x_, y_, z_;

protected:

    inline void copy_( const Pos < ValueType > & pos ) {
        valid_ = pos.valid(); assign( pos[ 0 ], pos[ 1 ], pos[ 2 ] );
    }

    bool valid_;

    ValueType mat_[ 3 ];

};

// template < class ValueType > const RVector3 RVector3::ZERO( 0.0, 0.0, 0.0 );
// template < class ValueType > const RVector3 RVector3::UNIT_X( 1.0, 0.0, 0.0 );
// template < class ValueType > const RVector3 RVector3::UNIT_Y( 0.0, 1.0, 0.0 );
// template < class ValueType > const RVector3 RVector3::UNIT_Z( 0.0, 0.0, 1.0 );
// template < class ValueType > const RVector3 RVector3::NEGATIVE_UNIT_X( -1.0,  0.0,  0.0 );
// template < class ValueType > const RVector3 RVector3::NEGATIVE_UNIT_Y(  0.0, -1.0,  0.0 );
// template < class ValueType > const RVector3 RVector3::NEGATIVE_UNIT_Z(  0.0,  0.0, -1.0 );
// template < class ValueType > const RVector3 RVector3::UNIT_SCALE( 1.0, 1.0, 1.0 );

inline bool operator == ( const RVector3 & a , const RVector3 & b ){
    if ( a.valid() != b.valid() ) return false;
    if ( a.distSquared( b ) < TOLERANCE ) return true; else return false;
}

inline bool operator != ( const RVector3 & a , const RVector3 & b ){
    return !( a == b );
}

inline bool operator < ( const RVector3 & a , const RVector3 & b ){
    std::cout << WHERE_AM_I << std::endl; return true;
}

#define DEFINE_POS_BIN_OPERATOR__( OP )                      \
inline Pos < double > operator OP ( const Pos< double > & a, const Pos< double > & b ){ \
    Pos < double > tmp( a ); return tmp OP##= b; } \
inline Pos < double > operator OP ( const Pos< double > & a, const double & b ){ \
    Pos < double > tmp( a ); return tmp OP##= b; } \
inline Pos < double > operator OP ( const double & a, const Pos< double > & b ){ \
    Pos < double > tmp( a, a, a ); return tmp OP##= b; } \
\

DEFINE_POS_BIN_OPERATOR__( + )
DEFINE_POS_BIN_OPERATOR__( - )
DEFINE_POS_BIN_OPERATOR__( * )
DEFINE_POS_BIN_OPERATOR__( / )

template < class ValueType > std::ostream & operator << ( std::ostream & str, const Pos< ValueType > & pos ){
  if ( pos.valid() ){
    str << pos[ 0 ] << "\t" << pos[ 1 ] << "\t" << pos[ 2 ];
  } else {
    str << " pos is not valid";
  }
  return str;
}

inline bool posLesserX( const RVector3 & a, const RVector3 & b ){
    if ( a[0] == b[0] ) {
        if ( a[1] == b[1] ) {
            return a[2] < b[2];
        } else {
            return a[1] < b[1];
        }
    } else return a[0] < b[0];
}

} // namespace GIMLI;

#endif // _GIMLI_POS__H
