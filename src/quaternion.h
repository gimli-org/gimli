/***************************************************************************
 *   Copyright (C) 2008-2011 by the resistivity.net development team       *
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

#ifndef _GIMLI_QUATERNION__H
#define _GIMLI_QUATERNION__H

#include <iostream>
#include "pos.h"

namespace GIMLI{

template < class ValueType > class Quaternion;

typedef Quaternion< double > RQuaternion;

DLLEXPORT RMatrix getRotation( const RVector3 & src, const RVector3 & dest );

template < class ValueType > class DLLEXPORT Quaternion{
public:
    Quaternion( ValueType w = 1.0, ValueType i = 0.0, ValueType j = 0.0, ValueType k = 0.0 )
        : re_( w ){
        im_[ 0 ] = i; im_[ 1 ] = j; im_[ 2 ] = k;
    }

    Quaternion( ValueType re, const Pos< ValueType > & im )
        : re_( re ), im_( im ){ }

    Quaternion( const Pos< ValueType > & xAxis, const Pos< ValueType > & yAxis,
                const Pos< ValueType > & zAxis ){
        THROW_TO_IMPL
    }

    Quaternion( const Quaternion < ValueType > & q ) { copy_( q ); }

    Quaternion & operator = ( const Quaternion < ValueType > & q ) {
        if ( this != & q ){
            copy_( q );
        }
        return *this;
    }

    const ValueType operator [] ( const size_t i ) const {
        if ( i == 0 ) return re_; else return im_[ i - 1 ];
    }

    ValueType & operator [] ( const size_t i ){ if ( i == 0 ) return re_; else return im_[ i - 1 ]; }

    Quaternion & operator /= ( const ValueType & s ) { re_ /= s; im_ /= s; return * this;}
    Quaternion & operator *= ( const ValueType & s ) { re_ *= s; im_ /= s; return * this;}

    void createFromAxisAngle( const Pos< ValueType > & axis, double angle ){
        //** The quaternion representing the rotation is
        //**   q = cos(a/2)+sin(a/2)*(x*i+y*j+z*k)
        double ah = 0.5 * angle;
        re_ = std::cos( ah );
        im_ = axis * std::sin( ah );
    }

    template < class Matrix > void rotMatrix( Matrix & rot ) const {
        ValueType x  = 2.0 * im_[ 0 ], y  = 2.0 * im_[ 1 ], z  = 2.0 * im_[ 2 ];

        ValueType wx = x * re_, wy = y * re_, wz = z * re_;
        ValueType xx = x * im_[ 0 ], xy = y * im_[ 0 ], xz = z * im_[ 0 ];
        ValueType yy = y * im_[ 1 ], yz = z * im_[ 1 ], zz = z * im_[ 2 ];

        rot[ 0 ][ 0 ] = 1.0 - ( yy + zz );
        rot[ 0 ][ 1 ] = xy - wz;
        rot[ 0 ][ 2 ] = xz + wy;
        rot[ 1 ][ 0 ] = xy + wz;
        rot[ 1 ][ 1 ] = 1.0 - ( xx + zz );
        rot[ 1 ][ 2 ] = yz - wx;
        rot[ 2 ][ 0 ] = xz - wy;
        rot[ 2 ][ 1 ] = yz + wx;
        rot[ 2 ][ 2 ] = 1.0 - ( xx + yy );
    }

    inline ValueType norm() const { return re_ * re_ + im_.distSquared( ); }

    inline ValueType length() const { return std::sqrt( norm() ); }

    inline void normalise(){ *this /= length(); }

    inline void setRe( const ValueType re ){ re_ = re; }
    inline ValueType re( ) const { return re_; }

    inline void setIm( const Pos < ValueType > & im ) {  im_ = im; }
    inline Pos < ValueType > im( ) const { return im_; }

protected:

    void copy_( const Quaternion < ValueType > & q ){ re_ = q.re(); im_ = q.im(); }

    ValueType         re_;
    Pos < ValueType > im_;
};

DLLEXPORT std::ostream & operator << ( std::ostream & str, const RQuaternion & q );

} // namespace GIMLI
#endif
