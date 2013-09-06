/***************************************************************************
 *   Copyright (C) 2008-2013 by the resistivity.net development team       *
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
    
#include "quaternion.h"
#include "matrix.h"

namespace GIMLI{

std::ostream & operator << ( std::ostream & str, const RQuaternion & q ){
    str << "Quaternion( " << q[ 0 ] << ", " << q[ 1 ] << ", " 
                            << q[ 2 ] << ", " << q[ 3 ] << " ) ";
    return str;
}
                                
RMatrix getRotation( const RVector3 & src, const RVector3 & dest ){
// Based on Stan Melax's article in Game Programming Gems
    RQuaternion q;
    RVector3 v0( src );
    RVector3 v1( dest );
    if ( v0.abs() < TOLERANCE || v1.abs() < TOLERANCE ) {
        q = RQuaternion( 1.0, 0.0, 0.0, 0.0 );
    } else {
        v0.normalise();
        v1.normalise();

        double d = v0.dot( v1 );
    
        if ( ::fabs( ( d - 1.0 ) ) < TOLERANCE ) { //** v1 == v2 
            q = RQuaternion( 1.0, 0.0, 0.0, 0.0 );
        } else if ( ::fabs( ( d + 1.0 ) ) < TOLERANCE ) { //** v1 == -v2
            RVector3 a( RVector3( 1.0, 0.0, 0.0 ).cross( v0 ) );
            if ( a.length() < TOLERANCE ){
                a = RVector3( 0.0, 1.0, 0.0 ).cross( v0 );
            }
            a.normalise();
            q.createFromAxisAngle( a, PI );
        } else {
            double s = std::sqrt( (1.0 + d) * 2.0 );
            RVector3 c = v0.cross( v1 ) / s;
            q = RQuaternion( s * 0.5, c );
            q.normalise();
        }
    }

    RMatrix rot( 3, 3 );
    q.rotMatrix( rot );
    return rot;
}   

                               
} // namespace GIMLI
