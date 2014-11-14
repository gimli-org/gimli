/***************************************************************************
 *   Copyright (C) 2008-2014 by the resistivity.net development team       *
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

#ifndef GIMLI_LINE__H
#define GIMLI_LINE__H

#include "gimli.h"
#include "pos.h"

namespace GIMLI{

//! A line
/*! A straight line in the geometrically point of view. */

class DLLEXPORT Line {
public:
    /*! Default constructor. Constructs an invalid straight line. */
    Line( );

    /*! Constructor, construct a line between p0 and the orign [0.0, 0.0, 0.0 ]. */
    Line( const RVector3 & p0 );

    /*! Constructor construct a line betweem p0, and p1.*/
    Line( const RVector3 & p0, const RVector3 & p1 );

    /*! Copycontructor */
    Line( const Line & line );

    /*! Default destructor. Destroys the object.*/
    ~Line();

    /*! Assignment operator*/
    Line & operator = ( const Line & line );

    /*! Functor to the line function. ret = p0_ + ( p1_ - p0_ ) * t; */
    inline const RVector3 operator()( double t ) const { return lineAt( t ); }

    /*! Equal_to operator */
    inline bool operator == ( const Line & line ) const { return this->compare( line ); }

    /*! Not_equal_to operator */
    inline bool operator != ( const Line & line ) const { return !(*this == line); }

    /*! Check validility of the straight line. False if [ p0-p1] < tol. */
    bool checkValidity( double tol = TOLERANCE );

    /*! Return the validity of this line.*/
    inline bool valid() const { return valid_; }

    /*! Compare two lines with a given tolerance. Check if both ref points from line touch this line..
    distance( this, line.p1 )< tol && distance( this, line.p2 )< tol */
    bool compare( const Line & line, double tol = TOLERANCE ) const;

    /*! Return const reference to ref. point p0 */
    const RVector3 & p0( ) const { return p0_; }

    /*! Return const reference to ref. point p1 */
    const RVector3 & p1( ) const { return p1_; }

    /*! Returns the intersection point between line and this line.
        Returns an invalid RVector3 of both do not touch, are equal or paralell. */
    RVector3 intersect( const Line & line) const;

    /*! Returns line parameter t that represents the nearest position to the line for the given point p. Feed it to lineAt and you get the nearest point. */
    double nearest( const RVector3 & p ) const;

    /*! Apply the line as function. ret = p0_ + ( p1_ - p0_ ) * t;*/
    inline RVector3 lineAt( double t ) const { return p0_ + ( p1_ - p0_ ) * t; }

    /*!  Return the distance between this line and pos,
    http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html;
        | x2-x1 cross x1-x0 | / | x2 -x1 | */
    double distance( const RVector3 & pos ) const;

    /*! Return position at line.
        solves: pos = p0_ + ( p1_ - p0_ ) * ret.
        ret = ( pos - p0_ ) / ( p1_ - p0_)
        Throws an exception if pos is not on Line. */
    double t( const RVector3 & pos, double tol = TOLERANCE ) const ;

    /*! Return true if pos touches the line. pFunIdx gives an identifier which shows the dependency between this Line and the RVector3 pos. Possible return values are: \n
    -1 -- this straight line don't touch the Position pos \n
    1 -- this straight line touch the Position pos and pos lies before p0_ \n
    2 -- this straight line touch the Position pos and pos lies at the same position like pos p0_\n
    3 -- this straight line touch the Position pos and pos lies within the definition positions p0_ and p1_\n
    4 -- this straight line touch the Position pos and pos lies at the same position like pos p1_\n
    5 -- this straight line touch the Position pos and pos lies behind p1_\n  */
    bool touch1( const RVector3 & pos, int & pFunIdx, double tol = TOLERANCE ) const ;

    /*! Return true if pos touches the line. */
    bool touch( const RVector3 & pos, double tol = TOLERANCE ) const ;

protected:

    /*! internal: copy content of line into this line */
    void copy_( const Line & line );

    RVector3 p0_;
    RVector3 p1_;
    bool valid_;
};

std::ostream & operator << ( std::ostream & str, const Line & l );

} // namespace GIMLI;

#endif // GIMLI_LINE__H
