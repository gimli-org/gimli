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

#include "numericbase.h"
#include "vector.h"
#include "pos.h"

namespace GIMLI{

RVector logTransDropTol(const RVector & data, double logdrop){
    RVector tmp(data);

    for (uint i = 0; i < tmp.size(); i ++) {
        tmp[i] =std::fabs(tmp[i] / logdrop);
        if (tmp[i] < 1.0) tmp[i] = 1.0;
    }

    tmp = log10(tmp);
    tmp /= max(tmp) * sign(data);
    return tmp;
}
    
void GaussLaguerre( uint n, RVector & x, RVector & w ){
//  taken from matlab-code (Thomas Guenther)
//  function [x, w] = gaulag(n)
//  GAULAG - Gauss-Laguerre Integration Points
//  [x,w] = gaulag(n)
//  Given alf = 0.0, the parameter alpha of the Laguerre polynomials, this routine
//  returns arrays x[1..n] and w[1..n] containing the abscissas and weights
//  of the n-point Gauss-Laguerre quadrature formula. The smallest abscissa
//  is returned in x[1], the largest in x[n].
//  For a description of the following routines see
//  Numerical Recipes, Press et a

  if ( x.size() != n ) x.resize( n );
  if ( w.size() != n ) w.resize( n );

  double epsilon = 3.0e-11;
  int maxiter = 20;

  double z = 0.0, z1 = 0.0, p1 = 1.0, p2 = 0.0, p3 = 0.0, pp = 0.0;
  int ai = 0;

  for ( size_t i = 1; i <=n; i ++ ){ //	Loop over desired roots
    if ( i == 1 ){
      z = 3.0 / ( 1.0 + 2.4 * n );
    } else if ( i == 2 ){
      z = z + 15.0 / ( 1.0 +2.5 * n );
    } else {
      ai = i - 2;
      z = z + ( 1.0 + 2.55 * ai ) / ( 1.9 * ai ) *( z - x[ ai - 1 ] );
    }

    for ( int its = 1; its <= maxiter; its ++ ){
      p1 = 1.0;
      p2 = 0.0;

      for ( size_t j = 1; j <= n; j ++ ){
	p3 = p2;
	p2 = p1;
	p1 = ( ( 2.0 * j - 1 - z ) * p2 - ( j - 1 ) * p3 ) / j;
      }
      pp = n *( p1 - p2 ) / z;
      z1 = z;
      z  = z1 - p1 / pp;

      if ( std::fabs( z - z1 ) <= epsilon ) break;
    }
    x[ i - 1 ] = z;
    w[ i - 1 ] = -1.0 / ( pp * n * p2 );
  }
}

void GaussLegendre( double x1, double x2, uint n, RVector & x, RVector & w ){
//  taken from matlab-code (Thomas Guenther)
//  function [x, w] = gauleg(x1, x2, n)
//  Given the lower and upper limits of integration x1 and x2 and given n,
//  this routine returns arrays x(1..n) and w(1..n) of length n,
//  containing the abscissas and weights of the Gauss-Legendre
//  n-point quadrature formula.
//  For a description of the following routines see
//  Numerical Recipes, Press et al.

  if ( x.size() != n ) x.resize( n );
  if ( w.size() != n ) w.resize( n );

  double epsilon = 3.0e-6;

  double m = ( n + 1.0 ) / 2.0 ;
  double xm = 0.5 * ( x2 + x1 );
  double xl = 0.5 * ( x2 - x1 );

  double z = 0.0, z1 = 0.0, p1 = 0.0, p2 = 0.0, p3 = 0.0, pp = 0.0;

  for ( int i = 1; i <= m; i ++ ){
    z = std::cos( PI * ( i - 0.25 ) / ( n + 0.5 ) );

   // Starting with the above approximation to the ith root, we enter
   // the main loop of refinements by Newton's method
   z1 = z + 2.0 * epsilon;

   while ( std::fabs( z - z1 ) > epsilon ){
     p1 = 1.0;
     p2 = 0.0;
     for ( size_t j = 1; j <= n; j ++ ){
       p3 = p2;
       p2 = p1;
       p1 = ( ( 2.0 * j - 1.0 ) * z * p2 - ( j - 1.0 ) * p3 ) / (double)j;
     }

     // p1 is now the desired Legendre polynomial. We next compute pp,
     // its derivative, by a standard relation involving also p2, the
     // polynomial of one lower order
     pp = (double)n * ( z * p1 - p2 ) / ( z * z - 1.0 );
     z1 = z;
     z = z1 - p1 / pp; // Newtons method
   }
   // Scale the root to the desired interval, and put in its
   // symmetric counterpart
   x[ i - 1 ] = xm - xl * z;
   x[ n - i ] = xm + xl * z;
   //   x[ n +1 - i ] = xm + xl * z;
   //Compute the weight and ist symmetric counterpart
   w[ i - 1 ] = 2.0 * xl / ( ( 1.0 -z * z ) *pp *pp );
   w[ n - i ] = w[ i - 1 ];
   //   w[ n + 1 - i ] = w[ i - 1 ];
  }
}

RVector3 sphTangential2Initerial(const RVector3 &V, double lat, double lon){
    
    double th = lat * PI/180.0;
    double ph = lon * PI/180.0;
    double ct = std::cos(th);
    double st = std::sin(th);
    double cp = std::cos(ph);
    double sp = std::sin(ph);
    return RVector3(
            (V[0]* ct + V[1] * -st) * cp + V[2] * -sp,
            (V[0]* ct + V[1] * -st) * sp + V[2] *  cp,
             V[0]* st + V[1] * ct);
}


} //namespace GIMLI;
