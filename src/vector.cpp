/***************************************************************************
 *   Copyright (C) 2007-2017 by the GIMLi development team       *
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

#include "gimli.h"
#include "vector.h"
#include "elementmatrix.h"

namespace GIMLI{

template<>
void Vector<double>::add(const ElementMatrix < double >& A){
    for (Index i = 0, imax = A.size(); i < imax; i++){
        data_[A.idx(i)] += A.row(0)[i];
    }
}

template <>
void Vector<double>::add(const ElementMatrix < double >& A, const double & a){
    for (Index i = 0, imax = A.size(); i < imax; i++){
        data_[A.idx(i)] += A.row(0)[i] * a;
    }
}

template <>
void Vector<double>::add(const ElementMatrix < double >& A, const RVector & a){
    for (Index i = 0, imax = A.size(); i < imax; i++){
        data_[A.idx(i)] += A.row(0)[i] * a[A.idx(i)];
    }
}



} // namespace GIMLI{


