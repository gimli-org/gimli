/******************************************************************************
 *   Copyright (C) 2007-2017 by the GIMLi development team                    *
 *   Carsten Rücker carsten@resistivity.net                                   *
 *                                                                            *
 *   Licensed under the Apache License, Version 2.0 (the "License");          *
 *   you may not use this file except in compliance with the License.         *
 *   You may obtain a copy of the License at                                  *
 *                                                                            *
 *       http://www.apache.org/licenses/LICENSE-2.0                           *
 *                                                                            *
 *   Unless required by applicable law or agreed to in writing, software      *
 *   distributed under the License is distributed on an "AS IS" BASIS,        *
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 *   See the License for the specific language governing permissions and      *
 *   limitations under the License.                                           *
 *                                                                            *
 ******************************************************************************/

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


