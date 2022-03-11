/******************************************************************************
 *   Copyright (C) 2007-2021 by the GIMLi development team                    *
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
#include "sparsemapmatrix.h"
#include "sparsematrix.h"

namespace GIMLI{

template <> void
SparseMapMatrix< double, Index >::copy_(const SparseMatrix< double > & S){
    clear();
    cols_ = S.cols();
    rows_ = S.rows();
    stype_ = S.stype();

    // std::vector < int > colPtr(S.vecColPtr());
    // std::vector < int > rowIdx(S.vecRowIdx());
    // Vector < ValueType > vals(S.vecVals());

    for (Index i = 0; i < S.rows(); i++){
        for (int j = S.vecColPtr()[i]; j < S.vecColPtr()[i + 1]; j ++){
            (*this)[i][S.vecRowIdx()[j]] = S.vecVals()[j];
        }
    }
}
template<>
void SparseMapMatrix< Complex, Index >::copy_(const SparseMatrix< Complex > & S){
    THROW_TO_IMPL
}

template <> void 
SparseMapMatrix< double, Index >::add(const ElementMatrix < double > & A, 
                                      const double & scale, bool neg){
    A.integrate();
    double tol = 1e-25;
    double b = scale;
    if (neg == true) b *= -1.0;

    for (Index i = 0, imax = A.rows(); i < imax; i++){
        for (Index j = 0, jmax = A.mat().cols(); j < jmax; j++){
            double v = A.getVal(i, j) * b;
            // if (isnan(v)){
            //     THROW_TO_IMPL
            // }
            // __MS(A.rowIDs()[i], A.colIDs()[j], v)
            this->addVal(A.rowIDs()[i], A.colIDs()[j], v);

            if (::fabs(v) > tol){
                // __MS(A.rowIDs()[i] << " " << A.colIDs()[j] << " " << v * scale)
            }
        }
    }
}
template <> void 
SparseMapMatrix< double, Index >::add(const ElementMatrix < double > & A, 
                                      const Pos & scale, bool neg){
    A.integrate();
    THROW_TO_IMPL
    // double tol = 1e-25;
    // for (Index i = 0, imax = A.rows(); i < imax; i++){
    //     for (Index j = 0, jmax = A.mat().cols(); j < jmax; j++){
    //         double v = A.getVal(i, j) * scale;
    //             this->addVal(A.rowIDs()[i], A.colIDs()[j], v);
    //         if (::fabs(v) > tol){
    //             // __MS(A.rowIDs()[i] << " " << A.colIDs()[j] << " " << v * scale)
    //         }
    //     }
    // }
}
template <> void 
SparseMapMatrix< double, Index >::add(const ElementMatrix < double > & A, 
                                      const RMatrix & scale, bool neg){
    A.integrate();
    THROW_TO_IMPL
    // double tol = 1e-25;
    // for (Index i = 0, imax = A.rows(); i < imax; i++){
    //     for (Index j = 0, jmax = A.mat().cols(); j < jmax; j++){
    //         double v = A.getVal(i, j) * scale;
    //             this->addVal(A.rowIDs()[i], A.colIDs()[j], v);
    //         if (::fabs(v) > tol){
    //             // __MS(A.rowIDs()[i] << " " << A.colIDs()[j] << " " << v * scale)
    //         }
    //     }
    // }
}

template <> void 
SparseMapMatrix< double, Index >::add(const ElementMatrix < double > & A, 
                                      const Vector < double > & scale, bool neg){
    A.integrate();
    __MS("inuse?")
    //neg unused
    double tol = 1e-25;
    for (Index i = 0, imax = A.rows(); i < imax; i++){
        for (Index j = 0, jmax = A.mat().cols(); j < jmax; j++){
            double v = A.getVal(i, j) * scale[i];
                this->addVal(A.rowIDs()[i], A.colIDs()[j], v);
            if (::fabs(v) > tol){
            }
        }
    }
}
template <> void 
SparseMapMatrix< Complex, Index >::add(const ElementMatrix < double > & A, 
                                       const Complex & scale, bool neg){
    A.integrate();
    // __MS("inuse?")
    //neg unused
    for (Index i = 0, imax = A.rows(); i < imax; i++){
        for (Index j = 0, jmax = A.mat().cols(); j < jmax; j++){
            double v = A.getVal(i, j);
            (*this)[A.rowIDs()[i]][A.colIDs()[j]] += v * scale;
        }
    }
}
template <> void 
SparseMapMatrix< Complex, Index >::add(const ElementMatrix < double > & A, 
                                       const Vector < Complex > & scale, bool neg){
    A.integrate();
    // __MS("inuse?")
    //neg unused
    for (Index i = 0, imax = A.rows(); i < imax; i++){
        for (Index j = 0, jmax = A.mat().cols(); j < jmax; j++){
            double v = A.getVal(i, j);
            (*this)[A.rowIDs()[i]][A.colIDs()[j]] += v * scale[i];
        }
    }
}
template <> void
SparseMapMatrix< Complex, Index >::add(const ElementMatrix < double > & A, 
                                       const CMatrix & scale, bool neg){
    A.integrate();
    THROW_TO_IMPL
}
template <> void 
SparseMapMatrix< Complex, Index >::add(const ElementMatrix < double > & A, 
                                       const Pos & scale, bool neg){
    A.integrate();
    THROW_TO_IMPL
}

template <> void SparseMapMatrix< double, Index >::
    addToCol(Index id, const ElementMatrix < double > & A, double scale, bool isDiag){
    A.integrate();
    for (Index i = 0, imax = A.size(); i < imax; i++){
            if (isDiag){
                (*this)[A.idx(i)][id] += A.getVal(i, i);
            } else {
                (*this)[A.idx(i)][id] += A.getVal(0, i);
            }
    }
}

template <> void SparseMapMatrix< double, Index >::
    addToRow(Index id, const ElementMatrix < double > & A, double scale, bool isDiag){
    A.integrate();

    for (Index i = 0, imax = A.size(); i < imax; i++){
        if (isDiag){
            (*this)[id][A.idx(i)] += A.getVal(i, i);
        } else {
            (*this)[id][A.idx(i)] += A.getVal(0, i);
        }
    }
}

template <> void SparseMapMatrix< Complex, Index >::
    addToCol(Index id, const ElementMatrix < double > & A,
             Complex scale, bool isDiag){
        A.integrate();
THROW_TO_IMPL
}
template <> void SparseMapMatrix< Complex, Index >::
    addToRow(Index id, const ElementMatrix < double > & A,
             Complex scale, bool isDiag){
        A.integrate();
THROW_TO_IMPL
}

template <class ValueType> void
mult_(const SparseMapMatrix< ValueType, Index > * A,
      const Vector < ValueType > & a, Vector < ValueType > & ret) {

    if (A->rows() != ret.size()) ret.resize(A->rows());
    ASSERT_GREATER_EQUAL(a.size(), A->cols())

    if (A->stype() == 0){
        for (auto it = A->begin(); it != A->end(); it ++){
            ret[it->first.first] += a[it->first.second] * it->second;
        }
    } else if (A->stype() == -1){

        for (auto it = A->begin(); it != A->end(); it ++){
            Index I = it->first.first;
            Index J = it->first.second;

            ret[I] += a[J] * conj(it->second);

            if (J > I){
                ret[J] += a[I] * it->second;
            }
        }
    } else if (A->stype() ==  1){
        for (auto it = A->begin(); it != A->end(); it ++){
            Index I = it->first.first;
            Index J = it->first.second;

            ret[I] += a[J] * conj(it->second);

            if (J < I){
                ret[J] += a[I] * it->second;
            }
        }
    }
}

template <> void SparseMapMatrix< double, Index >::
    mult(const Vector < double > & a, Vector < double > & ret) const{
    mult_(this, a, ret);
}
template <> void SparseMapMatrix< Complex, Index >::
    mult(const Vector < Complex > & a, Vector < Complex > & ret) const{
    mult_(this, a, ret);
}
template <> void SparseMapMatrix< double, Index >::
    mult(const Vector < double > & a, Vector < Pos > & ret) const {
    if (this->rows() != ret.size()) ret.resize(this->rows(), Pos(0.0, 0.0));

    Index nCoeff(a.size() / this->cols());
    Index dof(this->cols());
    ASSERT_GREATER_EQUAL(a.size(), nCoeff*dof)

    if (stype_ == 0){
        for (const_iterator it = this->begin(); it != this->end(); it ++){
            for (Index i = 0; i < nCoeff; i ++){
                ret[it->first.first][i] +=
                                a[it->first.second + i * dof] * it->second;
            }
        }
    } else if (stype_ == -1){
        THROW_TO_IMPL
    } else if (stype_ ==  1){
        THROW_TO_IMPL
    }
}
template <> void SparseMapMatrix< Complex, Index >::
    mult(const Vector < Complex > & a, Vector < Pos > & ret) const {
THROW_TO_IMPL
}


template < class ValueType >
void transMult_(const SparseMapMatrix< ValueType, Index > * A,
                const Vector < ValueType > & a,
                Vector < ValueType > & ret) {
    if (A->cols() != ret.size()) ret.resize(A->cols());

    ASSERT_GREATER_EQUAL(a.size(), A->rows())

    if (A->stype() == 0){
        for (auto it = A->begin(); it != A->end(); it ++){
            ret[it->first.second] += a[it->first.first] * it->second;
        }
    } else if (A->stype() == -1){
        THROW_TO_IMPL
    } else if (A->stype() ==  1){
        THROW_TO_IMPL
    }
}
template <> void SparseMapMatrix< double, Index >::
    transMult(const Vector < double > & a, Vector < double > & ret) const{
    transMult_(this, a, ret);
}
template <> void SparseMapMatrix< Complex, Index >::
    transMult(const Vector < Complex > & a, Vector < Complex > & ret) const{
    THROW_TO_IMPL
    //transMult_(this, a, ret);
}
template <> void SparseMapMatrix< double, Index >::
    transMult(const Vector < double > & a, Vector < Pos > & ret) const{
    THROW_TO_IMPL
}
template <> void SparseMapMatrix< Complex, Index >::
    transMult(const Vector < Complex > & a, Vector < Pos > & ret) const{
    THROW_TO_IMPL
}

template <> DLLEXPORT void SparseMapMatrix< double, Index >
::reduce(const IVector & ids, bool keepDiag) {

    for (auto it = begin(); it != end();){

        auto r1 = std::find(&ids[0], &ids[ids.size()], it->first.first);
        if (r1 != &ids[ids.size()]){

            if (it->first.first != it->first.second){
                it = C_.erase(it);
                continue;
            } else {
                if (!keepDiag){
                    it = C_.erase(it);
                    continue;
                }else {
                    ++it;
                    continue;
                }
            }
        }
        auto r2 = std::find(&ids[0], &ids[ids.size()], it->first.second);
        if (r2 != &ids[ids.size()]){
            it = C_.erase(it);
            continue;
        }
        ++it;
    }
}
template <> DLLEXPORT void SparseMapMatrix< Complex, Index >::
    reduce(const IVector & ids, bool keepDiag) {
    THROW_TO_IMPL
}


void mult(const std::vector < RSparseMapMatrix > & A,
          const RVector & b, std::vector< RVector > & ret){
    //!! implement with ompl .. refactor with above
    if (ret.size() != A.size()) ret.resize(A.size());
    for (Index i = 0; i < A.size(); i ++ ){
        ret[i] *= 0.0;
        A[i].mult(b, ret[i]);
    }
}

void mult(const std::vector < RSparseMapMatrix > & A,
          const RVector & b, std::vector< PosVector > & ret){
    //!! implement with ompl, refactor with funct above
    if (ret.size() != A.size()) ret.resize(A.size());
    for (Index i = 0; i < A.size(); i ++ ){
        ret[i] *= Pos(0.0, 0.0);
        A[i].mult(b, ret[i]);
    }
}

std::vector< RVector > mult(const std::vector < RSparseMapMatrix > & A,
                            const RVector & b){
    std::vector< RVector > ret;
    mult(A, b, ret);
    return ret;
}




} // namespace GIMLI{
