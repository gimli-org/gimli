/******************************************************************************
 *   Copyright (C) 2007-2022 by the GIMLi development team                    *
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
#include "sparsematrix.h"

namespace GIMLI{


// template <> void
// SparseMatrix(const IndexArray & colPtr,
//              const IndexArray & rowIdx,
//              const Vector < double > vals, int stype=0)
//         : SparseMatrixBase(){
//         colPtr_ = std::vector < int >(colPtr.size());
//         rowIdx_ = std::vector < int >(rowIdx.size());
//         for (Index i = 0; i < colPtr_.size(); i ++ ) colPtr_[i] = colPtr[i];
//         for (Index i = 0; i < rowIdx_.size(); i ++ ) rowIdx_[i] = rowIdx[i];
//         vals_   = vals;
//         stype_  = stype;
//         valid_  = true;
//         cols_ = max(rowIdx_) + 1;
//         rows_ = colPtr_.size() - 1;
//     }

// template <> void
// SparseMatrix(const std::vector < int > & colPtr,
//              const std::vector < int > & rowIdx,
//              const Vector < double > vals, int stype=0)
//         : SparseMatrixBase(){
//           colPtr_ = colPtr;
//           rowIdx_ = rowIdx;
//           vals_   = vals;
//           stype_  = stype;
//           valid_  = true;
//           cols_ = max(rowIdx_) + 1;
//           rows_ = colPtr_.size() - 1;
//     }



template <> void
SparseMatrix< double >::copy_(const SparseMapMatrix< double, Index > & S){
    this->clear();
    Index col = 0, row = 0;
    double val;
    cols_ = S.cols();
    rows_ = S.rows();
    std::vector < std::map < Index, double > > rowColMap(S.rows());

    for (typename SparseMapMatrix< double, Index>::const_iterator
        it = S.begin(); it != S.end(); it ++){
        row = S.idx1(it);
        col = S.idx2(it);
        val = S.val(it);
        rowColMap[row].insert(std::pair< Index, double >(col, val));
    }
    colPtr_.resize(S.rows() + 1);
    rowIdx_.resize(S.nVals());
    vals_.resize(S.nVals());
    stype_  = S.stype();

    this->colPtr_[0] = 0;

    Index colPtr = 0;
    row = 0;
    for (typename std::vector < std::map < Index, double > >::iterator
         it = rowColMap.begin(); it != rowColMap.end(); it++){
        for (typename std::map < Index, double >::iterator
             itR = (*it).begin(); itR != (*it).end(); itR++){
            rowIdx_[colPtr] = itR->first;
            vals_[colPtr] = (double)itR->second;
            colPtr ++;
        }
        row ++;
        this->colPtr_[row] = colPtr;
    }
    valid_ = true;
}

template<>
void SparseMatrix< Complex >::copy_(const SparseMapMatrix< Complex, Index > & S){
    this->clear();
    Index col = 0, row = 0;
    Complex val;
    cols_ = S.cols();
    rows_ = S.rows();
    std::vector < std::map < Index, Complex > > rowColMap(S.rows());

    for (typename SparseMapMatrix< Complex, Index>::const_iterator
        it = S.begin(); it != S.end(); it ++){
        row = S.idx1(it);
        col = S.idx2(it);
        val = S.val(it);
        rowColMap[row].insert(std::pair< Index, Complex >(col, val));
    }
    colPtr_.resize(S.rows() + 1);
    rowIdx_.resize(S.nVals());
    vals_.resize(S.nVals());
    stype_  = S.stype();

    this->colPtr_[0] = 0;

    Index colPtr = 0;
    row = 0;
    for (typename std::vector < std::map < Index, Complex > >::iterator
         it = rowColMap.begin(); it != rowColMap.end(); it++){
        for (typename std::map < Index, Complex >::iterator
             itR = (*it).begin(); itR != (*it).end(); itR++){
            rowIdx_[colPtr] = itR->first;
            vals_[colPtr] = (Complex)itR->second;
            colPtr ++;
        }
        row ++;
        this->colPtr_[row] = colPtr;
    }
    valid_ = true;
}

template <> void
SparseMatrix< double >::add(const ElementMatrix< double > & A,
                            const double & scale, bool neg){
    double b = scale;
    if (neg == true) b *= -1.0;

    if (A.oldStyle()){
        if (!valid_) SPARSE_NOT_VALID;
        for (Index i = 0, imax = A.size(); i < imax; i++){
            for (Index j = 0, jmax = A.size(); j < jmax; j++){
                addVal(A.idx(i), A.idx(j), b * A.getVal(i, j));
            }
        }
    } else {
        this->addS(A, scale, neg);
        // A.integrate();
        // for (Index i = 0, imax = A.rows(); i < imax; i++){
        //     for (Index j = 0, jmax = A.cols(); j < jmax; j++){
        //         // __MS(A.rowIDs()[i] << " " << A.colIDs()[j] << "  "
        //         //       << scale << " " << A.getVal(i, j))
        //         addVal(A.rowIDs()[i], A.colIDs()[j], b * A.getVal(i, j));
        //     }
        // }
    }
}
template <> void
SparseMatrix< double >::add(const ElementMatrix< double > & A, const Pos & scale, bool neg){
    THROW_TO_IMPL
}
template <> void
SparseMatrix< double >::add(const ElementMatrix< double > & A, const Matrix < double > & scale, bool neg){
    THROW_TO_IMPL
}

template <> void
SparseMatrix< Complex >::add(const ElementMatrix < double > & A, const Complex & scale, bool neg){
    if (!valid_) SPARSE_NOT_VALID;

    for (Index i = 0, imax = A.size(); i < imax; i++){
        for (Index j = 0, jmax = A.size(); j < jmax; j++){
            if (neg == true) {
                subVal(A.idx(i), A.idx(j), scale * A.getVal(i, j));
            } else {
                addVal(A.idx(i), A.idx(j), scale * A.getVal(i, j));
            }
        }
    }
}
template <> void
SparseMatrix< Complex >::add(const ElementMatrix < double > & A, const Pos & scale, bool neg){
    THROW_TO_IMPL
}
template <> void
SparseMatrix< Complex >::add(const ElementMatrix < double > & A, const CMatrix & scale, bool neg){
    THROW_TO_IMPL
}


template <> void SparseMatrix< double >
::setVal(Index i, Index j, const double & val){
    if (abs(val) > TOLERANCE || 1){
        for (int k = colPtr_[i]; k < colPtr_[i + 1]; k ++){
            if (rowIdx_[k] == (int)j) {
                vals_[k] = val; return;
            }
        }
        __M
        log(Error, " setVal: " ,i,  " ", j, " is not part of the sparsity pattern ");
    }
}

template <> void SparseMatrix< Complex >
::setVal(Index i, Index j, const Complex & val){
    if (abs(val) > TOLERANCE || 1){
        for (int k = colPtr_[i]; k < colPtr_[i + 1]; k ++){
            if (rowIdx_[k] == (int)j) {
                vals_[k] = val; return;
            }
        }
        __M
        log(Error, " setVal: " ,i,  " ", j, " is not part of the sparsity pattern ");
    }
}

template <> double SparseMatrix< double >
::getVal(Index i, Index j, bool warn) const{
    for (int k = colPtr_[i]; k < colPtr_[i + 1]; k ++){
        if (rowIdx_[k] == (int)j) {
            return vals_[k];
        }
    }
    if (warn == true){
        __M
        log(Error, " getVal: " , i, " ", j, " is not part of the sparsity pattern ");
    }
    return 0.0;
}

template <> Complex SparseMatrix< Complex >
::getVal(Index i, Index j, bool warn) const {
    for (int k = colPtr_[i]; k < colPtr_[i + 1]; k ++){
        if (rowIdx_[k] == (int)j) {
            return vals_[k];
        }
    }
    if (warn == true){
        __M
        log(Error, " getVal: " , i, " ", j, " is not part of the sparsity pattern ");
    }
    return Complex(0);
}

template <class ValueType > void
_T_fillMassMatrix(SparseMatrix< ValueType > * self,
                  const Mesh & mesh, const RVector & a, bool rebuildPattern){
    if (rebuildPattern == true || self->size() == 0){
        self->buildSparsityPattern(mesh);
    } else {
        self->clean();
    }
    ElementMatrix < double > A_l;

    for (uint i = 0; i < mesh.cellCount(); i ++){
        A_l.u2(mesh.cell(i));
        A_l *= a[mesh.cell(i).id()];
        *self += A_l;
    }
}

template <> void SparseMatrix< double >
::fillMassMatrix(const Mesh & mesh, const RVector & a, bool rebuildPattern){
    _T_fillMassMatrix(this, mesh, a, rebuildPattern);
}
template <> void SparseMatrix< Complex >
::fillMassMatrix(const Mesh & mesh, const RVector & a, bool rebuildPattern){
    _T_fillMassMatrix(this, mesh, a, rebuildPattern);
}

template <class ValueType > void
_T_fillStiffnessMatrix(SparseMatrix< ValueType > * self,
                       const Mesh & mesh, const RVector & a, bool rebuildPattern){
    if (rebuildPattern == true || self->size() == 0){
        self->buildSparsityPattern(mesh);
    } else {
        self->clean();
    }
    ElementMatrix < double > A_l;

    for (uint i = 0; i < mesh.cellCount(); i ++){
        A_l.ux2uy2uz2(mesh.cell(i));
        A_l *= a[mesh.cell(i).id()];
        *self += A_l;
    }
}

template <> void SparseMatrix< double >
::fillStiffnessMatrix(const Mesh & mesh, const RVector & a, bool rebuildPattern){
    _T_fillStiffnessMatrix(this, mesh, a, rebuildPattern);
}
template <> void SparseMatrix< Complex >
::fillStiffnessMatrix(const Mesh & mesh, const RVector & a, bool rebuildPattern){
    _T_fillStiffnessMatrix(this, mesh, a, rebuildPattern);
}

template <> void SparseMatrix< Complex >
::buildSparsityPattern(const Mesh & mesh){
    THROW_TO_IMPL
}

template <> void SparseMatrix< double >
::buildSparsityPattern(const Mesh & mesh){
    Stopwatch sw(true);

    colPtr_.resize(mesh.nodeCount() + 1);

        //*** much to slow
        //RSparseMapMatrix S(mesh.nodeCount(), mesh.nodeCount());

    Index col = 0, row = 0;

        // need unique(sort) maybe to slow
//        std::vector < std::vector< Index > > idxMap(mesh.nodeCount());
//         for (std::vector < std::vector< Index > >::iterator mIt = idxMap.begin(); mIt != idxMap.end(); mIt++){
//             (*mIt).reserve(100);
//         }

// using set is by a factor of approx 5 more expensive
        std::vector < std::set< Index > > idxMap(mesh.nodeCount());

        Cell *cell = 0;
        uint nc = 0;

        for (uint c = 0; c < mesh.cellCount(); c ++){
            cell = &mesh.cell(c);
            nc = cell->nodeCount();

            for (uint i = 0; i < nc; i ++){
                for (uint j = 0; j < nc; j ++){
                    row = cell->node(i).id();
                    col = cell->node(j).id();
                    //idxMap[col].push_back(row);
                    idxMap[col].insert(row);
                    //S[col][row] = 1;
                }
            }
        }

        rows_ = mesh.nodeCount();
        cols_ = mesh.nodeCount();

        __MS(sw.duration(true))
        this->buildSparsityPattern(idxMap);
        __MS(sw.duration())
        // __MS(vals_.size())
        return ;

//         int nVals = 0;
//         for (std::vector < std::set< Index > >::iterator mIt = idxMap.begin();
//              mIt != idxMap.end(); mIt++){
//             //std::sort((*mIt).begin(), (*mIt).end());
//             nVals += (*mIt).size();
//         }

// //         std::cout << "timwe: " << swatch.duration( true) << std::endl;
// //         exit(0);

//         rowIdx_.reserve(nVals);
//         rowIdx_.resize(nVals);
//         vals_.resize(nVals);

//         colPtr_[0] = 0;
//         Index k = 0;
//         row = 0;
//         for (std::vector < std::set< Index > >::iterator mIt = idxMap.begin(); mIt != idxMap.end(); mIt++){
//             for (std::set< Index >::iterator sIt = (*mIt).begin(); sIt != (*mIt).end(); sIt++){
//                 rowIdx_[k] = (*sIt);
//                 vals_[k] = 0.0;
//                 k++;
//             }
//             row++;
//             colPtr_[row] = k;
//         }
//         valid_ = true;

//         rows_ = colPtr_.size() - 1;
//         cols_ = max(rowIdx_) + 1;
//         //** freeing idxMap is expensive
}

template <class ValueType > void
_T_buildSparsityPattern_(SparseMatrix< ValueType > * self,
                         const std::vector < std::set< Index > > & idxMap){

    Index nVals = 0;
    for (auto & r: idxMap){
        nVals += r.size();
    }

    // for (std::vector < std::set< Index > >::iterator mIt = idxMap.begin();
    //     mIt != idxMap.end(); mIt++){
    //     //std::sort((*mIt).begin(), (*mIt).end());
    //     nVals += (*mIt).size();
    // }

//         std::cout << "timwe: " << swatch.duration( true) << std::endl;
//         exit(0);

    std::vector < int > & colPtr = self->vecColPtr();
    std::vector < int > & rowIdx = self->vecRowIdx();

    colPtr.resize(idxMap.size() + 1);
    rowIdx.resize(nVals);
    self->vecVals().resize(nVals);
    self->vecVals().clean();

    colPtr[0] = 0;
    Index k = 0;
    Index row = 0;

    for (auto & r: idxMap){
        // print(row, "\n");
        for (auto & c: r){
            // print(c);
            rowIdx[k] = c;
            
            //vals_[k] = (ValueType)0.0;
            k++;
        }
        row++;
        colPtr[row] = k;
    }
    self->setValid(true);

    if (colPtr.size() - 1 != self->rows() ||
        (Index)max(rowIdx) + 1 != self->cols()){
        __MS(colPtr.size() - 1, self->rows(), max(rowIdx) + 1, self->cols())
        log(Error, "build Sparsity Pattern failed!");
    }
    // rows_ = colPtr_.size() - 1;
    // cols_ = max(rowIdx_) + 1;
    //** freeing idxMap is expensive

}

template <> void SparseMatrix< double >
::buildSparsityPattern(const std::vector < std::set< Index > > & idxMap){
    _T_buildSparsityPattern_(this, idxMap);
}
template <> void SparseMatrix< Complex >
::buildSparsityPattern(const std::vector < std::set< Index > > & idxMap){
    _T_buildSparsityPattern_(this, idxMap);
}

template <class ValueType> void
mult_T_impl(const SparseMatrix< ValueType > & A,
            const Vector < ValueType > & b, Vector < ValueType > & c,
            const ValueType & alpha, const ValueType & beta,
            Index bOff, Index cOff, bool trans) {

    if (trans){
        ASSERT_GREATER_EQUAL(b.size() + bOff, A.rows())
        if (c.size() < A.cols() + cOff) c.resize(A.cols() + cOff);
    } else {
        ASSERT_GREATER_EQUAL(b.size() + bOff, A.cols())
        if (c.size() < A.rows() + cOff) c.resize(A.rows() + cOff);
    }
    THROW_TO_IMPL
    // c *= beta;

    // if (A.stype() == 0){ // non-symmetric
    //     if (trans){
    //         for (auto it = A.begin(); it != A.end(); it ++){
    //             c[it->first.second] += alpha * b[it->first.first] * it->second;
    //         }
    //     } else {
    //         for (auto it = A.begin(); it != A.end(); it ++){
    //             c[it->first.first] += alpha * b[it->first.second] * it->second;
    //         }
    //     }
    // } else if (A.stype() == -1){
    //     if (trans){
    //         THROW_TO_IMPL
    //     } else {
    //         for (auto it = A.begin(); it != A.end(); it ++){
    //             Index I = it->first.first;
    //             Index J = it->first.second;

    //             c[I] += alpha * b[J] * conj(it->second);

    //             if (J > I){
    //                 c[J] += alpha * b[I] * it->second;
    //             }
    //         }
    //     }
    // } else if (A.stype() ==  1){
    //     if (trans){
    //         THROW_TO_IMPL
    //     } else {
    //         for (auto it = A.begin(); it != A.end(); it ++){
    //             Index I = it->first.first;
    //             Index J = it->first.second;

    //             c[I] += alpha * b[J] * conj(it->second);

    //             if (J < I){
    //                 c[J] += alpha * b[I] * it->second;
    //             }
    //         }
    //     }
    // }
}



// /*! Return this * a  */
//     virtual Vector < ValueType > mult(const Vector < ValueType > & a) const {
//         if (a.size() < this->cols()){
//             throwLengthError(WHERE_AM_I + " SparseMatrix size(): " + str(this->cols()) + " a.size(): " +
//                                 str(a.size())) ;
//         }

//         Vector < ValueType > ret(this->rows(), 0.0);

//         if (stype_ == 0){
//             // for each row
//             for (Index i = 0; i < this->rows(); i++){
//             // iterate through compressed col
//                 for (int j = this->vecColPtr()[i]; j < this->vecColPtr()[i + 1]; j ++){
//                     ret[i] += a[this->vecRowIdx()[j]] * this->vecVals()[j];
//                 }
//             }
//         } else if (stype_ == -1){
//             Index J;
//             for (Index i = 0; i < ret.size(); i++){
//                 for (int j = this->vecColPtr()[i]; j < this->vecColPtr()[i + 1]; j ++){
//                     J = this->vecRowIdx()[j];

// //                     __MS( i << "  " << J << " " << this->vecVals()[j])
//                     ret[i] += a[J] * conj(this->vecVals()[j]);

//                     if (J > i){
// //                         __MS( J << "  " << i << " " << this->vecVals()[j])
//                         ret[J] += a[i] * this->vecVals()[j];
//                     }
//                 }
//             }

//             //#THROW_TO_IMPL
//         } else if (stype_ == 1){
//             Index J;
//             for (Index i = 0; i < ret.size(); i++){
//                 for (int j = this->vecColPtr()[i]; j < this->vecColPtr()[i + 1]; j ++){
//                     J = this->vecRowIdx()[j];

// //                     __MS( i << "  " << J << " " << this->vecVals()[j])
//                     ret[i] += a[J] * conj(this->vecVals()[j]);

//                     if (J < i){
// //                         __MS( J << "  " << i << " " << this->vecVals()[j])
//                         ret[J] += a[i] * this->vecVals()[j];
//                     }
//                 }
//             }
//         }
//         return ret;
//     }

//     /*! Return this.T * a */
//     virtual Vector < ValueType > transMult(const Vector < ValueType > & a) const {

//         if (a.size() < this->rows()){
//             throwLengthError(WHERE_AM_I + " SparseMatrix size(): " + str(this->rows()) + " a.size(): " +
//                                 str(a.size())) ;
//         }

//         Vector < ValueType > ret(this->cols(), 0.0);

//         if (stype_ == 0){
//             for (Index i = 0; i < this->rows(); i++){
//                 for (int j = this->vecColPtr()[i]; j < this->vecColPtr()[i + 1]; j ++){
//                     ret[this->vecRowIdx()[j]] += a[i] * this->vecVals()[j];
//                 }
//             }

//         } else if (stype_ == -1){
//             THROW_TO_IMPL
//         } else if (stype_ ==  1){
//             THROW_TO_IMPL
//         }
//         return ret;
//     }


























void mult(const SparseMatrix< double > & A,
          const RVector & b, RVector & c,
          const double & alpha, const double & beta,
          Index bOff, Index cOff){
    return mult_T_impl(A, b, c, alpha, beta, bOff, cOff, false);
}
void mult(const SparseMatrix< Complex > & A,
          const CVector & b, CVector & c,
          const Complex & alpha, const Complex & beta,
          Index bOff, Index cOff){
    return mult_T_impl(A, b, c, alpha, beta, bOff, cOff, false);
}
void transMult(const SparseMatrix< double > & A,
          const RVector & b, RVector & c,
          const double & alpha, const double & beta,
          Index bOff, Index cOff){
    return mult_T_impl(A, b, c, alpha, beta, bOff, cOff, true);
}
void transMult(const SparseMatrix< Complex > & A,
          const CVector & b, CVector & c,
          const Complex & alpha, const Complex & beta,
          Index bOff, Index cOff){
    return mult_T_impl(A, b, c, alpha, beta, bOff, cOff, true);
}

} // namespace GIMLI{


