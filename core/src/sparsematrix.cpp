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
SparseMatrix< double >::add(const ElementMatrix< double > & A, const double & scale, bool neg){
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
        A.integrate();
        for (Index i = 0, imax = A.rows(); i < imax; i++){
            for (Index j = 0, jmax = A.cols(); j < jmax; j++){
                // __MS(A.rowIDs()[i] << " " << A.colIDs()[j] << "  "
                //       << scale << " " << A.getVal(i, j))
                addVal(A.rowIDs()[i], A.colIDs()[j], b * A.getVal(i, j));
            }
        }
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

} // namespace GIMLI{


