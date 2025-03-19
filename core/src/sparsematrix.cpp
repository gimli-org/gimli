/******************************************************************************
 *   Copyright (C) 2007-2024 by the GIMLi development team                    *
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
#include <algorithm>

namespace GIMLI{

template <> void
SparseMatrix< double >::copy_(const SparseMapMatrix< double, Index > & S){
    this->clear();
    Index col = 0, row = 0;
    double val;
    _cols = S.cols();
    _rows = S.rows();
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
    _cols = S.cols();
    _rows = S.rows();
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
                            const double & f, const double & scale){
    double b = f*scale;

    if (A.oldStyle()){
        if (!valid_) SPARSE_NOT_VALID;
        for (Index i = 0, imax = A.size(); i < imax; i++){
            for (Index j = 0, jmax = A.size(); j < jmax; j++){
                addVal(A.idx(i), A.idx(j), b * A.getVal(i, j));
            }
        }
    } else {
        this->addS(A, f, scale);
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
SparseMatrix< double >::add(const ElementMatrix< double > & A, const Pos & f, const double & scale){
    THROW_TO_IMPL
}
template <> void
SparseMatrix< double >::add(const ElementMatrix< double > & A,
                            const RSmallMatrix & f, const double & scale){
    THROW_TO_IMPL
}

template <> void
SparseMatrix< Complex >::add(const ElementMatrix < double > & A,
                             const Complex & f, const double & scale){
    if (!valid_) SPARSE_NOT_VALID;
    A.integrate();

    for (Index i = 0, imax = A.size(); i < imax; i++){
        for (Index j = 0, jmax = A.size(); j < jmax; j++){
            addVal(A.idx(i), A.idx(j), f*scale * A.getVal(i, j));
        }
    }
}
template <> void
SparseMatrix< Complex >::add(const ElementMatrix < double > & A, const Pos & f, const double & scale){
    THROW_TO_IMPL
}
template <> void
SparseMatrix< Complex >::add(const ElementMatrix < double > & A,
                             const CSmallMatrix & f, const double & scale){
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


template <class ValueType > void
_T_buildSparsityPattern_(SparseMatrix< ValueType > * self,
                        const Mesh & mesh){
    Stopwatch sw(true);

    // self->colPtr_.resize(mesh.nodeCount() + 1);

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

        self->resize(mesh.nodeCount(), mesh.nodeCount());

        // __MS("collect pattern:", sw.duration(true))
        self->buildSparsityPattern(idxMap);
        // __MS("build pattern:", sw.duration())
        // __MS(vals_.size())
        return ;
}


template <> void SparseMatrix< Complex >
::buildSparsityPattern(const Mesh & mesh){
    _T_buildSparsityPattern_(this, mesh);
}


template <> void SparseMatrix< double >
::buildSparsityPattern(const Mesh & mesh){
    _T_buildSparsityPattern_(this, mesh);
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

    self->resize(colPtr.size() - 1, (Index)max(rowIdx) + 1);
    // if (colPtr.size() - 1 != self->rows() ||
    //     (Index)max(rowIdx) + 1 != self->cols()){

    //     __MS(colPtr.size() - 1, self->rows(), max(rowIdx) + 1, self->cols())
    //     log(Error, "build Sparsity Pattern failed!");
    // }
    // _rows = colPtr_.size() - 1;
    // _cols = max(rowIdx_) + 1;
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
template <class ValueType > void
_T_addSparsityPattern_(SparseMatrix< ValueType > * self,
                       const std::vector < std::set< Index > > & idxMap){
    // optimize if necessary
    // __M
    SparseMapMatrix< ValueType, Index > A1(*self);
    // print("self", self->rows(), self->cols());
    // print("A1", A1.rows(), A1.cols());
    SparseMatrix< ValueType > B2;
    B2.buildSparsityPattern(idxMap);
    // print("B2", B2.rows(), B2.cols());
    SparseMapMatrix< ValueType, Index > A2(B2);
    // print("A2", A2.rows(), A2.cols());
    A1 += A2;
    // print("A1", A1.rows(), A1.cols());
    self->copy_(A1);
    // print("self", self->rows(), self->cols());
}

template <> void SparseMatrix< double >
::addSparsityPattern(const std::vector < std::set< Index > > & idxMap){
    _T_addSparsityPattern_(this, idxMap);
}
template <> void SparseMatrix< Complex >
::addSparsityPattern(const std::vector < std::set< Index > > & idxMap){
    _T_addSparsityPattern_(this, idxMap);
}

template <class ValueType > void
_T_reduce_(SparseMatrix< ValueType > * self,
           const IVector & ids, bool keepDiag){

    for (Index i = 0; i < ids.size(); i ++){
        self->cleanRow(ids[i]);
        self->cleanCol(ids[i]);
    }

    // // optimize if necessary
    // SparseMapMatrix< ValueType, Index > A1(*self);
    // A1.reduce(ids, keepDiag);

    // if (keepDiag){
    //     // create diagonal entry, they will be needed probably
    //     for (Index i = 0; i < A1.rows(); i ++){
    //         A1[i][i] += .0;

    //     }
    // }
    // self->copy_(A1);
}

template <> void SparseMatrix< double >
::reduce(const IVector & ids, bool keepDiag){
    _T_reduce_(this, ids, keepDiag);
}
template <> void SparseMatrix< Complex >
::reduce(const IVector & ids, bool keepDiag){
    _T_reduce_(this, ids, keepDiag);
}

template <class ValueType >
IndexArray _T_createReduceMask_impl(SparseMatrix < ValueType > * self,
                                    const IVector & ids){
    IndexArray mask;

    {WITH_TICTOC("rows")
    for (auto row: ids){
        ASSERT_RANGE(row, 0, (int)self->rows())
        for (int ptr = self->colPtr()[row];
                 ptr < self->colPtr()[row + 1]; ptr ++){

            mask.push_back(ptr);
        }
    }
    }
    // TODO: optimize, try CRS->CCC

    {WITH_TICTOC("cols")
    for (auto rID: ids){
        int i = 0;
        for (auto row: self->vecRowIdx()){
            if (row == rID){
                mask.push_back(i);
            }
            i++;
        }
    }
    }

    // for (auto rID: ids){
    //     ASSERT_RANGE(rID, 0, (int)self->rows())
    //     ASSERT_RANGE(rID, 0, (int)self->cols())

    //     int j = 0;
    //     for (int i = 0; i < (int)self->rows(); i ++){
    //         for (int ptr = self->colPtr()[i];
    //                  ptr < self->colPtr()[i + 1]; ptr ++){

    //             j = self->vecRowIdx()[ptr];
    //             //print(i, j, self->values()[ptr]);
    //             if (rID == i || rID == j){
    //                 if (!keepDiag || i != j){
    //                     mask.push_back(ptr);
    //                 }
    //             }
    //         }
    //     }
    // }


    // std::sort(mask.begin(), mask.end());

    return mask;
}

template <> IndexArray SparseMatrix< double >
::createReduceMask(const IVector & ids){
    return _T_createReduceMask_impl(this, ids);
}

template <> IndexArray SparseMatrix< Complex >
::createReduceMask(const IVector & ids){
    return _T_createReduceMask_impl(this, ids);
}

template <class ValueType >
/**
 * @brief Creates a mask for the diagonal elements of a sparse matrix.
 *
 * This function iterates over the rows and columns of the sparse matrix
 * to identify the diagonal elements and stores their indices in a mask.
 *
 * Only diagonal values of the sparsity pattern are considered.
 *
 * @tparam ValueType The type of the values stored in the sparse matrix.
 * @param self Pointer to the sparse matrix object.
 * @return IndexArray A vector containing the indices of the diagonal elements.
 */
IndexArray _T_createDiagonalMask_impl(SparseMatrix < ValueType > * self){
    IndexArray mask;

    for (Index row = 0; row < self->rows(); row ++){
        for (int col = self->colPtr()[row];
                 col < self->colPtr()[row + 1]; col ++){
            if (self->vecRowIdx()[col] == (int)row){
                mask.push_back(col);
            }
        }
    }
    return mask;
}

template <> IndexArray SparseMatrix< double >
::createDiagonalMask(){
    return _T_createDiagonalMask_impl(this);
}

template <> IndexArray SparseMatrix< Complex >
::createDiagonalMask(){
    return _T_createDiagonalMask_impl(this);
}

template <class ValueType> void
mult_T_impl(const SparseMatrix< ValueType > & A,
            const Vector < ValueType > & b, Vector < ValueType > & c,
            const ValueType & alpha, const ValueType & beta,
            Index bOff, Index cOff, bool trans) {
        // c = alpha * (A*b) + beta * c

    if (trans){
        ASSERT_GREATER_EQUAL(b.size() + bOff, A.nRows())
        if (c.size() < A.nCols() + cOff) c.resize(A.nCols() + cOff);
    } else {
        ASSERT_GREATER_EQUAL(b.size() + bOff, A.nCols())
        if (c.size() < A.nRows() + cOff) c.resize(A.nRows() + cOff);
    }
    c *= beta;
    // Index count = 0;

    ValueType si = 0.0;
    ValueType bi = 0.0;

    if (A.stype() == 0){
        // for each row

        for (Index i = 0, iMax = A.rows(); i < iMax; i++){
            if (trans){
                bi = b[i];
                for (int j = A.vecColPtr()[i]; j < A.vecColPtr()[i + 1]; j ++){
                    c[A.vecRowIdx()[j]] += alpha * bi * A.vecVals()[j];
                }
            } else {
                si = c[i];
                for (int j = A.vecColPtr()[i], jMax=A.vecColPtr()[i+1];
                     j < jMax; j ++){
                    si += alpha * b[A.vecRowIdx()[j]] * A.vecVals()[j];
                    // count ++;
                }
                c[i] = si;
            }
        }
        return;
        // print('3', sw.duration(true), count);
    } else if (A.stype() == -1){
        if (trans){
            THROW_TO_IMPL // is needed?
        }
        for (Index i = 0; i < A.rows(); i++){
            si = c[i];
            for (int j = A.vecColPtr()[i]; j < A.vecColPtr()[i + 1]; j ++){
                Index J = A.vecRowIdx()[j];
                ValueType aij(A.vecVals()[j]);

                si += alpha * b[J] * conj(aij);

                if (J > i){
                    c[J] += alpha * b[i] * aij;
                }
            }
            c[i] = si;
        }
    } else if (A.stype() == 1){
        if (trans){
            THROW_TO_IMPL // is needed?
        }
        // symmetric upper part
        for (Index i = 0; i < A.rows(); i++){
            si = c[i];
            for (int j = A.vecColPtr()[i]; j < A.vecColPtr()[i + 1]; j ++){
                Index J = A.vecRowIdx()[j];
                ValueType aij(A.vecVals()[j]);

                si += alpha * b[J] * conj(aij);

                if (J < i){
                    c[J] += alpha * b[i] * aij;
                }
            }
            c[i] = si;
        }
    }
}


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


