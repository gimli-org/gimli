/******************************************************************************
 *   Copyright (C) 2005-2022 by the GIMLi development team                    *
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
#include "blockmatrix.h"

namespace GIMLI{

template < class ValueType, class Mat >
void mult_MV_T(const Mat & A, 
          const Vector < ValueType > & b, Vector < ValueType > & c,
          const ValueType & alpha, const ValueType & beta,
          Index bOff, Index cOff, bool trans=false){

    if (trans){
        ASSERT_GREATER_EQUAL(b.size() + bOff, A.rows())
        c.resize(A.cols() + cOff);
    } else {
        ASSERT_GREATER_EQUAL(b.size() + bOff, A.cols())
        c.resize(A.rows() + cOff);
    }

    c *= beta;

    for (auto &entry: A.entries()){
        if (entry.transpose) {
            THROW_TO_IMPL
        }
        const MatrixBase & mat = A.matRef(entry.matrixID);

        if (trans){
            c.addVal(alpha * mat.transMult(b.getVal(entry.rowStart,
                                            entry.rowStart + mat.rows())) *
                     entry.scale, entry.colStart, entry.colStart + mat.cols());
        } else {
            c.addVal(alpha * mat.mult(b.getVal(entry.colStart,
                                        entry.colStart + mat.cols())) *
                    entry.scale, entry.rowStart, entry.rowStart + mat.rows());
        }
    }
}

void mult(const RBlockMatrix & A, 
          const RVector & b, RVector & c,
          const double & alpha, const double & beta,
          Index bOff, Index cOff){
    mult_MV_T(A, b, c, alpha, beta, bOff, cOff, false);
}

void transMult(const RBlockMatrix & A, 
               const RVector & b, RVector & c, 
               const double & alpha, const double & beta,
               Index bOff, Index cOff){
    mult_MV_T(A, b, c, alpha, beta, bOff, cOff, true);
}

void mult(const CBlockMatrix & A, 
          const CVector & b, CVector & c,
          const Complex & alpha, const Complex & beta,
          Index bOff, Index cOff){
              THROW_TO_IMPL
    // mult_MV_T(A, b, c, alpha, beta, bOff, cOff, false);
}

void transMult(const CBlockMatrix & A, 
               const CVector & b, CVector & c, 
               const Complex & alpha, const Complex & beta,
               Index bOff, Index cOff){
                   THROW_TO_IMPL
    // mult_MV_T(A, b, c, alpha, beta, bOff, cOff, true);
}

template <> RSparseMapMatrix 
BlockMatrix< double >::sparseMapMatrix() const {

    RSparseMapMatrix ret(this->rows(), this->cols());

    for (Index i = 0; i < entries_.size(); i++){
        auto entry(entries_[i]);

        MatrixBase *mat = matrices_[entry.matrixID];

        RVector vals(0);
        IndexArray rows(0);
        IndexArray cols(0);

        switch (mat->rtti()){
            case GIMLI_SPARSE_CRS_MATRIX_RTTI:{
                RSparseMapMatrix S(*dynamic_cast< RSparseMatrix * >(mat));
                S.fillArrays(vals, rows, cols);
            }  break;
            case GIMLI_SPARSE_MAP_MATRIX_RTTI:
                dynamic_cast< RSparseMapMatrix * >(mat)->fillArrays(vals, rows, cols);
                break;
            default:
                log(Critical, "Matrix type need to be either SparseMatrix or SparseMapMatrix");
                return ret;
        }
        ret.add(rows + entry.rowStart, cols + entry.colStart, vals * entry.scale);
    }
    return ret;
}

template <> CSparseMapMatrix 
BlockMatrix< Complex >::sparseMapMatrix() const {
    THROW_TO_IMPL
    return CSparseMatrix();
}
template <> RSparseMatrix 
BlockMatrix< double >::sparseMatrix() const {
    log(Warning, "Efficiency check needed.");
    return RSparseMatrix(this->sparseMapMatrix());
}
template <> CSparseMatrix 
BlockMatrix< Complex >::sparseMatrix() const {
    log(Warning, "Efficiency check needed.");
    return CSparseMatrix(this->sparseMapMatrix());
}



} // namespace GIMLI

