/******************************************************************************
 *   Copyright (C) 2007-2019 by the GIMLi development team                    *
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

#include "matrix.h"
// #include "vector.h"
#include "stopwatch.h"

namespace GIMLI{
    
template < class ValueType > Vector < ValueType >
_mult(const Matrix< ValueType > & M, const Vector < ValueType > & b) {
    Index cols = M.cols();
    Index rows = M.rows();

    Vector < ValueType > ret(rows, 0.0);
    
    //ValueType tmpval = 0;
    if (b.size() == cols){
        for (Index i = 0; i < rows; ++i){
            ret[i] = sum(M.mat_[i] * b);
            // for (Index j = 0; j < cols; j++){
            //     ret[i] += M.mat_[i][j] * b[j];
            // }
        }
    } else {
        throwLengthError(1, WHERE_AM_I + " " + str(cols) + " != " + str(b.size()));
    }
    return ret;
}

template<> Vector < double > 
Matrix< double >::mult(const Vector < double > & b) const { return _mult((*this), b); }
template<> Vector < Complex > 
Matrix< Complex >::mult(const Vector < Complex > & b) const { return _mult((*this), b); }

template < class ValueType > Vector < ValueType >
_mult(const Matrix< ValueType > & M, const Vector < ValueType > & b, Index startI, Index endI) {
    Index cols = M.cols();
    Index rows = M.rows();
    Index bsize = Index(endI - startI);

    if (bsize != cols) {
        throwLengthError(1, WHERE_AM_I + " " + str(cols) + " < " + str(endI) + "-" + str(startI));
    }
    Vector < ValueType > ret(rows, 0.0);
    for (Index i = 0; i < rows; ++i){
        for (Index j = startI; j < endI; j++) {
            ret[i] += M.mat_[i][j] * b[j];
        }
    }
    return ret;
}

template<> Vector < double > 
Matrix< double >::mult(const Vector < double > & b, Index startI, Index endI) const {
    return _mult((*this), b, startI, endI);
}
template<> Vector < Complex > 
Matrix< Complex >::mult(const Vector < Complex > & b, Index startI, Index endI) const {
    return _mult((*this), b, startI, endI);
}

template < class ValueType > Vector < ValueType >
_transMult(const Matrix< ValueType > & M, const Vector < ValueType > & b) {
    Index cols = M.cols();
    Index rows = M.rows();
    Vector < ValueType > ret(cols, 0.0);

    if (b.size() == rows){
        for (Index i = 0; i < rows; i++){
            // ret += M.mat_[i] * b[i];
            for (Index j = 0; j < cols; j++){
                ret[j] += M.mat_[i][j] * b[i];
            }
        }
    } else {
        throwLengthError(1, WHERE_AM_I + " " + str(rows) + " != " + str(b.size()));
    }
    return ret;
}

template<> Vector< double > 
Matrix< double >::transMult(const Vector < double > & b) const {
    return _transMult((*this), b);
}
template<> Vector< Complex > 
Matrix< Complex >::transMult(const Vector < Complex > & b) const {
    return _transMult((*this), b);
}    

void matMultABA(const RMatrix & A, const RMatrix & B, RMatrix & C, RMatrix & AtB, double a){
    if (B.rows() != B.cols()){
        log(Error, "matMultABA B sizes mismatch.", B.cols(), "!=", B.rows());
        return;
    }
    AtB.resize(A.cols(), B.rows());
    AtB *= 0.0;
    matTransMult(A, B, AtB, 1.0);
    matMult(AtB, B, C, a);
}

void matMult(const RMatrix & A, const RMatrix & B, RMatrix & C, double a){
    // C += a * A*B || C += a * A*B.T
    // implement with openblas dgemm too and check performance
    if (A.cols() == B.rows()){
        C.resize(A.rows(), B.cols());

        for (Index i = 0; i < A.rows(); i ++){
            for (Index j = 0; j < B.cols(); j ++){
                double c = 0;
                for (Index k = 0; k < A.cols(); k ++){
                    c += A[i][k] * B[k][j];
                }
                C[i][j] += a * c;
            }
        }
    } else {
        log(Error, "matMult sizes mismatch. implement fallback A*.B.T", A.cols(), "!=", B.rows());
    }
}

void matTransMult(const RMatrix & A, const RMatrix & B, RMatrix & C, double a){
    // C += a * A.T*B || C += a * A.T*B.T
    // implement with openblas dgemm too and check performance
    if (A.rows() == B.rows()){
        C.resize(A.cols(), B.cols());

        for (Index i = 0; i < A.cols(); i ++){
            for (Index j = 0; j < B.cols(); j ++){
                double c = 0;
                for (Index k = 0; k < A.rows(); k ++){
                    c += A[k][i] * B[k][j];
                }
                C[i][j] += a * c;
            }
        }
    } else {
        log(Error, "matMult sizes mismatch. implement fallback A.T*B.T", A.rows(), "!=", B.rows());
    }
}



} // namespace GIMLI{
