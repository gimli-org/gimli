/******************************************************************************
 *   Copyright (C) 2007-2020 by the GIMLi development team                    *
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

#if OPENBLAS_CBLAS_FOUND
    #include <cblas.h>
#endif


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
        throwLengthError(WHERE_AM_I + " " + str(cols) + " != " + str(b.size()));
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
        throwLengthError(WHERE_AM_I + " " + str(cols) + " < " + str(endI) + "-" + str(startI));
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
        throwLengthError(WHERE_AM_I + " " + str(rows) + " != " + str(b.size()));
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

template < class ValueType > Matrix < ValueType > &
_transAdd(Matrix < ValueType > * a, const Matrix < ValueType > & b){
    if (a->rows() != b.cols() || a->cols() != b.rows()){
        __MS(a->rows() << " " << b.cols() << " " << a->cols() << " " << b.rows())
        log(Error, "Matrix trans add with wrong dimensions");
        return *a;
    }

    for (Index i = 0; i < a->rows(); i ++ ){
        for (Index j = 0; j < a->cols(); j ++ ){
            a->mat_[i][j] += b.mat_[j][i];
        }
    }
    return *a;
}
template <> Matrix < double > &
Matrix<double>::transAdd(const Matrix < double > & a){
    return _transAdd(this, a);
}
template <> Matrix < Complex > &
Matrix<Complex>::transAdd(const Matrix < Complex > & a){
    return _transAdd(this, a);
}


void matMultABA(const RMatrix & A, const RMatrix & B, RMatrix & C,
                RMatrix & AtB, double a, double b){
    // C = a A.T * B * A + b * C
    // __MS("matMultABA: "<< A.rows() << " " << A.cols() << " : " << B.rows() << " " << B.cols())

    if (A.rows() != B.rows()){
        log(Error, "matMultABA B sizes mismatch.", A.rows(), "!=", B.rows());
        return;
    }
    AtB.resize(A.cols(), B.rows());
    matTransMult(A, B, AtB, 1.0, 0);
    matMult(AtB, A, C, a, b);
}

void matMult(const RMatrix & A, const RMatrix & B, RMatrix & C, double a, double b){
    // C = a * A*B + b *C || C += a * A*B.T + b*C
    // __MS("matMult: "<< A.rows() << " " << A.cols() << " : " << B.rows() << " " << B.cols())
    Index m = A.rows(); // C.rows()
    Index n = B.cols(); // C.cols()
    
    if (A.cols() == B.rows()){ // A * B (k == k)
        C.resize(m, n);

#if OPENBLAS_CBLAS_FOUND
        Index k = A.cols(); // B.rows()

        double *A2 = new double[m * k];
        double *B2 = new double[k * n];
        double *C2 = new double[m * n];

        A.dumpData(A2);
        B.dumpData(B2);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, 
                    a, A2, m, B2, k, b, C2, m);
        
        C.fromData(C2, m, n);

        delete [] A2;
        delete [] B2;
        delete [] C2;
#else
    // __MS("matMult: "<< A.rows() << " " << A.cols() << " : " << B.rows() << " " << B.cols())

        for (Index i = 0; i < A.rows(); i ++){
            for (Index j = 0; j < B.cols(); j ++){
                double c = 0;
                for (Index k = 0; k < A.cols(); k ++){
                    c += A[i][k] * B[k][j];
                }
                C[i][j] = b*C[i][j] + a * c;
            }
        }
#endif
    } else { // A * B.T
        log(Error, "matMult sizes mismatch. implement fallback A*.B.T", A.cols(), "!=", B.rows());
    }
}

void matTransMult(const RMatrix & A, const RMatrix & B, RMatrix & C, double a, double b){
    //** C = a * A.T*B + b*C|| C = a * A.T*B.T  + b*C
    //** C = (a * A.T*B).T + b*C if C has the right dimension
    //** implement with openblas dgemm too and check performance
    // __MS("matTransMult: "<< A.rows() << " " << A.cols() << " : " << B.rows() << " " << B.cols())
    bool retTrans = false;

    // A(k, m).T * B(k, n) = C(m, n)

    Index k = A.rows(); // B.rows()
    Index m = A.cols(); // C.rows()
    Index n = B.rows(); // C.cols()

    if (A.rows() == B.rows()){ // A.T * B

        if (C.rows() != A.cols() || C.cols() != B.cols()){

            // __MS(C.rows() << " " << C.cols() << " " << A.cols() << " " << B.cols())
            //** Target array have wrong dimensions
            if (C.rows() == B.cols() && C.cols() == A.cols()){
                //** Target array seems needed to be transposed
                //** C += (a * A.T*B).T
                __MS("ret transmult")
                retTrans = true;
            } else {
                //** resize target array
                C.resize(A.cols(), B.cols());
            }
        }

#if OPENBLAS_CBLAS_FOUND_

        double *A2 = new double[k * m];
        double *B2 = new double[k * n];
        double *C2 = new double[m * n];

        A.dumpData(A2);
        B.dumpData(B2);
        //specifies row-major (C) or column-major (Fortran) data ordering.
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, k, n, m, 
                    a, A2, k, B2, k, b, C2, m);
            
        C.fromData(C2, m, n);

        delete [] A2;
        delete [] B2;
        delete [] C2;

    // cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, n, n, m, 1, A3, m, B3, m, 0, C3, n);

#else

        for (Index i = 0; i < A.cols(); i ++){
            for (Index j = 0; j < B.cols(); j ++){
                double c = 0;
                for (Index k = 0; k < A.rows(); k ++){
                    c += A[k][i] * B[k][j];
                }
                if (retTrans){
                    C[j][i] = a * c + C[j][i]*b;
                } else {
                    C[i][j] = a * c + C[j][i]*b;
                }
            }
        }
#endif
    } else { // A.T * B.T
        __MS(A)
        __MS(B)
        log(Error, "matTransMult sizes mismatch.", A.rows(), "!=", B.rows());
    }
}



} // namespace GIMLI{
