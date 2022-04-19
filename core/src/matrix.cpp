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
#include "vector.h"
#include "elementmatrix.h"

#include "matrix.h"
// #include "vector.h"
#include "stopwatch.h"

#if OPENBLAS_CBLAS_FOUND
    #include <cblas.h>
#endif


namespace GIMLI{

template class DenseMatrix< double >;

#define _MATRIX_WS_SIZE 8192

static double _wsA[_MATRIX_WS_SIZE];
static double _wsB[_MATRIX_WS_SIZE];
static double _wsC[_MATRIX_WS_SIZE];


void toEigenMatrix(const RMatrix & m, SmallMatrix & r){
    r.resize(m.rows(), m.cols());

    for (Index i=0; i < r.rows(); i ++){
#if USE_EIGEN3
        r(i, Eigen::all) = Eigen::Map <const Eigen::VectorXd>(&m[i][0],
                                                              r.cols());
#else
    for (Index j=0; j < r.cols(); j ++){
        r(i,j) = m(i,j);
    }
#endif
    }
}
void toRMatrix(const SmallMatrix & m, RMatrix & r){
    r.resize(m.rows(), m.cols());

    // optimize if in use
    for (Index i=0; i < r.rows(); i ++){
        for (Index j=0; j < r.cols(); j ++){
            r(i,j) = m(i,j);
        }
    }
}

#if USE_EIGEN3
void toRVector(const Eigen::VectorXd & m, RVector & r, double b){
    r.resize(m.size());

    // optimize if in use
    if (b == 0.0){
        for (Index i=0; i < m.size(); i ++){
            r[i] = m(i);
        }
    } else if (b == 1.0){
        for (Index i=0; i < m.size(); i ++){
            r[i] += m(i);
        }
    } else if (b == -1.0){
        for (Index i=0; i < m.size(); i ++){
            r[i] -= m(i);
        }
    } else if (b == -1.0){
        for (Index i=0; i < m.size(); i ++){
            r[i] = m(i) + r[i] * b;
        }
    }
}
#endif

template <> Vector< double >
DenseMatrix< double >::row(Index i) const {
    ASSERT_THIS_SIZE(i)
    return Vector< double >(this->_cols, _data, this->_cols * i);
}
template <> Vector< double >
DenseMatrix< double >::row(Index i) {
    ASSERT_THIS_SIZE(i)
    return Vector< double >(this->_cols, _data, this->_cols * i);
}
template <> void
DenseMatrix< double >::round(const double & tol){
    Vector < double > view(length(), this->_data, 0);
    view.round(tol);
}
template <> Vector< Complex >
DenseMatrix< Complex >::row(Index i) const {
    ASSERT_THIS_SIZE(i)
    return Vector< Complex >(this->_cols, _data, this->_cols * i);
}
template <> Vector< Complex >
DenseMatrix< Complex >::row(Index i) {
    ASSERT_THIS_SIZE(i)
    return Vector< Complex >(this->_cols, _data, this->_cols * i);
}
template <> void
DenseMatrix< Complex >::round(const Complex & tol){
    THROW_TO_IMPL
}

template < class ValueType, class Mat >
void mult_T(const Mat & A, const Vector < ValueType > & b,
          Vector < ValueType > & c,
          const ValueType & alpha, const ValueType & beta){
    Stopwatch sw;
    for (Index i = 0; i < A.rows(); i ++){
        if (beta == 0.0){
            c[i] = alpha * sum(A.row(i) * b);
        } else {
            c[i] = alpha * sum(A.row(i) * b) + beta * c[i];
        }
    }
    print("mult impl:", sw.duration());
}

void mult(const RDenseMatrix & A, const Vector < double > & b,
          Vector < double > & c,
          const double & alpha, const double & beta){
    Index m = A.cols();
    Index n = A.rows();
    ASSERT_VEC_SIZE(b, n)
    c.resize(m);

#if OPENBLAS_CBLAS_FOUND
    Stopwatch sw;
    if (noCBlas()){
        mult_T(A, b, c, alpha, beta);
    } else {
        cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, alpha,
                    A.pData(), n, &b[0],
                    1, beta, &c[0], 1);

        print("dgemv:", sw.duration());
    }
#else //#if OPENBLAS_CBLAS_FOUND
    mult_T(A, b, c, alpha, beta);
#endif
}
void mult(const CDenseMatrix & A, const Vector < Complex > & b,
          Vector < Complex > & c,
          const Complex & alpha, const Complex & beta){
    Index m = A.cols();
    Index n = A.rows();
    ASSERT_VEC_SIZE(b, n)
    c.resize(m);

#if OPENBLAS_CBLAS_FOUND
    Stopwatch sw;
    if (noCBlas()){
        mult_T(A, b, c, alpha, beta);
    } else {
        THROW_TO_IMPL
        // cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, alpha,
        //             A.pData(), n, &b[0],
        //             1, beta, &c[0], 1);

        // print("dgemv:", sw.duration());
    }
#else //#if OPENBLAS_CBLAS_FOUND
    mult_T(A, b, c, alpha, beta);
#endif
}

template <> DLLEXPORT RVector
DenseMatrix<double>::mult(const RVector & b, Index startI, Index endI) const{
    THROW_TO_IMPL
    return RVector();
}
template <> DLLEXPORT CVector
DenseMatrix<Complex>::mult(const CVector & b, Index startI, Index endI) const{
    THROW_TO_IMPL
    return CVector();
}


template <> Vector<double>
DenseMatrix<double>::transMult(const Vector < double > & b) const{
    ASSERT_VEC_SIZE(b, this->_rows)
    Vector < double > ret(this->_cols);

    for (Index i = 0; i < this->_cols; i ++){
        ret[i] = GIMLI::dot(col(i), b);
    }
    return ret;
}
template <> Vector<Complex>
DenseMatrix<Complex>::transMult(const Vector < Complex > & b) const{
    ASSERT_VEC_SIZE(b, this->_rows)
    Vector < double > ret(this->_cols);
    THROW_TO_IMPL
    return ret;
}
template <> DenseMatrix<double> &
DenseMatrix<double>::transAdd(const DenseMatrix < double > & a) {
    THROW_TO_IMPL
    return *this;
}
template <> DenseMatrix<Complex> &
DenseMatrix<Complex>::transAdd(const DenseMatrix < Complex > & a){
    THROW_TO_IMPL
    return *this;
}

//###########
// Matrix related implementations
//###########
void mult(const RMatrix & A, const RVector & b, RVector & c,
          const double & alpha, const double & beta){
    ASSERT_VEC_SIZE(b, A.rows())
    c.resize(A.cols());
    return mult_T(A, b, c, alpha, beta);
    // Bufferalloc for OPENBLAS to xpensive
}

void mult(const CMatrix & A, const CVector & b,
          CVector & c, const Complex & alpha, const Complex & beta){
    ASSERT_VEC_SIZE(b, A.rows())
    c.resize(A.cols());
    return mult_T(A, b, c, alpha, beta);
    // Bufferalloc for OPENBLAS to xpensive
}

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
    // a+=b.T
    if (a->rows() != b.cols() || a->cols() != b.rows()){
        __MS(a->rows(), a->cols(), ":",  b.rows(), b.cols())
        log(Error, "Matrix _transAdd with wrong dimensions");
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

template < class Mat >
void matMultABA_T(const Mat & A, const Mat & B,
                   Mat & C, Mat & AtB, double a, double b);
template < class Mat >
void matMult_T(const Mat & A, const Mat & B,
                Mat & C, double a, double b);
template < class Mat >
void matMult_T_impl(const Mat & A, const Mat & B, Mat & C,
                 double a, double b, bool bIsTrans, Index n);

template < class Mat >
void matTransMult_T(const Mat & A, const Mat & B,
                    Mat & C, double a, double b);
template < class Mat >
void matTransMult_T_impl(const Mat & A, const Mat & B,
                        Mat & C, double a, double b, bool bIsTrans, Index n);


void matMultABA(const RDenseMatrix & A, const RDenseMatrix & B,
                RDenseMatrix & C, RDenseMatrix & AtB, double a, double b){
    THROW_TO_IMPL
}
void matMultABA(const SmallMatrix & A, const SmallMatrix & B,
                SmallMatrix & C, SmallMatrix & AtB, double a, double b){

#if USE_EIGEN3
    THROW_TO_IMPL
#else
    return matMultABA_T(A, B, C, AtB, a, b);
#endif
}

template < class Mat >
void matMultABA_T(const Mat & A, const Mat & B,
                   Mat & C, Mat & AtB, double a, double b){
    // C = a A.T * B * A + b * C
    // __MS("matMultABA: ", A.rows(), A.cols(), B.rows(), B.cols())

    if (A.rows() != B.rows()){
        log(Error, "matMultABA B sizes mismatch.", A.rows(), "!=", B.rows());
        return;
    }
    AtB.resize(A.cols(), B.rows());
    matTransMult_T(A, B, AtB, 1.0, 0.0);
    matMult_T(AtB, A, C, a, b);
}

template < class Mat >
void matMult_T_impl(const Mat & A, const Mat & B, Mat & C,
                 double a, double b, bool bIsTrans, Index n){

    Stopwatch sw;

    for (Index i = 0; i < A.rows(); i ++){
        for (Index j = 0; j < n; j ++){
            double c = 0;
            for (Index k = 0; k < A.cols(); k ++){
                if (bIsTrans){
                    c += A[i][k] * B[j][k];
                } else {
                    c += A[i][k] * B[k][j];
                }
            }
            if (b == 0.0){
                C[i][j] = a * c;
            } else if (b == 1.0){
                C[i][j] += a * c;
            } else if (b == -1.0){
                C[i][j] -= a * c;
            } else {
                C[i][j] = b * C[i][j] + a * c;
            }
        }
    }

    print("matmult_impl:", sw.duration());

}
template < class Mat >
void matMult_T(const Mat & A, const Mat & B,
                Mat & C, double a, double b){
    // C = a * A*B + b*C || C = a * A*B.T + b*C

    // __MS("matMult: ", A.rows(), A.cols(), B.rows(), B.cols())

    Index m = A.rows(); // C.rows()
    Index k = A.cols();

    // B.rows()
    Index n = B.cols(); // C.cols()
    bool bIsTrans = false;
    Index bRows = n;

    if (k == B.rows()){ // A * B (k == k)

    } else if (k == B.cols()){ // A * B.T (k == k)
        bIsTrans = true;
        n = B.rows();
        bRows = k;
    } else {
        log(Error, "matMult sizes mismatch. ",
            A.cols(), "!=", B.rows());
    }
    C.resize(m, n);

#if OPENBLAS_CBLAS_FOUND

    if (noCBlas()){
        return matMult_T_impl(A, B, C, a, b, bIsTrans, n);
    }

    CBLAS_TRANSPOSE aTrans = CblasNoTrans;
    CBLAS_TRANSPOSE bTrans = CblasNoTrans;

    if (bIsTrans) bTrans = CblasTrans;

// not threadsafe at all
    if ((k * m) > _MATRIX_WS_SIZE || (k * n) > _MATRIX_WS_SIZE || (m * n) > _MATRIX_WS_SIZE){
        __MS(k * m, k * n, m * n )
        THROW_TO_IMPL
    }

    double * bA = A.toData(&_wsA[0], _MATRIX_WS_SIZE);
    double * bB = B.toData(&_wsB[0], _MATRIX_WS_SIZE);
    double * bC = C.toData(&_wsC[0], _MATRIX_WS_SIZE);
    // __MS(bC)
    // A.dumpData(&_wsA[0]);
    // B.dumpData(&_wsB[0]);
    // C.dumpData(&_wsC[0]);

    // double *A2 = new double[m * k];
    // double *B2 = new double[k * n];
    // double *C2 = new double[m * n];

    // A.dumpData(A2);
    // B.dumpData(B2);
    // C.dumpData(C2);

    // cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k,
    //             a, A2, m, B2, n, b, C2, n);

    // lda ## leading dimension for a, means column for CblasRowMajor
    Stopwatch sw;
    cblas_dgemm(CblasRowMajor, aTrans, bTrans,
                m, n, k,
                a, bA, k, bB, bRows,
                b, bC, n);

    // if (debug()){
        print("dgemm:", sw.duration());
    // }
    // __MS(bC)
    C.fromData(bC, m, n);


    // delete [] A2;
    // delete [] B2;
    // delete [] C2;
#else

    matMult_T_impl(A, B, C, a, b, bIsTrans, n);
    // __MS("\t: ", C.rows(), C.cols(), bIsTrans)

#endif
}

void matMult(const RDenseMatrix & A, const RDenseMatrix & B,
             RDenseMatrix & C, double a, double b){
    return matMult_T(A, B, C, a, b);
}

void matMult(const SmallMatrix & A, const SmallMatrix & B,
             SmallMatrix & C, double a, double b){
#if USE_EIGEN3
    if (A.cols() == B.rows()){

        if (C.rows() != A.rows() || C.cols() != B.cols()){
            C.resize(A.rows(), B.cols());
        }
        if (b == 0.0){
            C = a*(A*B);
        } else if (b == 1.0){
            C += a*(A*B);
        } else if (b == -1.0){
            C += a*(A*B);
        } else {
            C = b*C + a*(A*B);
        }
    } else if (A.cols() == B.cols()){
        if (C.rows() != A.rows() || C.cols() != B.rows()){
            C.resize(A.rows(), B.rows());
        }
        if (b == 0.0){
            C = a*(A*B.transpose());
        } else if (b == 1.0){
            C += a*(A*B.transpose());
        } else if (b == -1.0){
            C += a*(A*B.transpose());
        } else {
            C = b*C + a*(A*B.transpose());
        }
    } else {
        log(Error, "matMult sizes mismatch. ", A.cols(), "!=", B.rows());
    }

#else
    return matMult_T(A, B, C, a, b);
#endif
}
void matTransMult(const RDenseMatrix & A, const RDenseMatrix & B,
                  RDenseMatrix & C, double a, double b){
    matTransMult_T(A, B, C, a, b);
}

void matTransMult(const SmallMatrix & A, const SmallMatrix & B,
                  SmallMatrix & C, double a, double b){
//** C = a * A.T*B + b*C || C = a * A.T*B.T + b*C

#if USE_EIGEN3
    if (A.rows() == B.rows()){

        if (C.rows() != A.cols() || C.cols() != B.cols()){
            C.resize(A.cols(), B.cols());
        }
        if (b == 0.0){
            C = a*(A.transpose()*B);
        } else if (b == 1.0){
            C += a*(A.transpose()*B);
        } else if (b == -1.0){
            C += a*(A.transpose()*B);
        } else {
            C = b*C + a*(A.transpose()*B);
        }
    } else if (A.rows() == B.cols()){
        if (C.rows() != A.cols() || C.cols() != B.rows()){
            C.resize(A.cols(), B.rows());
        }
        if (b == 0.0){
            C = a*(A.transpose()*B.transpose());
        } else if (b == 1.0){
            C += a*(A.transpose()*B.transpose());
        } else if (b == -1.0){
            C += a*(A.transpose()*B.transpose());
        } else {
            C = b*C + a*(A.transpose()*B.transpose());
        }
    } else {
        log(Error, "matTransMult sizes mismatch. ", A.rows(), "!=", B.rows());
    }
#else
    matTransMult_T(A, B, C, a, b);
#endif
}

template < class Mat >
void matTransMult_T_impl(const Mat & A, const Mat & B,
                        Mat & C, double a, double b, bool bIsTrans, Index n){
    // private!! only use this after size checks

    for (Index i = 0; i < A.cols(); i ++){
        for (Index j = 0; j < n; j ++){
            double c = 0;

            for (Index k = 0; k < A.rows(); k ++){
                if (bIsTrans){
                    c += A[k][i] * B[j][k];
                } else {
                    c += A[k][i] * B[k][j];
                }
            }
            if (b == 0.0){
                C[i][j] = a * c;
            } else if (b == 1.0){
                C[i][j] += a * c;
            } else {
                C[i][j] = b * C[i][j] + a * c;
            }
        }
    }
}

template < class Mat >
void matTransMult_T(const Mat & A, const Mat & B,
                    Mat & C, double a, double b){

    // __MS("matTransMult: ", A.rows(), A.cols(), B.rows(), B.cols())

    // Mxk * kxN == MxN
    // (kxM).T * kxN == MxN
    // A (k,M), B(k,N) == C(M,N)
    //** C = a * A.T*B + b*C || C = a * A.T*B.T + b*C

    //** **MEH C = (a * A.T*B).T + b*C if C has the right dimension **  MEH
    // c-transpose check only for b != 0.0, else c is resized

    // bool retTrans = false;

    // A(k, m).T * B(k, n) = C(m, n)

    Index k = A.rows(); // B.rows()
    Index m = A.cols(); // C.rows()
    Index n = B.cols(); // C.cols()

    bool bIsTrans = false;
    Index bRows = n;

    if (k == B.rows()){ // A.T * B

        if (b == 0.0) C.resize(m, n);

        if (C.rows() != A.cols() || C.cols() != B.cols()){
            __MS("isinuse?")
            //__MS(C.rows(), C.cols(), A.cols(), B.cols())
            //** Target array have wrong dimensions
            if (C.rows() == B.cols() && C.cols() == A.cols()){
                // C = a * B.T*A + b*C
                //** Target array seems needed to be transposed
                //** C = a*(A.T*B).T + b * C
                // __MS("ret transmult")
                return matTransMult_T(B, A, C, a, b);

                // retTrans = true;
            } else {
                //** resize target array
                C.resize(m, n);
            }
        }
    } else if(k == B.cols()){
        bIsTrans = true;
        n = B.rows();
        bRows = k;

    } else {
        // if not A.T * B or A.T * B.T
        __MS(A.rows(), A.cols(), '\n', A)
        __MS(B.rows(), B.cols(), '\n', B)
        __MS(C.rows(), C.cols())
        throwLengthError("matTransMult sizes mismatch.");
    }

#if OPENBLAS_CBLAS_FOUND

    if (noCBlas()){
        matTransMult_T_impl(A, B, C, a, b, bIsTrans, n);
    }

// __MS("OPENBLAS")
    CBLAS_TRANSPOSE aTrans = CblasTrans;
    CBLAS_TRANSPOSE bTrans = CblasNoTrans;

    if (bIsTrans) bTrans = CblasTrans;

    if ((k * m) > _MATRIX_WS_SIZE || (k * n) > _MATRIX_WS_SIZE || (m * n) > _MATRIX_WS_SIZE){
        __MS(k * m, k * n, m * n )
        THROW_TO_IMPL
    }

    // double *A2 = new double[k * m];
    // double *B2 = new double[k * n];
    // double *C2 = new double[m * n];

    // not threadsafe at all
    double * bA = A.toData(&_wsA[0]);
    double * bB = B.toData(&_wsB[0]);
    double * bC = C.toData(&_wsC[0]);

    // A.dumpData(A2);
    // B.dumpData(B2);
    // C.dumpData(C2);

    // std::cout << "A" << std::endl;
    // for (Index i = 0; i < m*k; i ++ ){std::cout << A2[i] << " ";} std::cout << std::endl;
    // std::cout << "b" << std::endl;
    // for (Index i = 0; i < n*k; i ++ ){std::cout << B2[i] << " ";} std::cout << std::endl;
    // std::cout << "C1" << std::endl;
    // for (Index i = 0; i < m*n; i ++ ){std::cout << C2[i] << " ";} std::cout << std::endl;

    // Stopwatch sw;
    cblas_dgemm(CblasRowMajor, aTrans, bTrans, m, n, k,
                a, bA, m, bB, bRows,
                b, bC, n);
    // if (debug()){
    //     print("dgemm:", sw.duration());
    // }

    // std::cout << "C2" << std::endl;
    // for (Index i = 0; i < m*n; i ++ ){std::cout << C2[i] << " ";} std::cout << std::endl;

    C.fromData(bC, m, n);

    // std::cout << "C3" << std::endl;
    // std::cout << C << std::endl;

    // delete [] A2;
    // delete [] B2;
    // delete [] C2;

#else
    matTransMult_T_impl(A, B, C, a, b, bIsTrans, n);
#endif

}



} // namespace GIMLI{
