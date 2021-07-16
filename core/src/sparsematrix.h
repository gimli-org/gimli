/******************************************************************************
 *   Copyright (C) 2006-2021 by the GIMLi development team                    *
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

#ifndef GIMLI_SPARSEMATRIX__H
#define GIMLI_SPARSEMATRIX__H

#include "gimli.h"
#include "elementmatrix.h"
#include "vector.h"
#include "vectortemplates.h"
#include "matrix.h"
#include "mesh.h"
#include "meshentities.h"
#include "node.h"
#include "stopwatch.h"
#include "sparsemapmatrix.h"

#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <cmath>

namespace GIMLI{


#define SPARSE_NOT_VALID throwError(WHERE_AM_I + " no data/or sparsity pattern defined.");

//! Sparse matrix in compressed row storage (CRS) form
/*! Sparse matrix in compressed row storage (CRS) form.
* IF you need native CCS format you need to transpose CRS
* Symmetry type: 0 = nonsymmetric, -1 symmetric lower part, 1 symmetric upper part.*/

template < class ValueType > class SparseMatrix : public SparseMatrixBase{
public:

  /*! Default constructor. Builds invalid sparse matrix */
    SparseMatrix()
        : SparseMatrixBase(), valid_(false), stype_(0), rows_(0), cols_(0){ }

    /*! Copy constructor. */
    SparseMatrix(const SparseMatrix < ValueType > & S)
        : SparseMatrixBase(),
          colPtr_(S.vecColPtr()),
          rowIdx_(S.vecRowIdx()),
          vals_(S.vecVals()), valid_(true), stype_(0){
          rows_ = S.rows();
          cols_ = S.cols();
    }

    /*! Copy constructor. */
    #ifndef PYGIMLI_CAST // disallow automatic type conversion in python
    SparseMatrix(const SparseMapMatrix< ValueType, Index > & S)
        : SparseMatrixBase(), valid_(true){
        copy_(S);
    }
    #endif
    SparseMatrix(const IndexArray & colPtr,
                 const IndexArray & rowIdx,
                 const Vector < ValueType > vals, int stype=0)
        : SparseMatrixBase(){
        colPtr_ = std::vector < int >(colPtr.size());
        rowIdx_ = std::vector < int >(rowIdx.size());
        for (Index i = 0; i < colPtr_.size(); i ++ ) colPtr_[i] = colPtr[i];
        for (Index i = 0; i < rowIdx_.size(); i ++ ) rowIdx_[i] = rowIdx[i];
        vals_   = vals;
        stype_  = stype;
        valid_  = true;
        cols_ = max(rowIdx_) + 1;
        rows_ = colPtr_.size() - 1;
    }

    SparseMatrix(const std::vector < int > & colPtr,
                 const std::vector < int > & rowIdx,
                 const Vector < ValueType > vals, int stype=0)
        : SparseMatrixBase(){
          colPtr_ = colPtr;
          rowIdx_ = rowIdx;
          vals_   = vals;
          stype_  = stype;
          valid_  = true;
          cols_ = max(rowIdx_) + 1;
          rows_ = colPtr_.size() - 1;
    }

    /*! Destructor */
    virtual ~SparseMatrix(){ }

    /*! Copy assignment operator. */
    SparseMatrix < ValueType > & operator = (const SparseMatrix < ValueType > & S){
        if (this != &S){
            colPtr_ = S.vecColPtr();
            rowIdx_ = S.vecRowIdx();
            vals_   = S.vecVals();
            stype_  = S.stype();
            valid_  = true;
            cols_ = S.cols();
            rows_ = S.rows();

        } return *this;
    }

    SparseMatrix < ValueType > & operator = (const SparseMapMatrix< ValueType, Index > & S){
        this->copy_(S);
        return *this;
    }

    #define DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__(OP, FUNCT) \
        void FUNCT(int i, int j, ValueType val){ \
            if ((stype_ < 0 && i > j) || (stype_ > 0 && i < j)) return; \
            if (abs(val) > TOLERANCE || 1){ \
                for (int k = colPtr_[i]; k < colPtr_[i + 1]; k ++){ \
                    if (rowIdx_[k] == j) { \
                        vals_[k] OP##= val; return; \
                    } \
                } \
                std::cerr << WHERE_AM_I << " pos " << i << " " << j << " is not part of the sparsity pattern " << std::endl; \
            } \
        }\
        SparseMatrix< ValueType > & operator OP##= (const ValueType & v){\
            vals_ OP##= v; \
            return *this;\
        }\

        DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__(+, addVal)
        DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__(-, subVal)
        DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__(*, mulVal)
        DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__(/, divVal)

    #undef DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__

    SparseMatrix< ValueType > & operator += (const SparseMatrix< ValueType > & A){
        vals_ += A.vecVals();
        return *this;
    }
    SparseMatrix< ValueType > & operator -= (const SparseMatrix< ValueType > & A){
        vals_ -= A.vecVals();
        return *this;
    }

    SparseMatrix< ValueType > & operator += (const ElementMatrix< double > & A){
        if (!valid_) SPARSE_NOT_VALID;
        for (Index i = 0, imax = A.size(); i < imax; i++){
            for (Index j = 0, jmax = A.size(); j < jmax; j++){
                addVal(A.idx(i), A.idx(j), A.getVal(i, j));
            }
        }
        return *this;
    }

    virtual uint rtti() const { return GIMLI_SPARSE_CRS_MATRIX_RTTI; }

    /*! Return this * a  */
    virtual Vector < ValueType > mult(const Vector < ValueType > & a) const {
        if (a.size() < this->cols()){
            throwLengthError(WHERE_AM_I + " SparseMatrix size(): " + str(this->cols()) + " a.size(): " +
                                str(a.size())) ;
        }

        Vector < ValueType > ret(this->rows(), 0.0);

        if (stype_ == 0){
            // for each row
            for (Index i = 0; i < this->rows(); i++){
            // iterate through compressed col
                for (int j = this->vecColPtr()[i]; j < this->vecColPtr()[i + 1]; j ++){
                    ret[i] += a[this->vecRowIdx()[j]] * this->vecVals()[j];
                }
            }
        } else if (stype_ == -1){
            Index J;
            for (Index i = 0; i < ret.size(); i++){
                for (int j = this->vecColPtr()[i]; j < this->vecColPtr()[i + 1]; j ++){
                    J = this->vecRowIdx()[j];

//                     __MS( i << "  " << J << " " << this->vecVals()[j])
                    ret[i] += a[J] * conj(this->vecVals()[j]);

                    if (J > i){
//                         __MS( J << "  " << i << " " << this->vecVals()[j])
                        ret[J] += a[i] * this->vecVals()[j];
                    }
                }
            }

            //#THROW_TO_IMPL
        } else if (stype_ == 1){
            Index J;
            for (Index i = 0; i < ret.size(); i++){
                for (int j = this->vecColPtr()[i]; j < this->vecColPtr()[i + 1]; j ++){
                    J = this->vecRowIdx()[j];

//                     __MS( i << "  " << J << " " << this->vecVals()[j])
                    ret[i] += a[J] * conj(this->vecVals()[j]);

                    if (J < i){
//                         __MS( J << "  " << i << " " << this->vecVals()[j])
                        ret[J] += a[i] * this->vecVals()[j];
                    }
                }
            }
        }
        return ret;
    }

    /*! Return this.T * a */
    virtual Vector < ValueType > transMult(const Vector < ValueType > & a) const {

        if (a.size() < this->rows()){
            throwLengthError(WHERE_AM_I + " SparseMatrix size(): " + str(this->rows()) + " a.size(): " +
                                str(a.size())) ;
        }

        Vector < ValueType > ret(this->cols(), 0.0);

        if (stype_ == 0){
            for (Index i = 0; i < this->rows(); i++){
                for (int j = this->vecColPtr()[i]; j < this->vecColPtr()[i + 1]; j ++){
                    ret[this->vecRowIdx()[j]] += a[i] * this->vecVals()[j];
                }
            }

        } else if (stype_ == -1){
            THROW_TO_IMPL
        } else if (stype_ ==  1){
            THROW_TO_IMPL
        }
        return ret;
    }

    void add(const ElementMatrix< double > & A){
        return add(A, ValueType(1.0));
    }
    void add(const ElementMatrix< double > & A,
             ValueType scale);
    void add(const ElementMatrix< double > & A,
             const Pos & scale);
    void add(const ElementMatrix< double > & A,
             const Matrix < ValueType > & scale);

    void clean(){ for (Index i = 0, imax = nVals(); i < imax; i++) vals_[i] = (ValueType)(0); }

    void clear(){
        colPtr_.clear();
        rowIdx_.clear();
        vals_.clear();
        valid_ = false;
        cols_ = 0;
        rows_ = 0;
    }

    void setVal(int i, int j, ValueType val){
        if (abs(val) > TOLERANCE || 1){
            for (int k = colPtr_[i]; k < colPtr_[i + 1]; k ++){
                if (rowIdx_[k] == j) {
                    vals_[k] = val; return;
                }
            }
            std::cerr << WHERE_AM_I << " pos " << i << " " << j << " is not part of the sparsity pattern " << std::endl;
        }
    }

    /*!Get matrix value at i,j. If i and j is not part of the matrix
     * sparsity pattern return 0 and print a warning.
     * This warning can be disabled by setting warn to false.*/
    ValueType getVal(int i, int j, bool warn=true) const {
        for (int k = colPtr_[i]; k < colPtr_[i + 1]; k ++){
            if (rowIdx_[k] == j) {
                return vals_[k];
            }
        }
        if (warn) std::cerr << WHERE_AM_I << " pos " << i << " "
                    << j << " is not part of the sparsity pattern " << std::endl;
        return ValueType(0);
    }

    void cleanRow(int row){
        ASSERT_RANGE(row, 0, (int)this->rows())
        for (int col = colPtr_[row]; col < colPtr_[row + 1]; col ++){
            vals_[col] = ValueType(0);
        }
    }

    void cleanCol(int col){
        ASSERT_RANGE(col, 0, (int)this->cols())
        for (int i = 0; i < (int)this->rowIdx_.size(); i++){
            if (rowIdx_[i] == col) {
                vals_[i] = ValueType(0);
            }
        }
    }
    void copy_(const SparseMapMatrix< double, Index > & S);
    void copy_(const SparseMapMatrix< Complex, Index > & S);

    void buildSparsityPattern(const Mesh & mesh){
        Stopwatch swatch(true);

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

        int nVals = 0;
        for (std::vector < std::set< Index > >::iterator mIt = idxMap.begin();
             mIt != idxMap.end(); mIt++){
            //std::sort((*mIt).begin(), (*mIt).end());
            nVals += (*mIt).size();
        }

//         std::cout << "timwe: " << swatch.duration( true) << std::endl;
//         exit(0);

        rowIdx_.reserve(nVals);
        rowIdx_.resize(nVals);
        vals_.resize(nVals);

        colPtr_[0] = 0;
        Index k = 0;
        row = 0;
        for (std::vector < std::set< Index > >::iterator mIt = idxMap.begin(); mIt != idxMap.end(); mIt++){
            for (std::set< Index >::iterator sIt = (*mIt).begin(); sIt != (*mIt).end(); sIt++){
                rowIdx_[k] = (*sIt);
                vals_[k] = (ValueType)0.0;
                k++;
            }
            row++;
            colPtr_[row] = k;
        }
        valid_ = true;

        rows_ = colPtr_.size() - 1;
        cols_ = max(rowIdx_) + 1;
        //** freeing idxMap ist expensive
    }

    void fillStiffnessMatrix(const Mesh & mesh){
        RVector a(mesh.cellCount(), 1.0);
        fillStiffnessMatrix(mesh, a);
    }

    void fillStiffnessMatrix(const Mesh & mesh, const RVector & a){
        clean();
        buildSparsityPattern(mesh);
        ElementMatrix < double > A_l;

        for (uint i = 0; i < mesh.cellCount(); i ++){
            A_l.ux2uy2uz2(mesh.cell(i));
            A_l *= a[mesh.cell(i).id()];
            *this += A_l;
        }
    }
    void fillMassMatrix(const Mesh & mesh){
        RVector a(mesh.cellCount(), 1.0);
        fillMassMatrix(mesh, a);
    }

    void fillMassMatrix(const Mesh & mesh, const RVector & a){
        clean();
        buildSparsityPattern(mesh);
        ElementMatrix < double > A_l;

        for (uint i = 0; i < mesh.cellCount(); i ++){
            A_l.u2(mesh.cell(i));
            A_l *= a[mesh.cell(i).id()];
            *this += A_l;
        }
    }

    /*! symmetric type. 0 = nonsymmetric, -1 symmetric lower part, 1 symmetric upper part.*/
    inline int stype() const {return stype_;}

    inline int * colPtr() { if (valid_) return &colPtr_[0]; else SPARSE_NOT_VALID;  return 0; }
    inline const int & colPtr() const { if (valid_) return colPtr_[0]; else SPARSE_NOT_VALID; return colPtr_[0]; }
    inline const std::vector < int > & vecColPtr() const { return colPtr_; }

    inline int * rowIdx() { if (valid_) return &rowIdx_[0]; else SPARSE_NOT_VALID; return 0; }
    inline const int & rowIdx() const { if (valid_) return rowIdx_[0]; else SPARSE_NOT_VALID; return rowIdx_[0]; }
    inline const std::vector < int > & vecRowIdx() const { return rowIdx_; }

    inline ValueType * vals() { if (valid_) return &vals_[0]; else SPARSE_NOT_VALID; return 0; }
//     inline const ValueType * vals() const { if (valid_) return &vals_[0]; else SPARSE_NOT_VALID; return 0; }
//     inline const ValueType & vals() const { if (valid_) return vals_[0]; else SPARSE_NOT_VALID; return vals_[0]; }
    inline const Vector < ValueType > & vecVals() const { return vals_; }
    inline Vector < ValueType > & vecVals() { return vals_; }

    inline Index size() const { return rows(); }
    inline Index nVals() const { return vals_.size(); }
    inline Index cols() const { return cols_; }
    inline Index rows() const { return rows_; }
    inline Index nCols() const { return cols(); }
    inline Index nRows() const { return rows(); }

    void save(const std::string & fileName) const {
        if (!valid_) SPARSE_NOT_VALID;
        std::fstream file;
        openOutFile(fileName, & file);

        file.setf(std::ios::scientific, std::ios::floatfield);
        file.precision(14);

        for (Index i = 0; i < size(); i++){
            for (SIndex j = colPtr_[i]; j < colPtr_[i + 1]; j ++){
                file << i << "\t" << rowIdx_[j]
                          << "\t" << vals_[j] << std::endl;
            }
        }
        file.close();
    }

    bool valid() const { return valid_; }

protected:

    // int to be cholmod compatible!!!!!!!!

    std::vector < int > colPtr_;
    std::vector < int > rowIdx_;
    Vector < ValueType > vals_;

    bool valid_;
    int stype_;
    Index rows_;
    Index cols_;
};

template < class ValueType >
SparseMatrix< ValueType > operator + (const SparseMatrix< ValueType > & A,
                                      const SparseMatrix< ValueType > & B){
    SparseMatrix< ValueType > ret(A);
    return ret += B;
}

template < class ValueType >
SparseMatrix< ValueType > operator - (const SparseMatrix< ValueType > & A,
                                      const SparseMatrix< ValueType > & B){
    SparseMatrix< ValueType > ret(A);
    return ret -= B;
}

template < class ValueType >
SparseMatrix < ValueType > operator * (const SparseMatrix < ValueType > & A,
                                       const ValueType & b){
    SparseMatrix< ValueType > ret(A);
    return ret *= b;
}

template < class ValueType >
SparseMatrix < ValueType > operator * (const ValueType & b,
                                       const SparseMatrix < ValueType > & A){
    SparseMatrix< ValueType > ret(A);
    return ret *= b;
}

inline RVector operator * (const RSparseMatrix & A, const RVector & b){return A.mult(b);}
inline RVector transMult(const RSparseMatrix & A, const RVector & b){return A.transMult(b);}

inline CVector operator * (const CSparseMatrix & A, const CVector & b){return A.mult(b);}
inline CVector operator * (const CSparseMatrix & A, const RVector & b){return A.mult(toComplex(b));}
inline CVector transMult(const CSparseMatrix & A, const CVector & b){return A.transMult(b);}
inline CVector transMult(const CSparseMatrix & A, const RVector & b){return A.transMult(toComplex(b));}

inline CSparseMatrix operator + (const CSparseMatrix & A, const RSparseMatrix & B){
    CSparseMatrix ret(A);
    ret.vecVals() += toComplex(B.vecVals());
    return ret;
}

inline RSparseMatrix real(const CSparseMatrix & A){
    return RSparseMatrix(A.vecColPtr(), A.vecRowIdx(),
                         real(A.vecVals()), A.stype());
}
inline RSparseMatrix imag(const CSparseMatrix & A){
    return RSparseMatrix(A.vecColPtr(), A.vecRowIdx(),
                         imag(A.vecVals()), A.stype());
}

/*! SparseMatrix specialized type traits in sparsematrix.cpp */
template< typename ValueType >
void SparseMatrix< ValueType >::copy_(const SparseMapMatrix< double, Index > & S){THROW_TO_IMPL}
template< typename ValueType >
void SparseMatrix< ValueType >::copy_(const SparseMapMatrix< Complex, Index > & S){THROW_TO_IMPL}


template <> DLLEXPORT void SparseMatrix<double>::copy_(const SparseMapMatrix< double, Index > & S);
template <> DLLEXPORT void SparseMatrix<Complex>::copy_(const SparseMapMatrix< Complex, Index > & S);

template <> DLLEXPORT void SparseMatrix< double >::
    add(const ElementMatrix < double > & A, double scale);
template <> DLLEXPORT void SparseMatrix< double >::
    add(const ElementMatrix < double > & A, const Pos & scale);
template <> DLLEXPORT void SparseMatrix< double >::
    add(const ElementMatrix < double > & A, const RMatrix & scale);

template <> DLLEXPORT void SparseMatrix< Complex >::
    add(const ElementMatrix < double > & A, Complex scale);
template <> DLLEXPORT void SparseMatrix< Complex >::
    add(const ElementMatrix < double > & A, const Pos & scale);
template <> DLLEXPORT void SparseMatrix< Complex >::
    add(const ElementMatrix < double > & A, const CMatrix & scale);


} // namespace GIMLI

#endif //GIMLI_SPARSEMATRIX__H
