/***************************************************************************
 *   Copyright (C) 2007-2014 by the resistivity.net development team       *
 *   Thomas Günther  thomas@resistivity.net                                *
 *   Carsten Rücker carsten@resistivity.net                                *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef GIMLI_BLOCKMATRIX__H
#define GIMLI_BLOCKMATRIX__H

#include "gimli.h"
#include "sparsematrix.h"
#include "matrix.h"

//! block matrices for easier inversion, see appendix E in GIMLi tutorial
namespace GIMLI{
    
template < class ValueType > class DLLEXPORT BlockMatrix : public MatrixBase {
public:
    
    struct BlockMatrixEntry {
        Index rowStart;
        Index colStart;
        Index matrixID;
        ValueType scale;
    };
    
    BlockMatrix(bool verbose=false) 
        : MatrixBase(verbose){
        rows_ = 0;
        cols_ = 0;
    }
    
    virtual ~BlockMatrix(){
        this->clear();
    }
    
    /*! Return entity rtti value. */
    virtual uint rtti() const { return GIMLI_BLOCKMATRIX_RTTI; }
    
    virtual Vector < ValueType > operator [] (Index r) const { return row(r); }
    
    virtual Index rows() const { recalcMatrixSize(); return rows_; }
        
    virtual Index cols() const { recalcMatrixSize(); return cols_; }
    
    virtual Vector < ValueType > row(Index r) const {
        Vector < ValueType > b(rows(), 0.0);
        b[r]=1.0;
        return transMult(b);
    }
    
    virtual Vector < ValueType > col(Index r) const{
        Vector < ValueType > b(cols(), 0.0);
        b[r]=1.0;
        return mult(b);
    }
        
    virtual void clear(){
        matrieces_.clear();
        entries_.clear();
    }
    virtual void clean(){
        for (Index i = 0; i < matrieces_.size(); i ++ ){
            matrieces_[i]->clean();
        }
    }
    
    MatrixBase * mat(Index idx) {
        return matrieces_[idx];
    }
    
    MatrixBase & matRef(Index idx) const {
        return *matrieces_[idx];
    }

    Index addMatrix(MatrixBase * matrix){
//         __MS(matrix << " " << matrix->rtti())
        matrieces_.push_back(matrix);
        return matrieces_.size() - 1;
    }

    void addMatrixEntry(Index matrixID, Index rowStart, Index colStart,
                        ValueType scale = ValueType(1.0)){
        if (matrixID > matrieces_.size()){
            throwLengthError(1, WHERE_AM_I + " matrix entry to large: " +
            str(matrixID) + " " + str(matrieces_.size()));
        }
        BlockMatrixEntry entry;
        entry.rowStart = rowStart;
        entry.colStart = colStart;
        entry.matrixID = matrixID;
        entry.scale = scale;
        
        entries_.push_back(entry);
        recalcMatrixSize();
    }
        
    virtual Vector < ValueType > mult(const Vector < ValueType > & b) const{
        if (b.size() != this->cols()){
            throwLengthError(1, WHERE_AM_I + " wrong size of vector b (" +
            str(b.size()) + ") needed: " + str(this->cols()));
            
        }

        Vector < ValueType > ret(rows_);
        
        for (Index i = 0; i < entries_.size(); i ++ ){
            BlockMatrixEntry entry = entries_[i];
            
            MatrixBase &mat = *matrieces_[entry.matrixID];
            
            ret.addVal(mat.mult(b.getVal(entry.colStart, 
                                         entry.colStart + mat.cols())) * 
                       entry.scale,
                       entry.rowStart, entry.rowStart + mat.rows());
        } 
        
        return ret;
    }

    virtual Vector < ValueType > transMult(const Vector < ValueType > & b) const {
        if (b.size() != this->rows()){
            throwLengthError(1, WHERE_AM_I + " wrong size of vector b (" +
            str(b.size()) + ") needed: " + str(this->rows()));
            
        }
        
        Vector < ValueType > ret(cols_);
         for (Index i = 0; i < entries_.size(); i ++ ){
            BlockMatrixEntry entry = entries_[i];
            
            MatrixBase &mat = *matrieces_[entry.matrixID];
            
            ret.addVal(mat.transMult(b.getVal(entry.rowStart, 
                                         entry.rowStart + mat.rows())) * 
                       entry.scale,
                       entry.colStart, entry.colStart + mat.cols());
        } 
        return ret;
    }
    
    void recalcMatrixSize() const {
        
        for (Index i = 0; i < entries_.size(); i ++ ){
            BlockMatrixEntry entry = entries_[i];
//             __MS(matrieces_[entry.matrixID] << " " << matrieces_[entry.matrixID]->rtti())
            rows_ = max(rows_, 
                        entry.rowStart + matrieces_[entry.matrixID]->rows());
            cols_ = max(cols_, 
                        entry.colStart + matrieces_[entry.matrixID]->cols());
        } 
    }
    
    virtual void save(const std::string & filename) const {
        std::cerr << WHERE_AM_I << "WARNING " << " don't save blockmatrix."  << std::endl;
//         THROW_TO_IMPL
    }

protected:
    std::vector< MatrixBase * > matrieces_;
    std::vector< BlockMatrixEntry > entries_;
    
    /*! Max row size.*/ 
    mutable Index rows_;
    
    /*! Max col size.*/ 
    mutable Index cols_;
};


/*! Simple example for tutorial purposes. */
/*! Block Matrix consisting of two horizontally pasted sparse map matrices. */
class DLLEXPORT H2SparseMapMatrix : public MatrixBase{
public:
    
    H2SparseMapMatrix(){}
    
    virtual ~H2SparseMapMatrix(){}
    
    virtual Index rows() const { return H1_.rows(); }
    
    virtual Index cols() const { return H1_.cols() + H2_.cols(); }
    
    virtual void clear() { H1_.clear(); H2_.clear(); }
    
    /*! Return this * a . */
    virtual RVector mult(const RVector & a) const {
        return H1_ * a(0, H1_.cols()) + H2_ * a(H1_.cols(), cols());
    }
    
    /*! Return this.T * a = (a.T * this).T . */
    virtual RVector transMult(const RVector & a) const {
        return cat(H1_.transMult(a), H2_.transMult(a));
    }
    
    /*! Return references to the 2 matriced (const and non-const, why?). */
    inline const RSparseMapMatrix & H1() const { return H1_; }
    inline const RSparseMapMatrix & H2() const { return H2_; }
    inline RSparseMapMatrix & H1() { return H1_; }
    inline RSparseMapMatrix & H2() { return H2_; }

protected:
    //! create inplace (or better hold references of it?)
    RSparseMapMatrix H1_, H2_; 
}; // class H2SparseMapMatrix

inline void rank1Update(H2SparseMapMatrix & A, const RVector & u, const RVector & v) {
    CERR_TO_IMPL
    return;
}

inline bool save(const H2SparseMapMatrix & A, const std::string & filename, IOFormat format = Ascii){
    CERR_TO_IMPL
    return false;
}

// /*! Do we have to do this for every matrix type?? */  
// inline RVector operator * (const H2SparseMapMatrix & A, const RVector & x){
//     return A.H1() * x(0, A.H1().cols()) + A.H2() * x(A.H1().cols(), A.cols());
// }
// 
// inline RVector transMult(const H2SparseMapMatrix & A, const RVector & b){
//     return cat(transMult(A.H1(), b), transMult(A.H2(), b));
// }

/*! Block matrix with 2 arbitrary matrices pasted horizontally. */
template< class Matrix1, class Matrix2 > class H2Matrix{ 
public:
    H2Matrix(){}
    
    virtual ~H2Matrix(){}

    /*! Return rows and columns. */
    virtual Index rows() const { return H1_.rows(); }
    virtual Index cols() const { return H1_.cols() + H2_.cols(); }

    /*! Return this * a . */
    virtual RVector mult(const RVector & a) const {
        return H1_ * a(0, H1_.cols()) + H2_ * a(H1_.cols(), cols());
    }
    
    /*! Return this.T * a = (a.T * this).T . */
    virtual RVector transMult(const RVector & a) const {
        return cat(H1_.transMult(a), H2_.transMult(a));
    }
protected:
    Matrix1 H1_;
    Matrix2 H2_;

}; // H2Matrix

/*! Block matrix with 2 arbitrary matrices pasted vertically. */
template< class Matrix1, class Matrix2 > class V2Matrix : public MatrixBase{ };
/*! Block diagonal matrix with 2 arbitrary matrices as diagonals. */
template< class Matrix1, class Matrix2 > class D2Matrix : public MatrixBase{ };

/*! Block matrix with arbitrary number of matrices of the same type pasted horizontally. */
template< class Matrix > class HNMatrix : public MatrixBase{ 
public:
    HNMatrix(){}
    
    virtual ~HNMatrix(){}
    
    void push_back(Matrix Mat){ //** standard procedure to build up matrix
        if (Mats_.size() > 0 && Mats_[0].nrows() == Mat.rows()) {
            Mats_.push_back(Mat); 
        } else {
            throwLengthError(1, WHERE_AM_I + " matrix rows do not match " +
                                 toStr(Mats_[0].nrows()) + " " + toStr(Mat.nrows()));
        }
    }
    /*! Return rows and columns. */
    virtual Index rows() const { 
        if (Mats_.size() > 0) return Mats_[0].rows(); 
        return 0;
    }
    virtual Index cols() const { 
        Index ncols = 0;
        for (Index i=0 ; i < Mats_.cols() ; i++) ncols += Mats_[i].cols();
        return ncols; 
    }

protected:
    std::vector < Matrix > Mats_;
}; // HNMatrix

/*! Block matrix with arbitrary number of matrices of the same type pasted vertically. */
template< class Matrix > class VNMatrix : public MatrixBase{ };
/*! Block diagonal matrix with arbitrary number of matrices of the same type. */
template< class Matrix > class DNMatrix : public MatrixBase{ };

/*! Block matrix with one matrix repeatedly pasted horizontally. */
template< class Matrix > class HRMatrix : public MatrixBase{ 
public:
    /*! Constructors */
    HRMatrix(){} //useful?
    HRMatrix(const Matrix & Mat) : Mat_(Mat){}
    HRMatrix(const Matrix & Mat, Index nmats) : Mat_(Mat), nmats_(nmats){}
    
    ~HRMatrix(){}
    
    void setNumber(Index nmats){ nmats_ = nmats; }
    /*! Return rows and columns. */
    virtual Index rows() const { return Mat_.rows(); }
    virtual Index cols() const { return Mat_.cols() * nmats_; }

    inline const Matrix & Mat() const { return Mat_; }
    inline Matrix & Mat() { return Mat_; }
protected:
    Matrix Mat_;
    Index nmats_;
}; // HRMatrix

/*! Block matrix with one matrix repeatedly pasted vertically. */
template< class Matrix > class VRMatrix : public MatrixBase{ };
/*! Block diagonal matrix with one matrix repeatedly as diagonals. */
template< class Matrix > class DRMatrix : public MatrixBase{ };

//** specializations better moved into gimli.h
//! nomenclature: Type(R/S)+Alignment(H/V/D)+Number(2/N/R)+Matrix
typedef H2Matrix< RMatrix, RMatrix > RH2Matrix;
typedef H2Matrix< SparseMapMatrix< double, Index >, SparseMapMatrix< double, Index > > SH2Matrix;
typedef DRMatrix< RMatrix > RDRMatrix;
typedef DRMatrix< SparseMapMatrix< double, Index > > SDRMatrix;

//! Some examples useful for special inversions
typedef H2Matrix< IdentityMatrix, IdentityMatrix > TwoModelsCMatrix; // -I +I
typedef DRMatrix< TwoModelsCMatrix > ManyModelsCMatrix;  //** not really diagonal!
typedef DRMatrix< RSparseMapMatrix > ManyCMatrix;
typedef V2Matrix< ManyCMatrix, ManyModelsCMatrix > MMMatrix; // Multiple models (LCI,timelapse)

} // namespace GIMLI

#endif