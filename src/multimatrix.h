/***************************************************************************
 *   Copyright (C) 2012      by the resistivity.net development team       *
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

#ifndef GIMLI_MULTIMATRIX__H
#define GIMLI_MULTIMATRIX__H

#include "gimli.h"
#include "sparsematrix.h"
#include "matrix.h"

//! multi matrices 
namespace GIMLI{
    
   
template < class MatType > class MultiMatrixEntry{
public:
    
    MultiMatrixEntry( Index rowStart, Index colStart_ ) : rowStart_( rowStart ), colStart_( colStart ) { 
        mat_ = new MatType( ); 
    }
    
    ~MultiMatrixEntry(){ 
        delete mat_; 
    }
    
    Index rowStart() const { return rowStart_; }
    
    Index rowEnd() const { return rowStart_ + mat_->rows(); }
    
    Index colStart() const { return colStart_; }
    
    Index colEnd() const { return colStart_ + mat_->cols(); }
    
    const MatrixBase & mat() const { return *mat_; }
    
    MatrixBase * mat() { return mat_; }
    
protected:
   
    MatrixBase * mat_;
    //! starting row;
    Index rowStart_;
    //! starting col;
    Index colStart_;
};

template < class ValueType > class MultiMatrix : public MatrixBase{
public:
    
    MultiMatrix(){}
    
    virtual ~MultiMatrix(){
        for ( Index i = 0; i < mats_.size(); i ++ ){
            delete mats_[ i ];
        }
    }
    
    MatrixBase * createMatrix( uint8 matrixRtti, Index startRow, Index startCol ){
        
        MultiMatrixEntry< MatrixBase > * matEntry;
        
        switch( matrixRtti ){
            case GIMLI_MATRIX_RTTI :
                matEntry = new MultiMatrixEntry< Matrix< ValueType > >( startRow, startCol );
                break;
            case GIMLI_SPARSEMAPMATRIX_RTTI :
                matEntry = new MultiMatrixEntry< SparseMapMatrix< ValueType, Index > >( startRow, startCol );
                break;
        }
        
        mats_.push_back( matEntry );
        return matEntry->mat();
    }
    
    virtual Index rows() const { 
        Index maxRows = 0;
        for ( Index i = 0; i < mats_.size(); i ++ ){
            maxRows = max ( maxRows, mats_[ i ]->rowEnd() );
        }
        //return std::max( mats_.begin(), mats_.end(), 0, std::mem_fun( &MultiMatrixEntry::rows ) );
        return maxRows;
    }
    
    virtual Index cols() const {  
        Index maxCols= 0;
        for ( Index i = 0; i < mats_.size(); i ++ ){
            maxCols = max ( maxCols , mats_[ i ]->colEnd() );
        }
        return maxCols;
    }
    
    virtual void clear() { 
        
        
    }
    
    virtual RVector mult( const RVector & a ) const {
        THROW_TO_IMPL
        return RVector(0);
    }
    
    virtual RVector transMult( const RVector & a ) const {
        THROW_TO_IMPL
        return RVector(0);
    }
    
protected:
    std::vector < MultiMatrixEntry< MatrixBase > * > mats_;
    
}; 

} // namespace GIMLI

#endif