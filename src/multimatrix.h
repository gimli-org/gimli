/******************************************************************************
 *   Copyright (C) 2012-2018 by the GIMLi development team                    *
 *   Thomas Günther  thomas@resistivity.net                                   *
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
