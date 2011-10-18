/***************************************************************************
 *   Copyright (C) 2007-2011 by the resistivity.net development team       *
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

//! block matrices for easier inversion, see appendix E in GIMLi tutorial
namespace GIMLI{
    
/*! simple example for tutorial purposes */
/*! Sparse Block Matrix consisting of two horizontally pasted matrices */
class DLLEXPORT H2SparseMapMatrix {
public:
    H2SparseMapMatrix(){}
    ~H2SparseMapMatrix(){}
    
    inline const DSparseMapMatrix & H1() const { return H1_; }
    inline const DSparseMapMatrix & H2() const { return H2_; }

    inline DSparseMapMatrix & H1() { return H1_; }
    inline DSparseMapMatrix & H2() { return H2_; }

    inline size_t rows() const { return H1_.rows(); }
    inline size_t cols() const { return H1_.cols() + H2_.cols(); }
    
    inline void clear() { H1_.clear(); H2_.clear(); }
//    inline size_t size() const { return H1_.rows(); }
//    inline size_t size() { return H1_.size(); }
//    RVector operator[]( size_t n ){ return RVector( 0 ); }
protected:
    DSparseMapMatrix H1_, H2_; // create inplace (or better hold references of it?)
}; // class H2SparseMapMatrix

inline void rank1Update( H2SparseMapMatrix & A, const RVector & u, const RVector & v ) {
    CERR_TO_IMPL
    return;
}

inline bool save( const H2SparseMapMatrix & A, const std::string & filename, IOFormat format = Ascii ){
    CERR_TO_IMPL
    return false;
}
  
inline RVector operator * ( const H2SparseMapMatrix & A, const RVector & x ){
    return A.H1() * x + A.H2() * x;
}

inline RVector transMult( const H2SparseMapMatrix & A, const RVector & b ){
    return cat( transMult( A.H1(), b ), transMult( A.H2(), b ) );
}

template< class Matrix1, class Matrix2 > class H2Matrix{ };
template< class Matrix1, class Matrix2 > class V2Matrix{ };
template< class Matrix1, class Matrix2 > class D2Matrix{ };

template< class Matrix > class HNMatrix{ };
template< class Matrix > class VNMatrix{ };
template< class Matrix > class DNMatrix{ };

template< class Matrix > class HRMatrix{ };
template< class Matrix > class VRMatrix{ };
template< class Matrix > class DRMatrix{ };

//** specializations better moved into gimli.h
//! nomenclature: Type(R/S)+Alignment(H/V/D)+Number(2/N/R)+Matrix
typedef H2Matrix< RMatrix, RMatrix > RH2Matrix;
typedef H2Matrix< SparseMapMatrix< double, Index >, SparseMapMatrix< double, Index > > SH2Matrix;
typedef DRMatrix< RMatrix > RDRMatrix;
typedef DRMatrix< SparseMapMatrix< double, Index > > SDRMatrix;

} // namespace GIMLI

#endif