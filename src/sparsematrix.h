/***************************************************************************
 *   Copyright (C) 2006-2011 by the resistivity.net development team       *
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

#ifndef GIMLI_SPARSEMATRIX__H
#define GIMLI_SPARSEMATRIX__H

#include "gimli.h"
#include "elementmatrix.h"
#include "vector.h"
#include "vectortemplates.h"
#include "matrix.h"
#include "mesh.h"
#include "node.h"
#include "meshentities.h"

#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <cmath>

namespace GIMLI{

#define SPARSE_NOT_VALID throwError( EXIT_SPARSE_INVALID, WHERE_AM_I + " no data/or sparsity pattern defined.");

//! based on: Ulrich Breymann, Addison Wesley Longman 2000 , revised edition ISBN 0-201-67488-2, Designing Components with the C++ STL
template< class ValueType, class IndexType, class ContainerType > class MatrixElement {
public:
  typedef std::pair< IndexType, IndexType > IndexPair;
  typedef MatrixElement< ValueType, IndexType, ContainerType > & Reference;

  MatrixElement( ContainerType & Cont, IndexType r, IndexType c )
    : C( Cont ), I( C.find( IndexPair( r, c ) ) ), row( r ), column( c ) {
  }

  /* An assignment operator is required which in turn requires a
     reference to an object of type MatrixElement. When both the
     left- and right-hand side are identical, nothing has to happen.
     Otherwise, as above, it has to be checked whether the value of
     the right-hand element is 0 or not. The resulting behavior is
     described together with the above assignment operator, so that
     here it is simply called: */

    Reference operator = ( const Reference rhs ) {
        if ( this != & rhs ) {    // not identical?
            return operator = ( rhs.asValue() );  // see above
        } return *this;
    }

  /* The constructor initializes the private variables with all
     information that is needed. The container itself is located in
     the sparseMatrix class; here, the reference to it is entered.
     If the passed indices for row and column belong to an element
     not yet stored in the container, the iterator has the value
     C.end(). */

  ValueType asValue() const {
    if ( I == C.end() ) return ValueType( 0 ); else return ( *I ).second;
  }

  //** type conversion operator;
  operator ValueType () const { return asValue();  }

  /* According to the definition of the sparse matrix, 0 is returned
     if the element is not present in the container. Otherwise, the
     result is the second part of the object of type value_type
     stored in the container. */

  Reference operator = ( const ValueType & x ) {
    // not equal 0?
    if ( x != ValueType( 0 ) ) {
      /* If the element does not yet exist, it is put, together
	 with the indices, into an object of type value_type and
	 inserted with insert(): */

      if ( I == C.end() ) {
        assert( C.size() < C.max_size() );
        I = ( C.insert( typename ContainerType::value_type( IndexPair( row, column ), x ) ) ).first;
      } else ( *I ).second = x;
    }

    /* insert() returns a pair whose first part is an iterator
       pointing to the inserted object. The second part is of type
       bool and indicates whether the insertion took place because
       no element with this key existed. This is, however, not
       evaluated here because, due to the precondition (I ==
       C.end()), the second part must always have the value true.
       If, instead, the element already exists, the value is
       entered into the second part of the value_type object. If
       the value is equal 0, in order to save space the element is
       deleted if it existed. */

    else {                    // x = 0
      if ( I != C.end() ) {
        C.erase( I );
        I = C.end();
      }
    } return *this;
  }

  Reference operator += ( const ValueType & x ) {
    if ( x != ValueType( 0 ) ) {
      if ( I == C.end() ) {
        assert( C.size() < C.max_size() );
        I = ( C.insert( typename ContainerType::value_type( IndexPair( row, column ), x ) ) ).first;
      } else ( *I ).second += x;
    }
    return *this;
  }

  Reference operator -= ( const ValueType & x ) {
    if ( x != ValueType( 0 ) ) {
      if ( I == C.end() ) {
        assert( C.size() < C.max_size() );
        I = ( C.insert( typename ContainerType::value_type( IndexPair( row, column ), -x ) ) ).first;
      } else ( *I ).second -= x;
    }
    return *this;
  }

private:
  ContainerType & C;
  typename ContainerType::iterator I;
  IndexType row, column;

};  // class MatrixElement


//! based on: Ulrich Breymann, Addison Wesley Longman 2000 , revised edition ISBN 0-201-67488-2, Designing Components with the C++ STL
template< class ValueType, class IndexType > class SparseMapMatrix : public MatrixBase {
public:
    typedef std::pair< IndexType, IndexType > IndexPair;
    typedef std::map< IndexPair, ValueType, std::less< IndexPair > > ContainerType;
    typedef typename ContainerType::iterator          iterator;
    typedef typename ContainerType::const_iterator    const_iterator;
    typedef MatrixElement< ValueType, IndexType, ContainerType > MatElement;

    SparseMapMatrix( IndexType r = 0, IndexType c = 0 ) : rows_( r ), cols_( c ) {   }

    SparseMapMatrix( const std::string & filename ){
        this->load( filename );
    }

    SparseMapMatrix( const SparseMatrix< ValueType > & S ){
        this->copy( S );
    }

    SparseMapMatrix & operator = ( const SparseMatrix< ValueType > & S ){
        this->copy( S );
        return *this;
    }

    void copy( const SparseMatrix< ValueType > & S ){
        clear();
        cols_ = S.cols();
        rows_ = S.rows();

        std::vector < int > colPtr( S.vecColPtr() );
        std::vector < int > rowIdx( S.vecRowIdx() );
        Vector < ValueType > vals( S.vecVals() );

        for ( Index i = 0; i < S.size(); i++ ){
            for ( int j = colPtr[ i ]; j < colPtr[ i + 1 ]; j ++ ){
                (*this)[ i ][ rowIdx[ j ] ] = vals[ j ];
            }
        }
    }

    virtual void clear() {
        C_.clear();
        cols_ = 0; rows_ = 0;
    }

    inline void setRows( IndexType r ) { rows_ = r ; }
    inline IndexType rows()     const { return rows_; }

    inline void setCols( IndexType c ) { cols_ = c ; }
    inline IndexType cols()     const { return cols_; }

    inline IndexType columns()  const { return cols_; }

    inline IndexType size()     const { return C_.size(); }
    inline IndexType max_size() const { return C_.max_size(); }
    inline IndexType nVals()    const { return C_.size(); }

    inline iterator begin() { return C_.begin(); }
    inline iterator end() { return C_.end(); }

    inline const_iterator begin() const { return C_.begin(); }
    inline const_iterator end()   const { return C_.end(); }

    void addToCol( Index id, const ElementMatrix < double > & A ){
        for ( Index i = 0, imax = A.size(); i < imax; i++ ){
            (*this)[ A.idx( i ) ][ id ] += (ValueType) A.getVal( 0, i );
        }
    }

    void addToRow( Index id, const ElementMatrix < double > & A ){
        for ( Index i = 0, imax = A.size(); i < imax; i++ ){
            (*this)[ id ][ A.idx( i ) ] += (ValueType)A.getVal( 0, i );
        }
    }

    //SparseMatrix< T > & operator += ( const ElementMatrix < T > & A ){
    void operator += ( const ElementMatrix < double >& A ){
        for ( Index i = 0, imax = A.size(); i < imax; i++ ){
            for ( Index j = 0, jmax = A.size(); j < jmax; j++ ){
                (*this)[ A.idx( i ) ][ A.idx( j ) ] += (ValueType)A.getVal( i, j );
            }
        }
    }

    class Aux {  // for index operator below
    public:
        Aux( IndexType r, IndexType maxs, ContainerType & Cont )
            : Row( r ), maxColumns( maxs ), C( Cont ) { }

        MatElement operator [] ( IndexType c ) {
            if ( c < 0 || c >= maxColumns ){
                throwLengthError( EXIT_SPARSE_SIZE,
                                  WHERE_AM_I + " idx = " + toStr( c ) + " maxcol = "
                                  + toStr( maxColumns ) );
            }
            return MatElement( C, Row, c );
        }
    protected:
        IndexType Row, maxColumns;
        ContainerType & C;
    };

    Aux operator [] ( IndexType r ) {
        if ( r < 0 || r >= rows_ ){
            throwLengthError( EXIT_SPARSE_SIZE,
                              WHERE_AM_I + " idx = " + toStr( r ) + " maxrow = "
                              + toStr( rows_ ) );
        }
        return Aux( r, columns(), C_ );
    }

    inline IndexType idx1( const const_iterator & I ) const { return ( *I ).first.first; }
    inline IndexType idx2( const const_iterator & I ) const { return ( *I ).first.second; }

//     inline IndexType idx1( iterator & I ) { return ( *I ).first.first; }
//     inline IndexType idx2( iterator & I ) { return ( *I ).first.second; }

    inline const ValueType & val( const const_iterator & I ) const { return (*I).second;  }

    inline ValueType & val( const iterator & I ) { return (*I).second;  }

    inline ValueType getVal( IndexType i, IndexType j ) { return (*this)[ i ][ j ]; }

    inline void setVal( IndexType i, IndexType j, const ValueType & val ) {
        if ( ( i >= 0 && i < rows_ ) && ( j >=0 && j < cols_ ) ) {
            (*this)[ i ][ j ] = val;
        } else {
            throwLengthError( EXIT_SPARSE_SIZE,
                              WHERE_AM_I +
                              " i = " + toStr( i ) + " max_row = " + toStr( rows_ ) +
                              " j = " + toStr( j ) + " max_col = " + toStr( cols_ )
                              );
        }
    }
    inline void addVal( IndexType i, IndexType j, const ValueType & val ) {
        if ( ( i >= 0 && i < rows_ ) && ( j >=0 && j < cols_ ) ) {
            (*this)[ i ][ j ] += val;
        } else {
            throwLengthError( EXIT_SPARSE_SIZE,
                              WHERE_AM_I +
                              " i = " + toStr( i ) + " max_row = " + toStr( rows_ ) +
                              " j = " + toStr( j ) + " max_col = " + toStr( cols_ )
                              );
        }
    }

    
    /*! Return this * a  */
    virtual RVector mult( const RVector & a ) const {
        RVector ret( this->rows(), 0.0 );
        THROW_TO_IMPL
        return ret;
    }
    
    /*! Return this.T * a */
    virtual RVector transMult( const RVector & a ) const {
        RVector ret( this->cols(), 0.0 );
        THROW_TO_IMPL
        return ret;
    }
    
//     /*! Return this * a  */
//     virtual Vector < ValueType > mult( const Vector < ValueType > & a ) const {
//         Vector < ValueType > ret( this->rows(), 0.0 );
//         THROW_TO_IMPL
//         return ret;
//     }
//     
//     /*! Return this.T * a */
//     virtual Vector < ValueType > transMult( const Vector < ValueType > & a ) const {
//         Vector < ValueType > ret( this->cols(), 0.0 );
//         THROW_TO_IMPL
//         return ret;
//     }
    
    void save( const std::string & filename ) const {
        std::fstream file; openOutFile( filename, &file );

        for ( const_iterator it = begin(); it != end(); it ++ ){
            file <<  idx1( it ) << " " << idx2( it ) << " " << val( it ) << std::endl;
        }

        file.close();
    }

    void load( const std::string & filename ){
        std::fstream file; openInFile( filename, &file );
        std::vector < IndexType > vi, vj;
        std::vector < ValueType > vval;
        IndexType i, j;
        ValueType val;
        while ( file >> i >> j >> val ){
            vi.push_back( i );
            vj.push_back( j );
            vval.push_back( val );
        }
        file.close();
        setRows( (IndexType)max( vi ) + 1 );
        setCols( (IndexType)max( vj ) + 1 );

        for ( Index i = 0; i < vi.size(); i ++ ){
            (*this)[ vi[ i ] ][ vj[ i ] ] = vval[ i ];
        }
    }

protected:

  IndexType rows_, cols_;
  ContainerType C_;

};// class SparseMapMatrix

template < class ValueType, class IndexType >
int save( const SparseMapMatrix< ValueType, IndexType > & S, const std::string & fname ){
    return S.save( fname );
}

template < class ValueType, class IndexType >
int load( SparseMapMatrix< ValueType, IndexType > & S, const std::string & fname ){
    return S.load( fname );
}

/*! Scales a matrix A from left and right vectors such that A -> diag(l) * A * diag(r) */
template< class Vec >
void scaleMatrix ( SparseMapMatrix< double, Index > & S, const Vec & l, const Vec & r ) {

    if ( S.cols() != r.size() )
        throwLengthError(EXIT_SPARSE_SIZE, WHERE_AM_I + " " + toStr( S.cols() )
                                            + " != " + toStr( r.size() ) );
    if ( S.rows() != l.size() )
        throwLengthError(EXIT_SPARSE_SIZE, WHERE_AM_I + " " + toStr( S.rows() )
                                            + " != " + toStr( l.size() ) );

    for ( SparseMapMatrix< double, Index >::iterator it = S.begin(); it != S.end(); it ++ ){
                S.val( it ) *= l[ S.idx1( it ) ] * r[ S.idx2( it ) ];
//       int i = S.idx1( it );
//     int j = S.idx2( it );
//     S[ i ][ j ] *= ( l[ i ] * r[ j ] );
    }
    return;
}

/*! Performs a rank 1 update of a matrix such that A -> A + u * v^T */
template< class Vec >
void rank1Update ( SparseMapMatrix< double, Index > & S, const Vec & u, const Vec & v ) {

    if ( S.cols() != v.size() )
        throwLengthError(EXIT_SPARSE_SIZE, WHERE_AM_I + " " + toStr( S.cols() )
                                + " != " + toStr( v.size() ) );
    if ( S.rows() != u.size() )
        throwLengthError(EXIT_SPARSE_SIZE, WHERE_AM_I + " " + toStr( S.rows() )
                                + " != " + toStr( u.size() ) );

    for ( SparseMapMatrix< double, Index >::iterator it = S.begin(); it != S.end(); it ++ ){
        S.val( it ) += u[ S.idx1( it ) ] * v[ S.idx2( it ) ];
    }
    return;
}

// template < class ValueType, class IndexType, class V2, class T, class A >
// Vector < V2 > operator * ( const SparseMapMatrix< ValueType, IndexType > & S,
//                            const VectorExpr< T, A > & a ) {
//     return S * Vector< V2 >( a );
// }

template < class ValueType, class IndexType, class V2 >
Vector < V2 > operator * ( const SparseMapMatrix< ValueType, IndexType > & S,
                           const Vector< V2 > & a ) {

    typedef std::pair< IndexType, IndexType > IndexPair;
    typedef std::map< IndexPair, ValueType, std::less< IndexPair > > ContainerType;
    typedef typename ContainerType::iterator          iterator;
    typedef typename ContainerType::const_iterator    const_iterator;

    if ( S.cols() != a.size() ){
        throwLengthError(EXIT_SPARSE_SIZE, WHERE_AM_I + " S.cols() " + toStr( S.cols() ) + " != a.size() " + toStr( a.size() ) );
    }
    Vector < V2 > ret( S.rows() );

    for ( const_iterator it = S.begin(); it != S.end(); it ++ ){
        //tmp[ S.idx1( it ) ] += S.val( it ) * a[ S.idx2( it ) ];
        ret[ it->first.first ] += a[ it->first.second ] * it->second;
    }
    return ret;
}

// template < class ValueType, class IndexType, class V2, class T, class A >
// Vector < V2 > transMult( const SparseMapMatrix< ValueType, IndexType > & S,
//                            const VectorExpr< T, A > & a ) {
//     return transMult(S, Vector< V2 >( a ) );
// }

template < class ValueType, class IndexType, class V2 >
Vector < V2 > transMult( const SparseMapMatrix< ValueType, IndexType > & S,
                         const Vector < V2 > & a ){

    typedef std::pair< IndexType, IndexType > IndexPair;
    typedef std::map< IndexPair, ValueType, std::less< IndexPair > > ContainerType;
    typedef typename ContainerType::iterator          iterator;
    typedef typename ContainerType::const_iterator    const_iterator;

    if ( S.rows() != a.size() ){
        throwLengthError(EXIT_SPARSE_SIZE, WHERE_AM_I + " " + toStr( S.rows() ) + " != " + toStr( a.size() ) );
    }
    Vector < V2 > tmp( S.cols(), 0.0 );

    for ( const_iterator it = S.begin(); it != S.end(); it ++ ){
//         tmp[ S.idx2( it ) ] += S.val( it ) * a[ S.idx1( it ) ];
        tmp[ it->first.second ] += a[ it->first.first ] * it->second;

    }
    return tmp;
}

// template < class Vec > Vec transMult( const SparseMapMatrix< double, uint > & S, const Vec & a ){
//     if ( S.rows() != a.size() ){
//         throwLengthError(EXIT_SPARSE_SIZE, WHERE_AM_I + " " + toStr( S.cols() )
//                 + " != " + toStr( a.size() ) );
//     }
//     Vec tmp( S.cols() );
//
//     for ( SparseMapMatrix< double, uint >::const_iterator it = S.begin(); it != S.end(); it ++ ){
//         tmp[ S.idx2( it ) ] += S.val( it ) * a[ S.idx1( it ) ];
//     }
//     return tmp;
// }

//! Sparse matrix in compressed column storage (CCS) form
template < class ValueType > class SparseMatrix {
public:

  /*! Default constructor. Builds invalid sparse matrix */
    SparseMatrix() : valid_( false ) { }

    /*! Copy constructor. */
    SparseMatrix( const SparseMatrix < ValueType > & S )
        : colPtr_( S.vecColPtr() ), rowIdx_( S.vecRowIdx() ), vals_( S.vecVals() ), valid_( true ){
    }

    /*! Create Sparsematrix from c-arrays. Cant check for valid ranges, so please be carefull. */
    SparseMatrix( uint dim, Index * colPtr, Index nVals, Index * rowIdx, ValueType * vals ){
        colPtr_.reserve( dim + 1);
        colPtr_.resize( dim + 1);

        rowIdx_.reserve( nVals );
        rowIdx_.resize( nVals );

        vals_.resize( nVals );

//         std::copy( colPtr[ 0 ], colPtr[ dim ], colPtr_[ 0 ] );
//         std::copy( rowIdx[ 0 ], rowIdx[ nVals - 1 ], rowIdx_[ 0 ] );
//         std::copy( vals[ 0 ], vals[ nVals - 1 ], & vals_[ 0 ] );
    }

    /*! Destructor */
    ~SparseMatrix(){ }

    /*! Copy assignment operator. */
    SparseMatrix < ValueType > & operator = ( const SparseMatrix < ValueType > & S){
        if ( this != &S ){
            colPtr_ = S.vecColPtr();
            rowIdx_ = S.vecRowIdx();
            vals_   = S.vecVals();
            valid_  = true;
        } return *this;
    }

    SparseMatrix( const DSparseMapMatrix & S ){
        this->copy_( S );
    }

    SparseMatrix & operator = ( const DSparseMapMatrix & S ){
        this->copy_( S );
        return *this;
    }

    #define DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__( OP, FUNCT ) \
        void FUNCT( int i, int j, ValueType val ){ \
            if ( ::fabs( val ) > TOLERANCE ){ \
                for ( int k = colPtr_[ i ]; k < colPtr_[ i + 1 ]; k ++ ){ \
                    if ( rowIdx_[ k ] == j ) { \
                        vals_[ k ] OP##= val; return; \
                    } \
                } \
                std::cerr << WHERE_AM_I << " pos " << i << " " << j << " is not part of the sparsity pattern " << std::endl; \
            } \
        }\
        SparseMatrix< ValueType > & operator OP##= ( const ValueType & v ){\
            vals_ OP##= v; \
            return *this;\
        }\

        DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__( +, addVal )
        DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__( -, subVal )
        DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__( *, mulVal )
        DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__( /, divVal )

    #undef DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__

    SparseMatrix< ValueType > & operator += ( const SparseMatrix< ValueType > & A ){
        vals_ += A.vecVals();
        return *this;
    }

    SparseMatrix< ValueType > & operator += ( const ElementMatrix< double > & A ){
        if ( !valid_ ) SPARSE_NOT_VALID;
        for ( Index i = 0, imax = A.size(); i < imax; i++ ){
            for ( Index j = 0, jmax = A.size(); j < jmax; j++ ){
                addVal( A.idx( i ), A.idx( j ), A.getVal( i, j ) );
            }
        }
        return *this;
    }

    void clean(){ for ( Index i = 0, imax = nVals(); i < imax; i++ ) vals_[ i ] = (ValueType)( 0 ); }

    void clear(){
        colPtr_.clear();
        rowIdx_.clear();
        vals_.clear();
        valid_ = false;
    }

    void setVal( int i, int j, ValueType val ){
        if ( ::fabs( val ) > TOLERANCE ){
            for ( int k = colPtr_[ i ]; k < colPtr_[ i + 1 ]; k ++ ){
                if ( rowIdx_[ k ] == j ) {
                    vals_[ k ] = val; return;
                }
            }
            std::cerr << WHERE_AM_I << " pos " << i << " " << j << " is not part of the sparsity pattern " << std::endl;
        }
    }

    ValueType getVal( int i, int j ){
        for ( int k = colPtr_[ i ]; k < colPtr_[ i + 1 ]; k ++ ){
            if ( rowIdx_[ k ] == j ) {
                return vals_[ k ];
            }
        }
        std::cerr << WHERE_AM_I << " pos " << i << " "
                    << j << " is not part of the sparsity pattern " << std::endl;
        return (ValueType)-1.0;
    }

    void cleanRow( int i ){
         for ( int k = colPtr_[ i ]; k < colPtr_[ i + 1 ]; k ++ ){
            vals_[ k ] = (ValueType)0.0;
         }
    }

    void cleanCol( int i ){
        for ( int k = colPtr_[ i ]; k < colPtr_[ i + 1 ]; k ++ ){
            for ( int j = colPtr_[ rowIdx_[ k ] ]; j < colPtr_[ rowIdx_[ k ] + 1 ]; j ++ ){
                if ( rowIdx_[ j ] == i ) {
                    vals_[ j ] = ValueType(0);
                }
            }
        }
    }

    void copy_( const ISparseMapMatrix & S ){
        CERR_TO_IMPL
    }

    void copy_( const DSparseMapMatrix & S ){
        this->clear();
        Index col = 0, row = 0;
        double val;

        std::vector < std::map < Index, double > > idxMap( S.cols() );

        for ( DSparseMapMatrix::const_iterator it = S.begin(); it != S.end(); it ++ ){
            col = S.idx1( it );
            row = S.idx2( it );
            val = S.val( it );
            idxMap[ col ].insert( std::pair< Index, ValueType >( row, val ) );
        }

        colPtr_.reserve( S.cols() + 1);
        colPtr_.resize( S.cols() + 1);

        rowIdx_.reserve( S.nVals() );
        rowIdx_.resize( S.nVals() );

        vals_.resize( S.nVals() );

        colPtr_[ 0 ] = 0;

        Index colCounter = 0, rowCounter = 0;
        for ( std::vector < std::map < Index, double > >::iterator it = idxMap.begin(); it != idxMap.end(); it++ ){
            for ( std::map < Index, double >::iterator itR = (*it).begin(); itR != (*it).end(); itR++ ){
                rowIdx_[ rowCounter ] = itR->first;
                vals_[ rowCounter ] = (ValueType)itR->second;
                rowCounter ++;
            }
            colCounter ++;
            colPtr_[ colCounter ] = rowCounter;
        }
        valid_ = true;
    }

    void buildSparsityPattern( const Mesh & mesh ){
        colPtr_.resize( mesh.nodeCount() + 1 );

        Index col = 0, row = 0;
        std::vector < std::set < Index > > idxMap( mesh.nodeCount() );
        std::set < Index > tmp;
        Cell *cell;
        uint nc;
        for ( uint c = 0; c < mesh.cellCount(); c ++ ){
            cell = &mesh.cell( c );
            nc = cell->nodeCount();
            for ( uint i = 0; i < nc; i ++ ){
                for ( uint j = 0; j < nc; j ++ ){
                    row = cell->node( i ).id();
                    col = cell->node( j ).id();
                    idxMap[ col ].insert( row );
                }
            }
        }

        int nVals = 0;
        for ( std::vector < std::set < Index > >::iterator mIt = idxMap.begin(); mIt != idxMap.end(); mIt++ ){
            nVals += (*mIt).size();
        }

        rowIdx_.reserve( nVals );
        rowIdx_.resize( nVals );
        vals_.resize( nVals );

        colPtr_[ 0 ] = 0;
        Index k = 0;
        row = 0;
        for ( std::vector < std::set < Index > >::iterator mIt = idxMap.begin(); mIt != idxMap.end(); mIt++ ){
            for ( std::set< Index >::iterator sIt = (*mIt).begin(); sIt != (*mIt).end(); sIt++ ){
                rowIdx_[ k ] = (*sIt);
                vals_[ k ] = (ValueType)0.0;
                k++;
            }
            row++;
            colPtr_[ row ] = k;
        }
        valid_ = true;
    }

    inline int * colPtr() { if ( valid_ ) return &colPtr_[ 0 ]; else SPARSE_NOT_VALID;  return 0; }
    inline const int & colPtr() const { if ( valid_ ) return colPtr_[ 0 ]; else SPARSE_NOT_VALID; return colPtr_[ 0 ]; }
    inline const std::vector < int > & vecColPtr() const { return colPtr_; }

    inline int * rowIdx() { if ( valid_ ) return &rowIdx_[ 0 ]; else SPARSE_NOT_VALID; return 0; }
    inline const int & rowIdx() const { if ( valid_ ) return rowIdx_[ 0 ]; else SPARSE_NOT_VALID; return rowIdx_[ 0 ]; }
    inline const std::vector < int > & vecRowIdx() const { return rowIdx_; }

    inline ValueType * vals() { if ( valid_ ) return &vals_[ 0 ]; else SPARSE_NOT_VALID; return 0; }
    inline const ValueType & vals() const { if ( valid_ ) return vals_[ 0 ]; else SPARSE_NOT_VALID; return vals_[ 0 ]; }
    inline const Vector < ValueType > & vecVals() const { return vals_; }

    inline Index size() const { return colPtr_.size() - 1; }
    inline Index nVals() const { return vals_.size(); }
    inline Index cols() const { return size(); }
    inline Index rows() const { return size(); }

    int save( const std::string & fileName ){
        if ( !valid_ ) SPARSE_NOT_VALID;
        std::fstream file; if ( !openOutFile( fileName, & file ) ) return 0;

        file.setf( std::ios::scientific, std::ios::floatfield );
        file.precision( 14 );

        for ( uint i = 0; i < size(); i++ ){
            for ( Index j = colPtr_[ i ]; j < colPtr_[ i + 1 ]; j ++ ){
                file << i << "\t" << rowIdx_[ j ] << "\t" << vals_[ j ] << std::endl;
            }
        }
        file.close();
        return 1;
    }

    bool valid() const { return valid_; }

protected:

    // int to be cholmod compatible
    std::vector < int > colPtr_;
    std::vector < int > rowIdx_;
    Vector < ValueType > vals_;

    bool valid_;
};

template < class ValueType > SparseMatrix< ValueType > operator + ( const SparseMatrix< ValueType > & A, const SparseMatrix< ValueType > & B ){
    SparseMatrix< ValueType > ret( A );
    return ret += B;
}

template < class ValueType > Vector < ValueType > operator * ( const SparseMatrix < ValueType > & A
                                                             , const Vector < ValueType >& a ){

    if ( a.size() < A.size() ){
        throwLengthError( 1, WHERE_AM_I + " SparseMatrix size(): " + toStr( A.size() ) + " a.size(): " +
                                toStr( a.size() ) ) ;
    }
    Vector < ValueType > b( a.size() );
    for ( Index i = 0; i < A.size(); i++ ){
        for ( int j = A.vecColPtr()[ i ]; j < A.vecColPtr()[ i + 1 ]; j ++ ){
            b[ i ] += a[ A.vecRowIdx()[ j ] ] * A.vecVals()[ j ];
        }
    }
    return b;
}

} // namespace GIMLI

#endif //GIMLI__H
