/***************************************************************************
 *   Copyright (C) 2007-2011 by the resistivity.net development team       *
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

#ifndef _GIMLI_MATRIX__H
#define _GIMLI_MATRIX__H

#include "gimli.h"
#include "vector.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <cerrno>

#ifdef USE_THREADS
#ifdef HAVE_LIBBOOST_THREAD
#include <boost/thread.hpp>
#endif // HAVE_LIBBOOST_THREAD
#endif // USE_THREADS

namespace GIMLI{

//! Simple row-based dense matrix based on \ref Vector
/*! Simple row-based dense matrix based on \ref Vector */
template < class ValueType > class Matrix  {
public:
    /*! Constructs an empty matrix with the dimension rows x cols Content of the matrix is zero*/
    Matrix( size_t rows = 0, size_t cols = 0 ){
        resize( rows, cols );
    }

    /*! Copyconstructor */
    Matrix( const std::vector < Vector< ValueType > > & mat ){ copy_( mat ); }

    /*! Constructor, read matrix from file see \ref load( Matrix < ValueType > & A, const std::string & filename ). */
    Matrix( const std::string & filename ){ load( *this, filename ); }

    /*! Copyconstructor */
    Matrix( const Matrix < ValueType > & mat ){ copy_( mat ); }

    /*! Assignment operator */
    Matrix < ValueType > & operator = ( const Matrix< ValueType > & mat ){
        if ( this != & mat ){
            copy_( mat );
        } return *this;
    }

    /*! Destruct matrix and free memory. */
    ~Matrix(){}

    #define DEFINE_UNARY_MOD_OPERATOR__( OP, NAME ) \
    inline Matrix < ValueType > & operator OP##= ( const Matrix < ValueType > & A ) { \
        for ( register size_t i = 0; i < mat_.size(); i ++ ) mat_[ i ] OP##= A[ i ]; return *this; } \
    inline Matrix < ValueType > & operator OP##= ( const ValueType & val ) { \
        for ( register size_t i = 0; i < mat_.size(); i ++ ) mat_[ i ] OP##= val; return *this; } \

    DEFINE_UNARY_MOD_OPERATOR__( +, PLUS )
    DEFINE_UNARY_MOD_OPERATOR__( -, MINUS )
    DEFINE_UNARY_MOD_OPERATOR__( /, DIVID )
    DEFINE_UNARY_MOD_OPERATOR__( *, MULT )

    #undef DEFINE_UNARY_MOD_OPERATOR__

//     size_t col = cols();
//         for ( register size_t i = 0; i < mat_.size(); i ++ ) {
//             ValueType * Aj = &mat_[ i ][ 0 ];
//             ValueType * Aje = &mat_[ i ][ col ];
//             for ( ; Aj != Aje; ) *Aj++ OP##= val;
//         }   return *this; }

    /*! Readonly C style index operator, with boundary check.*/
    const Vector< ValueType > & operator [] ( size_t i ) const { return getVal( i ); }

    /*!*/
    Vector< ValueType > & operator [] ( size_t i ) {
        if ( i < 0 || i > this->rows()-1 ) {
            throwLengthError( 1, WHERE_AM_I + " row bounds out of range " +
                                toStr( i ) + " " + toStr( this->rows() ) ) ;
        }
        return mat_[ i ];
    }

    /*! Resize the matrix to rows x cols */
    inline void resize( size_t rows, size_t cols ){ allocate_( rows, cols ); }

    /*! Clear the matrix */
    inline void clear() { mat_.clear(); }

    /*! Return number of rows */
    inline size_t rows() const { return mat_.size(); }

    /*! Return number of cols */
    inline size_t cols() const { if ( mat_.size() > 0 ) return mat_[ 0 ].size(); return 0; }

    /*! Set a value. Throws out of range exception if index check fails. */
    inline void setVal( const Vector < ValueType > & val, size_t i ) {
        if ( i >= 0 && i < this->rows() ) {
            mat_[ i ] = val;
        } else {
            throwRangeError( 1, WHERE_AM_I, i, 0, this->rows() );
        }
    }

    /*! Readonly getter. */
    inline const Vector < ValueType > & getVal( size_t i ) const {
        if ( i < 0 || i > this->rows()-1 ) {
            throwLengthError( 1, WHERE_AM_I + " row bounds out of range " +
                                toStr( i ) + " " + toStr( this->rows() ) ) ;
        }
        return mat_[ i ];
    }

    /*! Return reference to row. Used for pygimli. */
    inline Vector < ValueType > & rowR( size_t i ) {
        if ( i < 0 || i > this->rows()-1 ) {
            throwLengthError( 1, WHERE_AM_I + " row bounds out of range " +
                                toStr( i ) + " " + toStr( this->rows() ) ) ;
        }
        return mat_[ i ];
    }

    /*! Readonly row of matrix, with boundary check.*/
    const Vector< ValueType > & row( size_t i ) const { return getVal( i ); }

    /*! Readonly col of matrix, with boundary check. Probably slow.*/
    Vector< ValueType > col( size_t i ) const {
        if ( i < 0 || i > this->cols()-1 ) {
            throwLengthError( 1, WHERE_AM_I + " col bounds out of range " +
                                toStr( i ) + " " + toStr( this->cols() ) ) ;
        }
        Vector < ValueType > col( this->rows() );
        for ( size_t j = 0, jmax = rows(); j < jmax; j ++ ) col[ j ] = mat_[ j ][ i ];
        return col;
    }

    inline void push_back( const Vector < ValueType > & vec ) {
        mat_.push_back( vec );
        rowFlag_.resize( rowFlag_.size() + 1 );
    }

    inline Vector< ValueType > & back() { return mat_.back(); }

    inline void setCol( uint col, const Vector < ValueType > & v ){
        if ( col < 0 || col > this->cols()-1 ) {
            throwLengthError( 1, WHERE_AM_I + " col bounds out of range " +
                                toStr( col ) + " " + toStr( this->cols() ) ) ;
        }
        if ( v.size() > this->rows() ) {
            throwLengthError( 1, WHERE_AM_I + " rows bounds out of range " +
                                toStr( v.size() ) + " " + toStr( this->rows() ) ) ;
        }
        for ( uint i = 0; i < v.size(); i ++ ) mat_[ i ][ col ] = v[ i ];
    }

    /*! Return reference to row flag vector. Maybee you can check if the rows are valid. Size is set automatic to the amount of rows. */
    BVector & rowFlag(){ return rowFlag_; }

protected:

    void allocate_( size_t rows, size_t cols ){
        if ( mat_.size() != rows ) mat_.resize( rows );
        for ( size_t i = 0; i < mat_.size(); i ++ ) mat_[ i ].resize( cols );
        rowFlag_.resize( rows );
    }

    void copy_( const Matrix < ValueType > & mat ){
        allocate_( mat.rows(), mat.cols() );
        for ( size_t i = 0; i < mat.rows(); i ++ ) mat_[ i ] = mat[ i ];
    }

    std::vector < Vector< ValueType > > mat_;

    /*! BVector flag(rows) for free use, e.g., check if rows are set valid. */
    BVector rowFlag_;
};

#define DEFINE_BINARY_OPERATOR__( OP, NAME ) \
template < class ValueType > \
Matrix < ValueType > operator OP ( const Matrix < ValueType > & A, const Matrix < ValueType > & B ) { \
Matrix < ValueType > tmp( A ); \
return tmp OP##= B; } \
template < class ValueType > \
Matrix < ValueType > operator OP ( const Matrix < ValueType > & A, const ValueType & v ) { \
Matrix < ValueType > tmp( A ); \
return tmp OP##= v; }

DEFINE_BINARY_OPERATOR__( +, PLUS )
DEFINE_BINARY_OPERATOR__( -, MINUS )
DEFINE_BINARY_OPERATOR__( /, DIVID )
DEFINE_BINARY_OPERATOR__( *, MULT )

#undef DEFINE_BINARY_OPERATOR__

template < class ValueType >
Vector < ValueType > operator * ( const Matrix < ValueType > & A, const Vector < ValueType > & b ){

    size_t cols = A.cols();
    size_t rows = A.rows();

    Vector < ValueType > ret( rows );

    register ValueType tmpval = 0;
    if ( b.size() == cols ){
        for ( register size_t i = 0; i < rows; ++i ){
            tmpval = 0;
            for ( register size_t j = 0; j < cols; ++j ) tmpval += A[ i ][ j ] * b[ j ];
            ret[ i ] = tmpval;
        }
//     for ( register size_t i = 0; i < rows; ++i ){
//       ret[ i ] = sum< T >( A[ i ] * v );
//     }
    } else {
        throwLengthError( 1, WHERE_AM_I + " " + toStr( cols ) + " != " + toStr( b.size() ) );
    }
    return ret;
}

template< class ValueType > class Mult{
public:
    Mult( Vector< ValueType > & x, const Vector< ValueType > & b, const Matrix < ValueType > & A, size_t start, size_t end ) :
        x_( &x ), b_( &b ), A_( &A), start_( start ), end_( end ){
    }
    void operator( )() {
        for ( register size_t i = start_; i < end_; i++ ) (*x_)[ i ] = sum( (*A_)[ i ] * *b_);
    }

    Vector< ValueType > * x_;
    const Vector< ValueType > * b_;
    const Matrix< ValueType > * A_;
    size_t start_;
    size_t end_;
};

template < class ValueType >
Vector < ValueType > multMT( const Matrix < ValueType > & A, const Vector < ValueType > & b ){
#ifdef USE_THREADS
    size_t cols = A.cols();
    size_t rows = A.rows();

    Vector < ValueType > ret( rows );
    boost::thread_group threads;
    size_t nThreads = 2;
    size_t singleCalcCount = size_t( ceil( (double)rows / (double)nThreads ) );
 //   CycleCounter cc;

    for ( size_t i = 0; i < nThreads; i ++ ){
// 	Vector < ValueType > *start = &A[ singleCalcCount * i ];
//         Vector < ValueType > *end   = &A[ singleCalcCount * ( i + 1 ) ];
        size_t start = singleCalcCount * i;
        size_t end   = singleCalcCount * ( i + 1 );
	if ( i == nThreads -1 ) end = A.rows();
//	cc.tic();
        threads.create_thread( Mult< ValueType >( ret, b, A, start, end ) );
  //      std::cout << cc.toc() << std::endl;
    }
    threads.join_all();
#else
    return mult( A, b );
#endif
}

template < class ValueType >
Vector < ValueType > mult( const Matrix < ValueType > & A, const Vector < ValueType > & b ){

    size_t cols = A.cols();
    size_t rows = A.rows();

    Vector < ValueType > ret( rows, 0.0 );

    register ValueType tmpval = 0;
    if ( b.size() == cols ){
//         for ( size_t i = 0; i < rows; ++i ){
//             tmpval = 0.0;
//             const ValueType * Aj = &A[i][0];
//             const ValueType * bj = &b[0];
//             const ValueType * bje = &b[ cols];
//            for ( ; bj != bje; ++ bj, ++ Aj) tmpval += *Aj * *bj;
//
//            //while (bj != bje) tmpval += *Aj++ * *bj++;
//
//             ret[ i ] = tmpval;
//         }
//           for (; __first != __last; ++__first)
// 	__init = __init + *__first;
//       return __init;
        for ( register size_t i = 0; i < rows; ++i ){
            ret[ i ] = sum( A[ i ] * b );
         }
    } else {
        throwLengthError( 1, WHERE_AM_I + " " + toStr( cols ) + " != " + toStr( b.size() ) );
    }
    return ret;
}

template < class ValueType >
bool operator == ( const Matrix< ValueType > & A, const Matrix< ValueType > & B ){
    if ( A.rows() != B.rows() || A.cols() != B.cols() ) return false;
    for ( size_t i = 0; i < A.rows(); i ++ ){
        if ( A[ i ] != B[ i ] ) return false;
    }
    return true;
}

template < class ValueType >
void scaleMatrix( Matrix < ValueType >& A,
                  const Vector < ValueType > & l, const Vector < ValueType > & r ){
    size_t rows = A.rows();
    size_t cols = A.cols();
    if ( rows != l.size() ){
        throwLengthError( 1, WHERE_AM_I + " " + toStr( rows ) + " != " + toStr( l.size() ) );
    };
    if ( cols != r.size() ){
        throwLengthError( 1, WHERE_AM_I + " " + toStr( cols ) + " != " + toStr( r.size() ) );
    }

    for ( size_t i = 0 ; i < rows ; i++ ) {
        //for ( size_t j = 0 ; j < cols ; j++ ) A[ i ][ j ] *= ( l[ i ] * r[ j ] );
        A[ i ] *= r * l[ i ];
    }
}

template < class ValueType >
void rank1Update( Matrix < ValueType > & A,
                  const Vector < ValueType > & u, const Vector < ValueType > & v ) {
    size_t rows = A.rows();
    size_t cols = A.cols();
    if ( rows != u.size() ){
        throwLengthError( 1, WHERE_AM_I + " " + toStr( rows ) + " != " + toStr( u.size() ) );
    };

    if ( cols != v.size() ){
        throwLengthError( 1, WHERE_AM_I + " " + toStr( cols ) + " != " + toStr( v.size() ) );
    }

    for ( size_t i = 0 ; i < rows ; i++ ) {
        //for ( size_t j = 0 ; j < ncols ; j++ ) A[ i ][ j ] += ( u[ i ] * v[ j ] );
        A[ i ] += v * u[ i ];
    }
    return;
}

template < class ValueType >
Matrix < ValueType > fliplr( const Matrix< ValueType > & m ){
    Matrix < ValueType > n;
    for ( size_t i = 0; i < m.rows(); i ++ ) n.push_back( fliplr( m[ i ] ) );
    return n;
}

//********************* MATRIX I/O *******************************

/*! Save matrix into a file (Binary).
    File suffix ($MATRIXBINSUFFIX) will be append if none given.
    Format: rows(uint32) cols(uint32) vals( rows*cols(ValueType) )
    If IOFormat == Ascii matrix will be saved in Ascii format, See: \ref saveMatrixRow
*/
template < class ValueType >
bool save( const Matrix < ValueType > & A, const std::string & filename, IOFormat format = Binary ){
    if ( format == Ascii ) return saveMatrixRow( A, filename );
    std::string fname( filename );
    if ( fname.rfind( '.' ) == std::string::npos ) fname += MATRIXBINSUFFIX;

    FILE *file; file = fopen( fname.c_str(), "w+b" );
    if ( !file ) {
        std::cerr << fname << ": " << strerror( errno ) << " " << errno << std::endl;
        return false;
    }

    uint32 rows = A.rows();
    uint ret = fwrite( & rows, sizeof(uint32), 1, file );
    uint32 cols = A.cols();
    ret = fwrite( & cols, sizeof(uint32), 1, file );

    for ( uint i = 0; i < rows; i ++ ){
        for ( uint j = 0; j < cols; j ++ ){
            ret = fwrite( &A[ i ][ j ], sizeof( ValueType ), 1, file );
        }
    }
    fclose( file );
    return true;
}

/*! Load matrix from a single or multiple files (Binary).
    File suffix (\ref MATRIXBINSUFFIX, ".matrix", ".mat") given or not -- loads single datafile, else try to load matrix from multiple binary vector files.
    Single format: see \ref save( const Matrix < ValueType > & A, const std::string & filename )
*/
template < class ValueType >
bool load( Matrix < ValueType > & A, const std::string & filename ){

    //!* First check if filename suffix is ".matrix", ".mat", \ref MATRIXBINSUFFIX;
    if ( filename.rfind( ".matrix" ) != std::string::npos ||
         filename.rfind( ".mat" ) != std::string::npos ||
         filename.rfind( MATRIXBINSUFFIX ) != std::string::npos ) {
        //!** yes, load \ref loadMatrixSingleBin( filename )
        return loadMatrixSingleBin( A, filename );
    }

    //!* no: check if filename is expandable with suffix ".matrix" or ".mat";
    if ( fileExist( filename + ".matrix" ) ) return loadMatrixSingleBin( A, filename + ".matrix" );
    if ( fileExist( filename + ".mat" ) )    return loadMatrixSingleBin( A, filename + ".mat" );
    if ( fileExist( filename + MATRIXBINSUFFIX ) )
        //!** yes , load \ref loadMatrixSingleBin( filename + \ref MATRIXBINSUFFIX )
        return loadMatrixSingleBin( A, filename + MATRIXBINSUFFIX );

    //!* no: try to load matrix from multiple binary vectors;
    return loadMatrixVectorsBin( A, filename );
}

/*! Force to load single matrix binary file.
    Format: see \ref save( const Matrix < ValueType > & A, const std::string & filename ). */
template < class ValueType >
bool loadMatrixSingleBin( Matrix < ValueType > & A,
                          const std::string & filename ){

    FILE *file;  file = fopen( filename.c_str(), "r+b" );
    if ( !file ) {
        throwError( EXIT_OPEN_FILE, WHERE_AM_I + " " + filename + ": " + strerror( errno ) );
    }

    uint ret = 0;
    int rows = 0; ret = fread( &rows, sizeof(uint32), 1, file );
    int cols = 0; ret = fread( &cols, sizeof(uint32), 1, file );

    A.resize( rows, cols );
    for ( int i = 0; i < rows; i ++ ){
        for ( int j = 0; j < cols; j ++ ){
            ret = fread( (char*)&A[ i ][ j ], sizeof( ValueType ), 1, file );
        }
    }
    fclose( file );
    A.rowFlag().fill( 1 );
    return true;
}

/*! Force to load multiple binary vector files into one matrix (row-based). File name will be determined from filenamebody + successive increased number (read while files exist). \n
e.g. read "filename.0.* ... filename.n.* -> Matrix[0--n)[ 0..vector.size() )\n
kCount can be given to use as subcounter. \n
e.g. read "filename.0_0.* ... filename.n_0.* ... filename.0_kCount-1.* ... filename.n_kCount-1.* ->
Matrix[0--n*kCount)[ 0..vector.size() )
*/
template < class ValueType >
bool loadMatrixVectorsBin( Matrix < ValueType > & A,
                            const std::string & filenameBody, uint kCount = 1 ){

    A.clear();
    Vector < ValueType > tmp;
    std::string filename;

    for ( uint i = 0; i < kCount; i++ ){
        uint count = 0;
        while ( 1 ){ // load as long as posible
            if ( kCount > 1 ){
                filename = filenameBody + "." + toStr( count ) + "_" + toStr( i ) + ".pot";
            } else {
                filename = filenameBody + "." + toStr( count ) + ".pot";
            }

            if ( !fileExist( filename ) ){
                filename = filenameBody + "." + toStr( count );
                if ( !fileExist( filename ) ){
                    if ( count == 0 ) {
	               std::cerr << " can not found: " << filename << std::endl;
                    }
                    break;
                }
            }
            if ( load( tmp, filename, Binary  ) ) A.push_back( tmp );
            count ++;
        } // while files exist
    } // for each k count
    return true;
}

/*! Save Matrix into Ascii File (column based) */
template < class ValueType >
bool saveMatrixCol( const Matrix < ValueType > & A, const std::string & filename ){
    return saveMatrixCol( A, filename, "" );
}

/*! Save Matrix into Ascii File (column based)  with optional comments header line*/
template < class ValueType >
bool saveMatrixCol( const Matrix < ValueType > & A, const std::string & filename,
                    const std::string & comments ){
    std::fstream file; openOutFile( filename, & file, true );
    if ( comments.length() > 0 ){
        file << "#" << comments << std::endl;
    }

    for ( uint i = 0; i < A.cols(); i ++ ){
        for ( uint j = 0; j < A.rows(); j ++ ){
            file << A[ j ][ i ] << "\t";
        }
        file << std::endl;
    }
    file.close();
    return true;
}

/*! Load Matrix from Ascii File (column based) */
template < class ValueType >
bool loadMatrixCol( Matrix < ValueType > & A, const std::string & filename ){
    std::vector < std::string > comments;
    return loadMatrixCol( A, filename, comments );
}

/*! Load Matrix from Ascii File (column based), with optional comments header line. */
template < class ValueType >
bool loadMatrixCol( Matrix < ValueType > & A, const std::string & filename,
                    std::vector < std::string > & comments ){

    size_t commentCount = 0;
    size_t cols = countColumnsInFile( filename, commentCount );
//     size_t rows = countRowsInFile( filename );
//     // get length of file:
//     std::fstream file; openInFile( filename, & file, true );
//     size_t length = fileLength( file );
//
//     // allocate memory:
//     char * buffer = new char[ length ];
//     file.read( buffer, length );
//     file.close();
//
//     delete buffer;
//     return true;

    std::vector < double > tmp;
    Vector < ValueType > row( cols );
    std::fstream file; openInFile( filename, & file, true );
    for ( uint i = 0; i < commentCount; i ++ ) {
        std::string str;
        getline( file, str );
        comments = getSubstrings( str.substr( str.find( '#' ), -1 ) );
    }

    double val;
    while( file >> val ) tmp.push_back( val );

    file.close();
    size_t rows = tmp.size() / cols ;
    A.resize( cols, rows );

    for ( uint i = 0; i < rows; i ++ ){
        for ( uint j = 0; j < cols; j ++ ){
            A[ j ][ i ] = tmp[ i * cols + j ];
        }
    }
    return true;
}

/*! Save Matrix into Ascii File (row based) */
template < class ValueType >
bool saveMatrixRow( const Matrix < ValueType > & A, const std::string & filename ){
    return saveMatrixRow( A, filename, "" );
}

/*! Save Matrix into Ascii File (row based)  with optional comments header line*/
template < class ValueType >
bool saveMatrixRow( const Matrix < ValueType > & A, const std::string & filename,
                    const std::string & comments ){
    std::fstream file; openOutFile( filename, & file, true );
    if ( comments.length() > 0 ){
        file << "#" << comments << std::endl;
    }

    for ( uint i = 0; i < A.rows(); i ++ ){
        for ( uint j = 0; j < A.cols(); j ++ ){
            file << A[ i ][ j ] << "\t";
        }
        file << std::endl;
    }
    file.close();
    return true;
}

/*! Load Matrix from Ascii File (row based) */
template < class ValueType >
bool loadMatrixRow( Matrix < ValueType > & A, const std::string & filename ){

    std::vector < std::string > comments;
    return loadMatrixRow( A, filename, comments );
}

/*! Load Matrix from Ascii File (row based), with optional comments header line. */
template < class ValueType >
bool loadMatrixRow( Matrix < ValueType > & A,
                    const std::string & filename,
                    std::vector < std::string > & comments ){

    size_t commentCount = 0;
    size_t cols = countColumnsInFile( filename, commentCount );

    Vector < ValueType > row( cols );
    std::fstream file; openInFile( filename, & file, true );
    for ( uint i = 0; i < commentCount; i ++ ) {
        std::string str;
        getline( file, str );
        comments = getSubstrings( str.substr( str.find( '#' ), -1 ) );
    }

    double val;
    std::vector < double > tmp;
    while( file >> val ) tmp.push_back( val );

    file.close();
    size_t rows = tmp.size() / cols ;
    A.resize( rows, cols );

    for ( uint i = 0; i < rows; i ++ ){
        for ( uint j = 0; j < cols; j ++ ){
            A[ i ][ j ] = tmp[ i * cols + j ];
        }
    }
    return true;
}

/*! Return determinant for Matrix(2 x 2)*/
template < class T > inline T det( const T & a, const T & b, const T & c, const T & d ){
    return a * d - b * c;
}

/*! Return determinant for Matrix A. This function is a stub. Only Matrix dimension of 2 and 3 are considered. */
template < class Matrix > double det( const Matrix & A ){
  //** das geht viel schoener, aber nicht mehr heute.;
  double det = 0.0;
  switch ( A.rows() ){
  case 2: det = A[ 0 ][ 0 ] * A[ 1 ][ 1 ] - A[ 0 ][ 1 ] * A[ 1 ][ 0 ];
    break;
  case 3:
    det = A[ 0 ][ 0 ] * ( A[ 1 ][ 1 ] * A[ 2 ][ 2 ] - A[ 1 ][ 2 ] * A[ 2 ][ 1 ] ) -
      A[ 0 ][ 1 ] * ( A[ 1 ][ 0 ] * A[ 2 ][ 2 ] - A[ 1 ][ 2 ] * A[ 2 ][ 0 ] ) +
      A[ 0 ][ 2 ] * ( A[ 1 ][ 0 ] * A[ 2 ][ 1 ] - A[ 1 ][ 1 ] * A[ 2 ][ 0 ] );
    break;
  default:
    std::cerr << WHERE_AM_I << " matrix determinant of dim not yet implemented -- dim: " << A.rows() << std::endl;
    break;
  }
  return det;
}
    
} //namespace GIMLI

#endif // _GIMLI_MATRIX__H

