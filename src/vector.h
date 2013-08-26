/***************************************************************************
 *   Copyright (C) 2007-2012 by the resistivity.net development team       *
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

//** Idea taken from
//
// @Article{Veldhuizen95b,
//   author =       "Todd L. Veldhuizen",
//   title =        "Expression templates",
//   journal =      "C++ Report",
//   volume =       "7",
//   number =       "5",
//   pages =        "26--31",
//   month =        jun,
//   year =         "1995",
//   note =         "Reprinted in C++ Gems, ed. Stanley Lippman"
// }

#ifndef GIMLI_VECTOR__H
#define GIMLI_VECTOR__H

#define EXPRVEC_USE_TEMPORARY_EXPRESSION
//#define EXPRVEC_USE_STD_ALGORITHM
#define EXPRVEC_USE_INDIRECTION

#include "gimli.h"
#include "expressions.h"

#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstring>
#include <fstream>
#include <cerrno>
#include <iterator>

#ifdef USE_BOOST_BIND
#include <boost/bind.hpp>
#endif

namespace GIMLI{

typedef std::vector < Index > IndexArray;template < class ValueType, class A > class __VectorExpr;

template < class ValueType > class DLLEXPORT VectorIterator {
public:
    VectorIterator( ) : val_( NULL ), maxSize_( 0 ){ }

    VectorIterator(const VectorIterator < ValueType > & iter) : val_( iter.val_ ), maxSize_( iter.maxSize_ ){ }

    VectorIterator < ValueType > & operator = ( const VectorIterator < ValueType > & iter ){
        if ( this != & iter ){
            val_ = iter.val_;
            maxSize_ = iter.maxSize_ ;
        }
        return *this;
    }

    inline const ValueType & operator * () const { return * val_; }
    inline ValueType & operator * () { return * val_; }

    inline const ValueType & operator [] ( const Index i ) const { return val_[ i ]; }
    inline ValueType & operator [] ( const Index i ) { return val_[ i ]; }

    inline VectorIterator< ValueType > & operator ++ () { ++val_; return *this; }
    inline VectorIterator< ValueType > & operator -- () { --val_; return *this; }
	
    inline bool operator == ( const VectorIterator< ValueType > & a ){ return val_ == a.val_; }
    inline bool operator != ( const VectorIterator< ValueType > & a ){ return val_ != a.val_; }
 
    inline Index size() const { return maxSize_; }

    inline ValueType * ptr() const { return val_; }
    inline ValueType * ptr() { return val_; }

    ValueType * val_;
    Index maxSize_;
	
	#ifdef _MSC_VER
	// basic iterator traits needed by msvc
	typedef ValueType value_type;
	typedef Index difference_type;
	typedef value_type & pointer;    
	typedef value_type & reference;
	typedef std::random_access_iterator_tag iterator_category;
	bool operator < (const VectorIterator< ValueType > & a) const { return val_ < a.val_; }
	#endif
};

//! One dimensional array aka Vector of limited size.
/*!
One dimensional array aka Vector of limited size. Size limit depends on platform (32bit system maxsize = 2^32, 64bit system, maxsize=2^64)
*/
template< class ValueType > class DLLEXPORT Vector {
public:
    typedef ValueType ValType;
    /*!
     * Construct one-dimensional array of size n. The vector is cleaned (filled with zero)
     */
    Vector( Index n = 0 ) : data_( NULL ), begin_( NULL ), end_( NULL ) {
    // explicit Vector( Index n = 0 ) : data_( NULL ), begin_( NULL ), end_( NULL ) {
        allocate_( n );
        clean();
    }

    /*!
     * Construct one-dimensional array of size n, and fill it with val
     */
    Vector( Index n, const ValueType & val ) : data_( NULL ), begin_( NULL ), end_( NULL ) {
        allocate_( n );
        fill( val );
    }

    /*!
     * Construct vector from file. Shortcut for Vector::load
     */
    Vector( const std::string & filename, IOFormat format = Ascii )
        : size_( 0 ), data_( NULL ), begin_( NULL ), end_( NULL ) {
        this->load( filename, format );
    }

    /*!
     * Copy constructor. Create new vector as a deep copy of v.
     */
    Vector( const Vector< ValueType > & v ) : data_( NULL ), begin_( NULL ), end_( NULL ) {
        allocate_( v.size() );
        copy_( v );
    }

    /*!
     * Copy constructor. Create new vector as a deep copy of the slice v[start, end )
     */
    Vector( const Vector< ValueType > & v, Index start, Index end )
        : data_( NULL ), begin_( NULL ), end_( NULL ) {
        allocate_( end - start );
        std::copy( &v[ start ], &v[ end ], data_ );
    }

    /*!
     * Copy constructor. Create new vector from expression
     */
    template < class A > Vector( const __VectorExpr< ValueType, A > & v ) : data_( NULL ), begin_( NULL ), end_( NULL ) {
        allocate_( v.size() );
        assign_( v );
    }

    /*!
     * Copy constructor. Create new vector as a deep copy of std::vector( Valuetype )
     */
    template < class U > Vector( const std::vector< U > & v ) : data_( NULL ), begin_( NULL ), end_( NULL ) {
        allocate_( v.size() );
        for ( Index i = 0; i < v.size(); i ++ ) data_[ i ] = ValueType( v[ i ] );
        //std::copy( &v[ 0 ], &v[ v.size() ], data_ );
    }

    /*! Default destructor. */
    ~Vector() { free_(); }

    /*! Assignment operator. Creates a new vector as copy of v */
    Vector< ValueType > & operator = ( const Vector< ValueType > & v ) {
        if ( this != &v ) {
            resize( v.size() );
            copy_( v );
        }
        return *this;
    }

    /*! Assignment operator. Creates a new vector as from expression. */
    template < class A > Vector< ValueType > & operator = ( const __VectorExpr< ValueType, A > & v ) {
        assign_( v );
        return *this;
    }

    /*! Assignment operator. Fill the existing vector with val. Shortcut for fill. */
    Vector< ValueType > & operator = ( const ValueType & val ) {
        fill( val );
        return *this;
    }

    inline const ValueType & operator[]( const Index i ) const { return data_[ i ]; }

    inline ValueType & operator[]( const Index i ) { return data_[ i ]; }

    inline const Vector < ValueType > operator[]( const IndexArray & i ) const { return (*this)( i ); }
    
    inline Vector < ValueType > operator[]( const IndexArray & i ) { return (*this)( i ); }

     /*! Return a new vector that match the slice [start, end).  end == -1 or larger size() sets end = size.
        Throws exception on violating boundaries. */
    Vector < ValueType > operator () ( Index start, SIndex end ) const {
        Index e = (Index) end;
        if ( end == -1 || end > (SIndex)size_ ) e = size_;

        Vector < ValueType > v( end-start );
        if ( start >= 0 && start < e ){
            std::copy( & data_[ start ], & data_[ e ], &v[ 0 ] );
        } else {
            throwLengthError(1, WHERE_AM_I + " bounds out of range " +
                                str(start) + " " + str(end) + " " + str(size_));
        }
        return v;
    }

    /*!
     * Return a new vector that based on indices's. Throws exception if indices's are out of bound
     */
    Vector < ValueType > operator () ( const IndexArray & idx ) const {
        Vector < ValueType > v( idx.size() );
        Index id;
        for ( Index i = 0; i < idx.size(); i ++ ){
           id = idx[ i ];
           if ( id >= 0 && id < size_ ){
                v[ i ] = data_[ id ];
           } else {
                throwLengthError( 1, WHERE_AM_I + " idx out of range " +
                                     str( id ) + " [" + str( 0 ) + " " + str( size_ ) + ")" );
           }
        }
        return v;
    }

    Vector < ValueType > operator () ( const std::vector < int > & idx ) const {
        Vector < ValueType > v( idx.size() );
        Index id;
        for ( Index i = 0; i < idx.size(); i ++ ){
           id = idx[ i ];
           if ( id >= 0 && id < size_ ){
                v[ i ] = data_[ (Index)id ];
           } else {
                throwLengthError( 1, WHERE_AM_I + " idx out of range " +
                                     str( id ) + " [" + str( 0 ) + " " + str( size_ ) + ")" );
           }
        }
        return v;
    }
    
    /*!
    Implicite type converter
     */
    template < class T > operator Vector< T >( ){
        COUTMARKER
        Vector< T > f( this->size() );
        for ( uint i = 0; i < this->size(); i ++ ){ f[i] = T( data_[ i ] ); }
        return f;
    }

    /*! Set a value. Same as fill(val) */
    inline Vector< ValueType > & setVal( const ValueType & val ) {
        this->fill( val );
        return *this;
    }

    /*! Set a value at index i. Throws out of range exception if index check fails. */
    inline Vector< ValueType > & setVal( const ValueType & val, Index i ) {
        if ( i >= 0 && i < this->size() ) {
            data_[ i ] = val;
        } else {
            throwRangeError( 1, WHERE_AM_I, i, 0, this->size() );
        }
        return *this;
    }

    /*! Set a value at slice index [start, end). Throws out of range exception if index check fails.
        end will set to this->size() if larger or -1. */
    inline Vector< ValueType > & setVal( const ValueType & val, Index start, SIndex end ) {
        Index e = (Index)end;
        if ( end == -1 ) e = this->size();
        else if ( e > this->size() ) e = this->size();

        if ( start > e ) start = e;

        if ( start >= 0 && start < this->size() ) {
            std::fill( data_+ start, data_ + e, val );
        } else {
            throwRangeError( 1, WHERE_AM_I, start, 0, this->size() );
        }
        return *this;
    }

    /*! Set multiple values. Throws out of range exception if index check fails. */
    inline Vector< ValueType > & setVal( const ValueType & val, const std::vector < Index > & idx ) {
        for ( Index i = 0; i < idx.size(); i ++ ) setVal( val, idx[ i ] );
        return *this;
    }

    /*! Set multiple values. Throws out of range exception if index check fails. */
    inline Vector< ValueType > & setVal( const Vector < ValueType > & vals,
                                            const std::vector < Index > & idx ) {
        if ( idx.size() != vals.size() ){
            throwLengthError( 1, WHERE_AM_I + " idx.size() != vals.size() " +
                                toStr( idx.size() ) + " " + toStr( vals.size() ) );
        }
        for ( Index i = 0; i < idx.size(); i ++ ){
           data_[ idx[ i ] ] = vals[ i ];
//            setVal( vals[ i ], idx[ i ] );
        }
        return *this;
    }

    /*! Set values from slice. If vals.size() == this.size() copy vals[start, end) -> this[start, end) else
        assume vals is a slice itsself, so copy vals[0, end-start] -> this[start, end)
         if end larger than this size() sets end = size. Throws exception on violating boundaries. */
    inline Vector< ValueType > & setVal( const Vector < ValueType > & vals, Index start, Index end ) {
        if ( start > this->size() ){
            throwLengthError( 1, WHERE_AM_I + " vals.size() < start " +
                                toStr( vals.size() ) + " " + toStr( start ) + " " + toStr( end ) ) ;
        }
        
        if ( end > this->size() ) end = this->size();
        if ( start > end ) start = end;

        if ( vals.size() < ( end - start ) ){
            throwLengthError( 1, WHERE_AM_I + " vals.size() < ( end-start ) " +
                                toStr( vals.size() ) + " " + toStr( start ) + " " + toStr( end ) ) ;
        }

        if ( this->size() == vals.size() ){
            std::copy( &vals[ start ], &vals[ end ], &data_[ start ] );
        } else {
            std::copy( &vals[ 0 ], &vals[ end - start ], &data_[ start ] );
        }
        return *this;
    }

    /*! Like setVal( vals, start, end ) instead copy use += */
    inline Vector< ValueType > & addVal( const Vector < ValueType > & vals, Index start, Index end ) {
        if ( end > this->size() ) end = this->size();
        if ( start > end ) start = end;

        if ( vals.size() < end - start){
            throwLengthError( 1, WHERE_AM_I + " vals.size() < ( end-start ) " +
                                toStr( vals.size() ) + " " + toStr( start ) + " " + toStr( end ) ) ;
        }

        if ( this->size() == vals.size() ){
            for ( Index i = start; i < end; i ++ ) data_[ i ] += vals[ i ];
        } else {
            for ( Index i = start; i < end; i ++ ) data_[ i ] += vals[ i - start ];
        }

        return *this;
    }

    /*! Add values from vals id index idx. Throws out of range exception if index check fails. */
    inline Vector< ValueType > & addVal( const Vector < ValueType > & vals, const std::vector < Index > & idx ) {
        if ( idx.size() != vals.size() ){
            throwLengthError( 1, WHERE_AM_I + " idx.size() != vals.size() " +
                                toStr( idx.size() ) + " " + toStr( vals.size() ) );
        }
        for ( Index i = 0; i < idx.size(); i ++ ) data_[ idx[ i ] ] += vals[ i ];
        
        return *this;
    }
    
    /*! Get a value. Throws out of range exception if index check fails. */
    inline const ValueType & getVal( Index i ) const {
        if ( i >= 0 && i < this->size() ) {
            return data_[ i ];
        } else {
            throwRangeError( 1, WHERE_AM_I, i, 0, this->size() );
        }
        return data_[ 0 ];
    }

#ifdef PYGIMLI
//    needed for: /usr/include/boost/python/def_visitor.hpp
    bool operator < ( const Vector< ValueType > & v ) const { return false; }
#else
    BVector operator < ( const Vector< ValueType > & v ) const {
	BVector ret( this->size(), 0 );
	    if ( this->size() != v.size() ) {
		throwLengthError( 1, WHERE_AM_I + " array size unequal" +
                                toStr( this->size() ) + " != " + toStr( v.size() ) );
	    }
	for ( Index i = 0; i < v.size(); i ++ ) ret[ i ] = isLesser( data_[ i ], v[ i ] );
	return ret;
    }
#endif

#define DEFINE_COMPARE_OPERATOR_VEC__( OP, FUNCT ) \
    BVector operator OP ( const Vector< ValueType > & v ) const { \
        BVector ret( this->size(), 0 ); \
        if ( this->size() != v.size() ) { \
            throwLengthError( 1, WHERE_AM_I + " array size unequal " + \
                                toStr( this->size() ) + " != " + toStr( v.size() ) ); \
        } \
        for ( Index i = 0; i < v.size(); i ++ ) ret[ i ] = FUNCT( data_[ i ], v[ i ] ); \
        return ret; \
    } \

DEFINE_COMPARE_OPERATOR_VEC__( <=, isLesserEqual )
DEFINE_COMPARE_OPERATOR_VEC__( >=, isGreaterEqual )
DEFINE_COMPARE_OPERATOR_VEC__( >, isGreater )

#undef DEFINE_COMPARE_OPERATOR_VEC__
   
#define DEFINE_COMPARE_OPERATOR__( OP, FUNCT ) \
    inline BVector operator OP ( const ValueType & v ) const { \
        BVector ret( this->size(), 0 ); \
        for ( Index i = 0; i < this->size(); i ++ ){ ret[ i ] = FUNCT( data_[ i ], v ); } \
        return ret;\
    } \

DEFINE_COMPARE_OPERATOR__( <, isLesser )
DEFINE_COMPARE_OPERATOR__( <=, isLesserEqual )
DEFINE_COMPARE_OPERATOR__( >=, isGreaterEqual )
DEFINE_COMPARE_OPERATOR__( ==, isEqual )
DEFINE_COMPARE_OPERATOR__( !=, isNonEqual )
DEFINE_COMPARE_OPERATOR__( >, isGreater )

#undef DEFINE_COMPARE_OPERATOR__

//     inline BVector operator > ( const ValueType & v ) const { 
//         BVector ret( this->size(), 0 ); 
//         std::transform( data_, data_ + size_, &ret[ 0 ], boost::bind( isGreater< ValueType >, _1, v ) );
//         return ret;
//     } 


#define DEFINE_UNARY_MOD_OPERATOR__( OP, FUNCT ) \
  inline Vector< ValueType > & operator OP##= ( const Vector < ValueType > & v ) { \
        std::transform( data_, data_ + size_, &v[ 0 ], data_, FUNCT() ); return *this; } \
  inline Vector< ValueType > & operator OP##= ( const ValueType & val ) { \
        for ( register Index i = 0; i < size_; i ++ ) data_[ i ] OP##= val; return *this; } \

DEFINE_UNARY_MOD_OPERATOR__( +, PLUS )
DEFINE_UNARY_MOD_OPERATOR__( -, MINUS )
DEFINE_UNARY_MOD_OPERATOR__( /, DIVID )
DEFINE_UNARY_MOD_OPERATOR__( *, MULT )

#undef DEFINE_UNARY_MOD_OPERATOR__

    /*! Negation operator thats return a copy of this with negative values. */
    inline Vector < ValueType > operator - () const { return *this * -1.0; }

    /*! Resize if n differs size() and fill new with val. Old data are preserved. */
    void resize( Index n, ValueType val = 0 ){
        if ( n != size_ ) {
            Vector < ValueType > tmp( *this );
            free_();
            allocate_( n );
            fill( val );
            std::copy( &tmp[ 0 ], &tmp[ min( tmp.size(), n ) ], data_ );
        }
    }

    /*! Fill the whole vector from the pointer of val. CAUTION!! There is no boundary check. val must be properly ( [val, val+this->size() ) )assigned.  */
    template< class V > Vector< ValueType > & fill( V * val ) {
        for ( register Index i = 0; i < size_; i ++ ) data_[ i ] = ValueType(val[i]);
        //std::copy( val, val + size_, data_ );
        return *this;
    }

    /*! Fill the whole vector with val. */
    Vector< ValueType > & fill( const ValueType & val ) { std::fill( data_, data_ + size_, val ); return *this; }

    /*! Fill the whole vector with function expr( i ) */
    template< class Ex > Vector< ValueType > & fill( Expr< Ex > expr ){
        for ( register Index i = 0; i < size_; i ++ ) data_[ i ] = expr( (ValueType)i );
        return *this;
    }

    
//     /*! Fill the whole vector with function expr( i ) */
//     template< class V > void fill( const V & val ){
//         for ( register Index i = 0; i < size_; i ++ ) data_[ i ] = ValueType( val );
//     }

    /*! Fill Vector with 0.0. Don't change size.*/
    void clean( ){ std::memset( data_, '\0', sizeof( ValueType ) * size_ ); }

    /*! Empty the vector. Frees memory and resize to 0.*/
    void clear( ){ free_(); }

    /*! Round all values of this array to a given tolerance. */
    Vector< ValueType > & round( const ValueType & tolerance ){
#ifdef USE_BOOST_BIND
        std::transform( data_, data_ + size_, data_, boost::bind( roundTo< ValueType >, _1, tolerance ) );
#else
        for ( register Index i = 0; i < size_; i ++ ) data_[ i ] = roundTo( data_[ i ], tolerance );
#endif
        return *this;
    }
        
//   ValueType min()

    const VectorIterator< ValueType > & begin() const { return *begin_; }

   // VectorIterator< ValueType > * begin() { return begin_; }

    const VectorIterator< ValueType > & end() const { return *end_; }

    //VectorIterator< ValueType > * end() { return end_; }

    inline bool empty() const { return size() == 0; }

    inline Index size() const { return size_; }

    inline Index nThreads() const { return nThreads_; }

    inline Index singleCalcCount() const { return singleCalcCount_; }

    ValueType * data() { return data_; }

    /*! Save the object to file. Returns true on success and in case of trouble an exception is thrown.
     * The IOFormat flag will be overwritten if the filename have a proper file suffix. Ascii format is forced if \ref VECTORASCSUFFIX is given.
     * Binary format is forced if \ref VECTORBINSUFFIX is set. If no suffix is provided \ref VECTORASCSUFFIX or \ref VECTORBINSUFFIX will be append. \n\n
    Binary format is: \n
    Unsigned int64 [8Byte] - length of the Array [0.. length) \n
    ValueType [sizeof(ValueType)] - 0th value \n
    ... \n
    ValueType [sizeof(ValueType)] - length -1 value \n\n
    Ascii format is a simple list of the values \n\n
    \param filename string of the file name
    \param IOFormat enum, either Ascii for human readable format, or Binary for fast save and load.
    */
    bool save( const std::string & filename, IOFormat format = Ascii ) const {

        if ( filename.rfind( VECTORASCSUFFIX ) != std::string::npos ) format = Ascii;
        else if ( filename.rfind( VECTORBINSUFFIX ) != std::string::npos ) format = Binary;
        std::string fname( filename );

        if ( format == Ascii ){
            if ( fname.rfind( "." ) == std::string::npos ) fname += VECTORASCSUFFIX;

            std::ofstream file; file.open( fname.c_str() );
            if ( !file ) {
                throwError( 1, filename + ": " + strerror( errno ) );
            }

            file.setf( std::ios::scientific, std::ios::floatfield );
            file.precision( 14 );

            for ( Index i = 0, imax = size_; i < imax; i ++ ) file << data_[ i ] << std::endl;
            file.close();
        } else {
            if ( fname.rfind( "." ) == std::string::npos ) fname += VECTORBINSUFFIX;
        // so haett ich das gern //
    //     std::ofstream file( filename.c_str(), std::ofstream::binary );
    //     std::copy( &a[ 0 ], &a[ a.size()-1 ], ostream_iterator< double >( &file ) );
    //     file.close();
            FILE *file; file = fopen( fname.c_str(), "w+b" );
            if ( !file ) {
                throwError( 1, filename + ": " + strerror( errno ) );
            }

            int64 count = (int64)size_;
            Index ret = 0; ret = fwrite( (char*)&count, sizeof( int64 ), 1, file );
            if ( ret == 0) return false;
            for ( Index i = 0; i < size_; i++ ) ret = fwrite( (char*)&data_[ i ], sizeof( ValueType ), 1, file );
            fclose( file );
        }
        return true;
    }

    /*!
     * Load the object from file. Returns true on success and in case of trouble an exception is thrown.
     * The IOFormat flag will be overwritten if the filename have a proper file suffix. Ascii format is forced if \ref VECTORASCSUFFIX is given.
     * Binary format is forced if \ref VECTORBINSUFFIX is set.
     * See Vector< ValueType >::save for file format.
     */
    bool load( const std::string & filename, IOFormat format = Ascii ){

        if ( filename.rfind( VECTORASCSUFFIX ) != std::string::npos ) format = Ascii;
        else if ( filename.rfind( VECTORBINSUFFIX ) != std::string::npos ) format = Binary;

        if ( !fileExist( filename ) ){
            if ( fileExist( filename + VECTORBINSUFFIX ) )
                return this->load( filename + VECTORBINSUFFIX, Binary );
            if ( fileExist( filename + VECTORASCSUFFIX ) )
                return this->load( filename + VECTORASCSUFFIX, Ascii );
        }

        if ( format == Ascii ){
            std::vector < ValueType > tmp;

            std::fstream file; openInFile( filename.c_str(), &file );
            ValueType val; while( file >> val ) {
            // check !!! if ( isinfnan( val ) ) throwLengthError(
                tmp.push_back( val );
            }

    //so haett ich das gern
//     std::ifstream file( filename.c_str() );
//     std::copy(  std::istream_iterator<double>( file ),
//                 std::istream_iterator<double>(),
//                 std::back_inserter( tmp ) );

//std::back_inserter< double > (tmp) );
    //std::copy( file.begin(), file.end(), back_inserter< double >( & tmp[ 0 ] ) );

            this->resize( tmp.size() );
            std::copy( tmp.begin(), tmp.end(), &data_[ 0 ] );
            file.close();

        } else {
            FILE *file;
            file = fopen( filename.c_str(), "r+b" );
            if ( !file ) {
                throwError( 1, filename +  ": " + strerror( errno ) );
            }
            Index ret = 0;
            int64 size; ret = fread( &size, sizeof( int64 ), 1, file );
            if (ret) this->resize( size );

            ret = fread( &data_[ 0 ], sizeof( ValueType ), size, file );
            if (!ret) {
            }
            fclose( file );
        }
        return true;
    }

protected:

    void free_( ){
        //    std::cout << "free: " << begin_ << std::endl;
        size_ = 0;
        if ( data_ )  delete [] data_;
        if ( begin_ ) delete begin_;
        if ( end_ )   delete end_;

        begin_ = NULL;
        end_   = NULL;
        data_  = NULL;
    }

    void allocate_( Index n ){
        size_  = n;
        data_  = new ValueType[ size_ ];
        begin_ = new VectorIterator< ValueType >;
        end_   = new VectorIterator< ValueType >;
        begin_->val_ = data_;
        begin_->maxSize_ = size_;
        end_->val_ = data_ + size_;
        nThreads_ = std::min( (int)ceil( (float)size_ / (float)minSizePerThread ), (int)maxThreads );
        singleCalcCount_ = (Index)ceil( (double)size_ / (double)nThreads_ );
        //std::cout << "alloc: " << begin_ << std::endl;
    }

    void copy_( const Vector< ValueType > & v ){
        if ( v.size() ) {
            resize( v.size() );
            std::copy( &v[ 0 ], &v[ v.size() ], data_ );
        }
    }

    template < class ExprOP > inline void assign_( const ExprOP & v ) {
        if ( v.size() ) {
            //std::cout << "assign_" << v.size() << std::endl;
            resize( v.size() );
            v.assign( *this );
        }
    }

    Index size_;
    ValueType * data_;

    VectorIterator< ValueType > * begin_;
    VectorIterator< ValueType > * end_;

    static const Index minSizePerThread = 10000;
    static const int maxThreads = 8;
    int nThreads_;
    Index singleCalcCount_;
};

template< class ValueType, class Iter > class AssignResult{
public:
    AssignResult( Vector< ValueType > & a, const Iter & result, Index start, Index end ) :
        a_( &a ), iter_( result ), start_( start ), end_( end ){
    }
    void operator()() {
        ValueType * iter = a_->begin().ptr();
        //std::cout << start_ << " " << end_ << std::endl;
        for ( register Index i = start_; i < end_; i++ ) iter[ i ] = iter_[ i ];
    }

    Vector< ValueType > * a_;
    Iter iter_;
    Index start_;
    Index end_;
};

struct BINASSIGN { template < class T > inline T operator()( const T & a, const T & b ) const { return b; } };

template< class ValueType, class Iter > void assignResult( Vector< ValueType > & v, const Iter & result) {
#ifdef EXPRVEC_USE_TEMPORARY_EXPRESSION
  // Make a temporary copy of the iterator.  This is faster on segmented
  // architectures, since all the iterators are in the same segment.
    Iter result2 = result;
#else
    // Otherwise, cast away const (eek!).  No harmful side effects.
    Iter& result2 = (Iter&)result;
#endif


#ifdef EXPRVEC_USE_BOOST_THREAD
//DEPRECATED
    if ( v.nThreads() == 1 ) {
      AssignResult< ValueType, Iter >( v, result, 0, v.size() )();
    } else {
      boost::thread_group threads;
      for ( Index i = 0; i < v.nThreads(); i ++ ){
	Index start = v.singleCalcCount() * i;
	Index end   = v.singleCalcCount() * ( i + 1 );
	if ( i == v.nThreads() -1 ) end = v.size();
	threads.create_thread( AssignResult< ValueType, Iter >( v, result, start, end ) );
      }
      threads.join_all();
    }
#else
#ifdef EXPRVEC_USE_STD_ALGORITHM

    std::transform( v.begin().ptr(), v.end().ptr(), result, v.begin().ptr(), BINASSIGN() );

#else
#ifdef EXPRVEC_USE_INDIRECTION
    ValueType * iter = v.begin().ptr();

    // Inlined expression
    for ( register Index i = v.size(); i--; ) iter[ i ] = result2[ i ];
#else
    ValueType * iter = v.begin().ptr();
    ValueType * end  = v.end().ptr();

    do {
        *iter = *result2;       // Inlined expression
        ++result2;
    } while ( ++iter != end );
#endif
#endif
#endif
}

template< class ValueType, class A > class __VectorExpr {
public:
    __VectorExpr( const A & a ) : iter_( a ) { }

    inline ValueType operator [] ( Index i ) const { return iter_[ i ]; }

    inline ValueType operator * () const { return *iter_; }

    inline void operator ++ () { ++iter_; }

    void assign( Vector< ValueType > & x ) const { assignResult( x, *this ); }

    inline Index size() const { return iter_.size(); }

    A * begin() { return iter_.begin(); }
    A * end() { return iter_.end(); }

private:
    A iter_;
};

template< class ValueType, class A, class Op > class __VectorUnaryExprOp {
public:
    __VectorUnaryExprOp( const A & a ) : iter_( a ) { }

    inline ValueType operator [] ( Index i ) const { return Op()( iter_[ i ] ); }

    inline ValueType operator * () const { return Op()( *iter_ ); }

    inline void operator ++ () { ++iter_;  }

    inline Index size() const { return iter_.size(); }

private:
    A iter_;
};

template< class ValueType, class A, class B, class Op > class __VectorBinaryExprOp {
public:
  __VectorBinaryExprOp( const A & a, const B & b) : iter1_( a ), iter2_( b ) { }

  inline ValueType operator [] ( Index i ) const { return Op()( iter1_[ i ], iter2_[ i ] ); }

  inline ValueType operator * () const { return Op()( *iter1_, *iter2_ ); }

  inline void operator ++ () { ++iter1_; ++iter2_; }

  inline Index size() const { return iter2_.size(); }

private:
  A iter1_;
  B iter2_;
};

template< class ValueType, class A, class Op > class __VectorValExprOp {
public:
  __VectorValExprOp( const A & a, const ValueType & val ) : iter_( a ), val_( val ) { }

  inline ValueType operator [] ( Index i ) const { return Op()( iter_[ i ], val_ ); }

  inline ValueType operator * () const { return Op()( *iter_, val_ ); }

  inline void operator ++ () { ++iter_; }

  inline Index size() const { return iter_.size(); }

private:
  A iter_;
  ValueType val_;
};

template< class ValueType, class A, class Op > class __ValVectorExprOp {
public:
  __ValVectorExprOp( const ValueType & val, const A & a ) : iter_( a ), val_( val ) { }

  inline ValueType operator [] ( Index i ) const { return Op()( val_, iter_[ i ] ); }

  inline ValueType operator * () const { return Op()( val_, *iter_ ); }

  inline void operator ++ () { ++iter_; }

  inline Index size() const { return iter_.size(); }

private:
  A iter_;
  ValueType val_;
};

#define DEFINE_UNARY_EXPR_OPERATOR__( OP, FUNCT )\
\
template < class T > \
__VectorExpr< T, __VectorUnaryExprOp< T, VectorIterator< T >, FUNCT > > \
OP( const Vector< T > & a ){ \
    typedef __VectorUnaryExprOp< T, VectorIterator< T >, FUNCT > ExprT; \
    return __VectorExpr< T, ExprT >( ExprT( a.begin() ) ); } \
\
template < class T, class A > \
__VectorExpr< T, __VectorUnaryExprOp< T, __VectorExpr< T, A >, FUNCT > > \
OP( const __VectorExpr< T, A > & a ){ \
    typedef __VectorUnaryExprOp< T, __VectorExpr< T, A >, FUNCT > ExprT; \
    return __VectorExpr< T, ExprT >( ExprT( a ) ); \
} \

DEFINE_UNARY_EXPR_OPERATOR__( abs,   ABS_ )
DEFINE_UNARY_EXPR_OPERATOR__( acot,  ACOT )
DEFINE_UNARY_EXPR_OPERATOR__( atan,  ATAN )
DEFINE_UNARY_EXPR_OPERATOR__( cos,   COS )
DEFINE_UNARY_EXPR_OPERATOR__( cot,   COT )
DEFINE_UNARY_EXPR_OPERATOR__( exp,   EXP )
DEFINE_UNARY_EXPR_OPERATOR__( exp10, EXP10 )
DEFINE_UNARY_EXPR_OPERATOR__( fabs,  ABS_ )
DEFINE_UNARY_EXPR_OPERATOR__( log,   LOG )
DEFINE_UNARY_EXPR_OPERATOR__( log10, LOG10 )
DEFINE_UNARY_EXPR_OPERATOR__( sign, SIGN )
DEFINE_UNARY_EXPR_OPERATOR__( sin,   SIN )
DEFINE_UNARY_EXPR_OPERATOR__( sqrt,  SQRT )
DEFINE_UNARY_EXPR_OPERATOR__( square, SQR )
DEFINE_UNARY_EXPR_OPERATOR__( tan,   TAN )
DEFINE_UNARY_EXPR_OPERATOR__( tanh,  TANH )

#undef DEFINE_UNARY_EXPR_OPERATOR__

#define DEFINE_EXPR_OPERATOR__( OP, FUNCT )				\
template < class T >							\
__VectorExpr< T, __VectorBinaryExprOp< T, VectorIterator< T >, VectorIterator< T >, FUNCT > > \
operator OP (const Vector< T > & a, const Vector< T > & b ){		\
    typedef __VectorBinaryExprOp< T, VectorIterator< T >, VectorIterator< T >, FUNCT > ExprT; \
    return __VectorExpr< T, ExprT >( ExprT( a.begin(), b.begin() ) );	\
}									\
                                                                        \
template < class T >							\
__VectorExpr< T, __VectorValExprOp< T, VectorIterator< T >, FUNCT > >	\
operator OP ( const Vector< T > & a, const T & val ){			\
  typedef __VectorValExprOp< T, VectorIterator< T >, FUNCT > ExprT;        \
  return __VectorExpr< T, ExprT >( ExprT( a.begin(), val ) );		\
}                                                                       \
                                                                        \
template < class T >							\
__VectorExpr< T, __ValVectorExprOp< T, VectorIterator< T >, FUNCT > >	\
operator OP ( const T & val, const Vector< T > & a ){			\
  typedef __ValVectorExprOp< T, VectorIterator< T >, FUNCT > ExprT;        \
  return __VectorExpr< T, ExprT >( ExprT( val, a.begin() ) );		\
}									\
    									\
template< class T, class A >						\
__VectorExpr< T, __VectorBinaryExprOp< T, __VectorExpr< T, A >, VectorIterator< T >, FUNCT > > \
operator OP ( const __VectorExpr< T, A > & a, const Vector< T > & b ){	\
  typedef __VectorBinaryExprOp< T, __VectorExpr< T, A >, VectorIterator< T >, FUNCT > ExprT; \
  return __VectorExpr< T, ExprT >( ExprT( a, b.begin() ) );		\
}									\
    									\
template< class T, class A >					\
__VectorExpr< T, __VectorBinaryExprOp< T, VectorIterator< T >, __VectorExpr< T, A >, FUNCT > > \
operator OP ( const Vector< T > & a, const __VectorExpr< T, A > & b ){	\
  typedef __VectorBinaryExprOp< T, VectorIterator< T >, __VectorExpr< T, A >, FUNCT > ExprT; \
  return __VectorExpr< T, ExprT >( ExprT( a.begin(), b ) );		\
}									\
									\
template< class T, class A >					\
__VectorExpr< T, __VectorValExprOp< T, __VectorExpr< T, A >, FUNCT > >	\
operator OP ( const __VectorExpr< T, A > & a, const T & val ){	\
  typedef __VectorValExprOp< T, __VectorExpr< T, A >, FUNCT > ExprT;	\
  return __VectorExpr< T, ExprT >( ExprT( a, val ) );		\
}								\
        \
template< class T, class A >				\
__VectorExpr< T, __ValVectorExprOp< T, __VectorExpr< T, A >, FUNCT > >	\
operator OP ( const T & val, const __VectorExpr< T, A > & a ){	\
  typedef __ValVectorExprOp< T, __VectorExpr< T, A >, FUNCT > ExprT;	\
  return __VectorExpr< T, ExprT >( ExprT( val, a ) );		\
}								\
                                                                \
template< class T, class A, class B >				\
__VectorExpr< T, __VectorBinaryExprOp< T, __VectorExpr< T, A >, __VectorExpr< T, B >, FUNCT > > \
operator OP ( const __VectorExpr< T, A > & a, const __VectorExpr< T, B > & b ){ \
  typedef __VectorBinaryExprOp< T, __VectorExpr< T, A >, __VectorExpr< T, B >, FUNCT > ExprT; \
  return __VectorExpr< T, ExprT >( ExprT( a, b ) );			\
}									\
\
template< class T, class T2, class A >					\
        __VectorExpr< T, __VectorValExprOp< T, __VectorExpr< T, A >, FUNCT > >	\
        operator OP ( const __VectorExpr< T, A > & a, const T2 & val ){	\
        typedef __VectorValExprOp< T, __VectorExpr< T, A >, FUNCT > ExprT;	\
        return __VectorExpr< T, ExprT >( ExprT( a, (T)val ) );		\
}								\
        \
template< class T, class T2, class A >				\
        __VectorExpr< T, __ValVectorExprOp< T, __VectorExpr< T, A >, FUNCT > >	\
        operator OP ( const T2 & val, const __VectorExpr< T, A > & a ){	\
        typedef __ValVectorExprOp< T, __VectorExpr< T, A >, FUNCT > ExprT;	\
        return __VectorExpr< T, ExprT >( ExprT( (T)val, a ) );		\
}								\
        \
template < class T, class T2 >					        \
        __VectorExpr< T, __ValVectorExprOp< T, VectorIterator< T >, FUNCT > >	\
        operator OP ( const T2 & val, const Vector< T > & a ){			\
        typedef __ValVectorExprOp< T, VectorIterator< T >, FUNCT > ExprT;        \
        return __VectorExpr< T, ExprT >( ExprT( (T)val, a.begin() ) );		\
}									\
        \
template < class T, class T2 >						\
        __VectorExpr< T, __VectorValExprOp< T, VectorIterator< T >, FUNCT > >	\
        operator OP ( const Vector< T > & a, const T2 & val ){			\
        typedef __VectorValExprOp< T, VectorIterator< T >, FUNCT > ExprT;      \
        return __VectorExpr< T, ExprT >( ExprT( a.begin(), (T)val ) );	\
}                                                                       \
        \

DEFINE_EXPR_OPERATOR__( +, PLUS  )
DEFINE_EXPR_OPERATOR__( -, MINUS )
DEFINE_EXPR_OPERATOR__( *, MULT )
DEFINE_EXPR_OPERATOR__( /, DIVID )

#undef DEFINE_EXPR_OPERATOR__

//********************************************************************************
//** define some utility functions

// inline bool operator < (const GIMLI::Vector<double>&a, const GIMLI::Vector<double> &b) {
//     return false;
// }

template < class ValueType >
bool operator == ( const Vector< ValueType > & v1, const Vector< ValueType > & v2 ){
    if ( v1.size() != v2.size() ) return false;
    for ( Index i = 0; i < v1.size(); i ++ ){
        if ( !isEqual( v1[ i ], v2[ i ] ) ) return false;
    }
    return true;
}

template < class ValueType, class A >
bool operator == ( const Vector< ValueType > & v1, const __VectorExpr< ValueType, A > & v2 ){
    return v1 == Vector< ValueType >( v2 );
}

/*! Return true if at least one value is greater than TOLERANCE. */
template < class ValueType >
bool nonZero(const Vector< ValueType > & v){
    return (min(abs(v)) > TOLERANCE);
}



template < class ValueType >
bool operator != ( const Vector< ValueType > & v1, const Vector< ValueType > & v2 ){
    return !(v1==v2);
}

/*! Find function. Return index vector of true values */
inline IndexArray find( const BVector & v ){
    IndexArray idx;
    idx.reserve( v.size() );
    for ( Index i = 0; i < v.size(); i ++ ){
        if ( v[ i ] ) idx.push_back( i );
    }
    return idx;
}

/*! Refactor with expression templates */
inline BVector operator ~ ( const BVector & a ){
    BVector ret( a.size() );
    for ( Index i = 0; i < ret.size(); i ++ ) ret[ i ] = !a[ i ];
    return ret;
}

/*! Refactor with expression templates */
inline BVector operator & ( const BVector & a, const BVector & b ){
    BVector ret( a.size() );
    for ( Index i = 0; i < ret.size(); i ++ ) ret[ i ] = a[ i ] && b[ i ];
    return ret;
}

/*! Refactor with expression templates */
inline BVector operator | ( const BVector & a, const BVector & b ){
    BVector ret( a.size() );
    for ( Index i = 0; i < ret.size(); i ++ ) ret[ i ] = a[ i ] || b[ i ];
    return ret;
}

inline RVector operator * ( const BVector & a, const RVector & b ){
    RVector ret( a.size() );
    for ( Index i = 0; i < ret.size(); i ++ ) ret[ i ] = a[ i ] * b[ i ];
    return ret;
}

inline RVector operator * ( const RVector & a, const BVector & b ){
    RVector ret( a.size() );
    for ( Index i = 0; i < ret.size(); i ++ ) ret[ i ] = a[ i ] * b[ i ];
    return ret;
}

// /*! Refactor for speed */
// template < class ValueType, class A > BVector
// operator == ( const __VectorExpr< ValueType, A > & vec, const ValueType & v ){
//     return Vector< ValueType >( vec ) == v;
// }

#define DEFINE_COMPARE_OPERATOR__( OP ) \
template < class ValueType, class A > BVector \
operator OP ( const __VectorExpr< ValueType, A > & vec, const ValueType & v ){ \
    BVector ret( vec.size(), 0 ); \
    for ( Index i = 0; i < ret.size(); i ++ ) ret[ i ] = vec[ i ] OP v; \
    return ret;\
} \
template < class ValueType > BVector \
operator OP ( const std::vector < ValueType > & vec, const ValueType & v ){ \
    BVector ret( vec.size(), 0 ); \
    for ( Index i = 0; i < ret.size(); i ++ ) ret[ i ] = vec[ i ] OP v; \
    return ret;\
} \

DEFINE_COMPARE_OPERATOR__( < )
DEFINE_COMPARE_OPERATOR__( <= )
DEFINE_COMPARE_OPERATOR__( >= )
DEFINE_COMPARE_OPERATOR__( == )
DEFINE_COMPARE_OPERATOR__( != )
DEFINE_COMPARE_OPERATOR__( > )

#undef DEFINE_COMPARE_OPERATOR__

#define DEFINE_UNARY_COMPARE_OPERATOR__( OP, FUNCT ) \
template < class ValueType, class A > \
BVector OP( const __VectorExpr< ValueType, A > & vec ){ \
    BVector ret( vec.size(), 0 ); \
    for ( Index i = 0; i < ret.size(); i ++ ) ret[ i ] = FUNCT()( vec[ i ] ); \
    return ret; \
}\
template < class ValueType > \
BVector OP ( const Vector< ValueType > & vec ){ \
    BVector ret( vec.size(), 0 ); \
    std::transform( vec.begin().ptr(), vec.end().ptr(), ret.begin().ptr(), FUNCT() ); \
    return ret; \
} \

DEFINE_UNARY_COMPARE_OPERATOR__( isInf, ISINF)
DEFINE_UNARY_COMPARE_OPERATOR__( isNaN, ISNAN)
DEFINE_UNARY_COMPARE_OPERATOR__( isInfNaN, ISINFNAN)

#undef DEFINE_UNARY_COMPARE_OPERATOR__

template < class T > Vector < T > cat( const Vector< T > & a, const Vector< T > & b ){
    Vector < T > c ( a.size() + b.size() );
    std::copy( &a[ 0 ], &a[ a.size() ], &c[ 0 ] );
    std::copy( &b[ 0 ], &b[ b.size() ], &c[ a.size() ] );
    return c;
}

template < class T, class A > T sum( const __VectorExpr< T, A > & a ){
    //std::cout << "sum(vectorExpr)" << std::endl;
    T tmp( 0.0 );
    for ( register Index i = 0, imax = a.size(); i < imax; i++ ) tmp += a[ i ];
    return tmp;

//     T tmp( 0.0 );
//     __VectorExpr< T, A > al = a;
//     for ( register Index i = 0; i < a.size(); i++, ++al ) {
//         tmp += *al;
//     }
//     return tmp;

//     T tmp( 0.0 );
//     __VectorExpr< T, A > al = a;
//     for (; al != a.end(); ++al){
//         tmp += *al;
//     }
//     return tmp;


//return std::accumulate( a[0], a[a.size()], T() );
}

template < class T > T sum( const Vector < T > & v ){
     //std::cout << "sum(vector)" << std::endl;
//     std::cout << *v.begin() << " "  << *v.end() << std::endl;

    return std::accumulate( v.begin(), v.end(), T(0) );
    //return std::accumulate( v.begin(), v.end(), (T)0.0 );
}

template < class T, class A > T min( const __VectorExpr< T, A > & a ){ return min( Vector< T >( a ) ); }
template < class T, class A > T max( const __VectorExpr< T, A > & a ){ return max( Vector< T >( a ) ); }

template < class T > T min( const Vector < T > & v ){
    return *std::min_element( &v[ 0 ], &v[ 0 ] + v.size() );
}

template < class T > T max( const Vector < T > & v ){
    return *std::max_element( &v[ 0 ], &v[ 0 ] + v.size() );
}

template < class ValueType > bool haveInfNaN( const Vector < ValueType > & v ){
    for ( VectorIterator < ValueType > it = v.begin(); it != v.end(); ++it ){
        if ( isInfNaN( *it ) ) return true;
    }
    return false;
}

template < class ValueType > Vector < ValueType > fixZero( const Vector < ValueType > & v, const ValueType tol = TOLERANCE ){
    Vector < ValueType > ret( v );
    for ( VectorIterator < ValueType > it = ret.begin(); it != ret.end(); ++it ){
        if ( ::fabs( *it ) < TOLERANCE ) *it = tol;
    }
    return ret;
}
 
template < class ValueType > Vector < ValueType > round( const Vector < ValueType > & v, ValueType tol ){
    return Vector< ValueType >( v ).round( tol );
}

template < class T >
Vector < T > fliplr( const Vector < T > & v ){
    Vector < T > n( v.size() );
    for ( Index i = 0; i < v.size(); i ++ ) n[ i ] = v[ v.size() - 1 - i ];
    return n;
}

template < class T, class A, class T2 > Vector < T > pow( const __VectorExpr< T, A > & a, T2 power ){
    return pow( Vector< T >( a ), power );
}

template < class T, class T2 > Vector < T > pow( const Vector < T > & v, T2 npower ){
    Vector < T > r( v.size() );
    for ( Index i = 0; i < v.size(); i ++ ) r[ i ] = std::pow( v[ i ], T( npower ) );
    return r;
}

template < class T > Vector< T > sort( const Vector < T > & a ){
    std::vector < T > tmp( a.size(), 0.0 ) ;
    for ( Index i = 0; i < a.size(); i ++ ) tmp[ i ] = a[ i ];
    std::sort( tmp.begin(), tmp.end() );

    Vector < T > ret( tmp  );
    return ret;
//     Vector < T > t( a );
//     std::sort( t.begin(), t.end() );
//     return t;
}

/*! Returning a copy of the vector and replacing all consecutive occurrences of a value by a single instance of that value. e.g. [0 1 1 2 1 1] -> [0 1 2 1 ]. To remove all double values from the vector use an additionally sorting. e.g. unique( sort( v ) ) gets you [ 0 1 2 ]. */
template < class T > Vector< T > unique( const Vector < T > & a ){
    std::vector < T > tmp( a.size() ), u;
    for ( Index i = 0; i < a.size(); i ++ ) tmp[ i ] = a[ i ];
    std::unique_copy( tmp.begin(), tmp.end(), back_inserter( u ) );

    Vector < T > ret( u  );
    return ret;
}

// template < template < class T > class Vec, class T >
// std::ostream & operator << ( std::ostream & str, const Vec < T > & vec ){
//     for ( Index i = 0; i < vec.size(); i ++ ) str << vec[ i ] << " ";
//     return str;
// }

template < class T >
std::ostream & operator << ( std::ostream & str, const std::vector < T > & vec ){
    for ( Index i = 0; i < vec.size(); i ++ ) str << vec[ i ] << " ";
    return str;
}

template < class T >
std::ostream & operator << ( std::ostream & str, const Vector < T > & vec ){
    for ( Index i = 0; i < vec.size(); i ++ ) str << vec[ i ] << " ";
    return str;
}

/*!
Return a RVector with increasing values of size(n+1) filled with : 0, first, ... ,last
*/
template < class ValueType >
Vector< ValueType > increasingRange( const ValueType & first, const ValueType & last, Index n ){
    if ( sign( first ) != sign( last ) ){
        throwError( 1, "cant increase range from [0 " + toStr( first) + " to " + toStr( last ) + " ]" );
    }
    
    Placeholder x__;
    RVector y( n + 1 ); y.fill( x__ );

    ValueType dy = ( last - first * n ) / ( sum( y ) - ValueType( n ) );

    if ( dy < 0.0 ){
        std::cout << "decreasing number of layers: " << n << " " << dy << std::endl;
        return increasingRange( first, last, n-1);
    }

    ValueType yval = 0.0;
    for ( Index i = 0; i < n; i ++ ){
        yval = yval + first + dy * i;
        y[ i + 1 ] = yval;
    }
    //y = fliplr( y ) * -1.0 //#+ g.max( y )
    return y;
    //return fliplr( y ) * -1.0;
}

template < class ValueType >
Vector < std::complex < ValueType > > toComplex( const Vector < ValueType > & re,
                                                 const Vector < ValueType > & im ){
    Vector < std::complex < ValueType > > cv( re.size() );
    for ( Index i = 0; i < cv.size(); i ++ ) cv[ i ] = std::complex < ValueType >( re[ i ], im[ i ] );
    return cv;
}

inline CVector toComplex( double re, const RVector & im ){
     return toComplex( RVector( im.size(), re ), im );
}

// template < class ValueType >
// Vector < std::complex < ValueType > > toComplex( ValueType re, const Vector < ValueType > & im ){
//     return toComplex( Vector < ValueType > ( im.size(), re ), im );
// }

template < class ValueType >
Vector < std::complex < ValueType > > toComplex( const Vector < ValueType > & re, ValueType im = 0 ){
    return toComplex( re, Vector < ValueType > ( re.size(), im ) );
}

template < class ValueType >
Vector < std::complex < ValueType > > operator * ( const Vector < std::complex< ValueType > > & cv,
                                                   const Vector < ValueType > & v ){
    return cv * toComplex( v );
}

template < class ValueType >
Vector < std::complex < ValueType > > operator * ( const Vector < ValueType > & v,
                                                   const Vector < std::complex< ValueType > > & cv ){
    return cv * toComplex( v );
}

template < class ValueType, class A >
Vector < ValueType > real( const __VectorExpr< std::complex< ValueType >, A > & a ){
    return real( Vector < std::complex< ValueType > >( a ) );
}

template < class ValueType >
Vector < ValueType > real( const Vector < std::complex< ValueType > > & cv ){
    Vector < ValueType > v( cv.size() );
    for ( Index i = 0; i < cv.size(); i ++ ) v[ i ] = cv[ i ].real();
    return v;
}

template < class ValueType, class A >
Vector < ValueType > imag( const __VectorExpr< std::complex< ValueType >, A > & a ){
    return imag( Vector < std::complex< ValueType > >( a ) );
}

template < class ValueType >
Vector < ValueType > imag( const Vector < std::complex< ValueType > > & cv ){
    Vector < ValueType > v( cv.size() );
    for ( Index i = 0; i < cv.size(); i ++ ) v[ i ] = cv[ i ].imag();
    return v;
}

template < class ValueType, class A >
Vector < ValueType > angle( const __VectorExpr< std::complex< ValueType >, A > & a ){
    return angle( Vector < std::complex< ValueType > >( a ) );
}

template < class ValueType >
Vector < ValueType > angle( const Vector < std::complex< ValueType > > & cv ){
    return imag( log( cv ) );
}

template < class ValueType, class A >
Vector < ValueType > abs( const __VectorExpr< std::complex< ValueType >, A > & a ){
    return abs( Vector < std::complex< ValueType > >( a ) );
}

template < class ValueType >
Vector < ValueType > abs( const Vector < std::complex< ValueType > > & cv ){
    return sqrt( real( cv * conj( cv ) ) );
}

template < class ValueType, class A >
Vector < std::complex< ValueType > > conj( const __VectorExpr< std::complex< ValueType >, A > & a ){
    return conj( Vector < std::complex< ValueType > >( a ) );
}

template < class ValueType >
Vector < std::complex< ValueType > > conj( const Vector < std::complex< ValueType > > & cv ){
    Vector < std::complex< ValueType > > v( cv.size() );
    for ( Index i = 0; i < cv.size(); i ++ ) v[ i ] = Complex( cv[ i ].real(), -cv[ i ].imag() );
    return v;
}

#define DEFINE_SCALAR_COMPLEX_BINARY_OPERATOR( OP ) \
template <class T, class U > \
inline std::complex< T > operator OP ( const std::complex< T > & lhs, const U & rhs ) { \
    std::complex< T > ret; \
    return ret OP##= rhs; \
} \

DEFINE_SCALAR_COMPLEX_BINARY_OPERATOR( * )
DEFINE_SCALAR_COMPLEX_BINARY_OPERATOR( / )
DEFINE_SCALAR_COMPLEX_BINARY_OPERATOR( - )
DEFINE_SCALAR_COMPLEX_BINARY_OPERATOR( + )

#undef DEFINE_SCALAR_COMPLEX_BINARY_OPERATOR

template < class ValueType >
bool save( const Vector< ValueType > & a, const std::string & filename, IOFormat format = Ascii ){
    return saveVec( a, filename, format );
}

template < class ValueType >
bool load( Vector< ValueType > & a, const std::string & filename, IOFormat format = Ascii,
                                    bool verbose = true ){
    return loadVec( a, filename, format, verbose );
}

/*!
 Save vector to file. See Vector< ValueType >::save.
*/
template < class ValueType >
    bool saveVec( const Vector< ValueType > & a, const std::string & filename, IOFormat format = Ascii ){ return a.save( filename, format ); }

/*!
 Load vector from file. See Vector< ValueType >::load.
*/
template < class ValueType >
    bool loadVec(Vector < ValueType > & a, 
                 const std::string & filename, 
                 IOFormat format = Ascii ){ return a.load( filename, format ); }



} // namespace GIMLI

#endif
