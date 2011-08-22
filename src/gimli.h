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

#ifndef _GIMLI_GIMLI__H
#define _GIMLI_GIMLI__H

// we don't want the Old compatibility names for C types: typedef unsigned int uint
#undef __USE_MISC


#ifdef HAVE_CONFIG_H
    #include <config.h>
#else
    #ifndef PACKAGE_NAME
        #define PACKAGE_NAME "gimli"
        #define PACKAGE_VERSION "0.7.0-win"
        #define PACKAGE_BUGREPORT "carsten@resistivity.net"
    #endif // PACKAGE_NAME
#endif

#define PACKAGE_AUTHORS "carsten@resistivity.net, thomas@resistivity.net"
#include "platform.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <stdint.h>
#include <complex>
#include <limits>
#include <algorithm>

#include "exitcodes.h"

//! GIMLi main namespace for the Geophyiscal Inversion and Modelling Library
namespace GIMLI{

typedef size_t  uint;
typedef uint8_t  uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

typedef int8_t  int8;
typedef int16_t int16;
typedef int32_t int32;
typedef int64_t int64;

#ifndef __ASSERT_FUNCTION
#define __ASSERT_FUNCTION "__ASSERT_FUNCTION"
#endif

#define WHERE std::string( __FILE__ ) + ": " + GIMLI::toStr( __LINE__ )+ "\t"
#define WHERE_AM_I std::string( WHERE + "\t" + std::string( __ASSERT_FUNCTION ) + " " )
#define THROW_TO_IMPL { std::stringstream str; str << WHERE_AM_I << " not yet implemented\n " << GIMLI::versionStr() << "\nPlease send the messages above, the commandline and all necessary data to the author." << std::endl;\
                         throwToImplement( str.str() ); }
#define CERR_TO_IMPL std::cerr << WHERE_AM_I << " not yet implemented\n " << GIMLI::versionStr() << "\nPlease send the messages above, the commandline and all necessary data to the author." << std::endl;
#define DEPRECATED std::cerr << WHERE_AM_I << " is deprecated " << std::endl;

#define TOLERANCE 1e-12
#define TOUCH_TOLERANCE 1e-12
#define MAX_DOUBLE std::numeric_limits<double>::max()
#define MAX_INT std::numeric_limits< int >::max()

#define MESHBINSUFFIX   ".bms"
#define MATRIXBINSUFFIX ".bmat"
#define MATRIXASCSUFFIX ".matrix"
#define VECTORBINSUFFIX ".bvec"
#define VECTORASCSUFFIX ".vector"
#define NOT_DEFINED "notDefined"

static const int MARKER_BOUND_HOMOGEN_NEUMANN = -1;
static const int MARKER_BOUND_MIXED = -2;
static const int MARKER_BOUND_HOMOGEN_DIRICHLET = -3;
static const int MARKER_BOUND_DIRICHLET = -4;
static const int MARKER_BOUND_ELECTRODE = -10000;
static const int MARKER_CELL_PARAMETER = 2;

static const uint8 MESH_BASEENTITY_RTTI      = 000;
static const uint8 MESH_MESHENTITY_RTTI      = 001;
static const uint8 MESH_NODE_RTTI            = 010;
static const uint8 MESH_BOUNDARY_RTTI        = 020;
static const uint8 MESH_BOUNDARY_NODE_RTTI   = 021;
static const uint8 MESH_EDGE_RTTI            = 022;
static const uint8 MESH_EDGE3_RTTI           = 023;
static const uint8 MESH_TRIANGLEFACE_RTTI    = 024;
static const uint8 MESH_TRIANGLEFACE6_RTTI   = 025;
static const uint8 MESH_QUADRANGLEFACE_RTTI  = 026;
static const uint8 MESH_QUADRANGLEFACE8_RTTI = 027;
static const uint8 MESH_CELL_RTTI            = 030;
static const uint8 MESH_EDGE_CELL_RTTI       = 031;
static const uint8 MESH_EDGE3_CELL_RTTI      = 032;
static const uint8 MESH_TRIANGLE_RTTI        = 033;
static const uint8 MESH_TRIANGLE6_RTTI       = 034;
static const uint8 MESH_QUADRANGLE_RTTI      = 035;
static const uint8 MESH_QUADRANGLE8_RTTI     = 036;
static const uint8 MESH_TETRAHEDRON_RTTI     = 041;
static const uint8 MESH_TETRAHEDRON10_RTTI   = 042;
static const uint8 MESH_HEXAHEDRON_RTTI      = 043;
static const uint8 MESH_HEXAHEDRON20_RTTI    = 044;

static const uint8 MESH_SHAPE_NODE_RTTI        = 210;
static const uint8 MESH_SHAPE_EDGE_RTTI        = 211;
static const uint8 MESH_SHAPE_TRIANGLE_RTTI    = 221;
static const uint8 MESH_SHAPE_QUADRANGLE_RTTI  = 222;
static const uint8 MESH_SHAPE_TETRAHEDRON_RTTI = 231;
static const uint8 MESH_SHAPE_HEXAHEDRON_RTTI  = 232;

/*! Flag load/save Ascii or binary */
enum IOFormat{ Ascii, Binary };

//** start forward declaration
class Boundary;
class Cell;
class DataContainer;
class DataMap;
class Line;
class Mesh;
class MeshEntity;
class ModellingBase;
class Node;
class Plane;
class Region;
class RegionManager;
class Shape;

template < class ValueType > class Pos;
typedef Pos< int >          IntPos;
typedef Pos< double >       RVector3;

// DEPRECATED ??
// template < class ValueType > class STLMatrix;
// typedef STLMatrix< double > RSTLMatrix;

template < class ValueType >        class SparseMatrix;
typedef SparseMatrix< int >         ISparseMatrix;
typedef SparseMatrix< double >      DSparseMatrix;

template< class ValueType, class IndexType > class SparseMapMatrix;
typedef SparseMapMatrix< int, size_t >     ISparseMapMatrix;
typedef SparseMapMatrix< double, size_t >  DSparseMapMatrix;

template < class ValueType > class Vector;
//template <> class Vector< double >;

typedef std::complex < double > Complex;
typedef Vector< double > RVector;
typedef Vector< Complex > CVector;
typedef Vector< int > BVector;
typedef Vector< long > LVector;
//#typedef Vector< unsigned char > BVector;

template < class ValueType > class Matrix;

typedef Matrix < double > RMatrix;
//typedef Matrix < Complex > CMatrix;
// typedef std::vector < RVector > RMatrix;
// typedef std::vector < CVector > CMatrix;

template < class ValueType > class ElementMatrix;

template < class Vec > class Trans;

//** end forward declaration

DLLEXPORT extern bool __SAVE_PYTHON_GIL__;

/*! */
inline void savePythonGIL( bool s ){
    __SAVE_PYTHON_GIL__ = s;
}


class PythonGILSave {
public:
    PythonGILSave( ): saved_( false ) { save(); }
    ~PythonGILSave() { restore(); }
    void save() { if( !saved_ ) {
#ifdef PYGIMLI
//#warning  "save_ = PyEval_SaveThread();"
        if ( __SAVE_PYTHON_GIL__ ) save_ =  PyEval_SaveThread();
#endif
        saved_ = true; } }
    void restore() { if ( saved_ ) {
#ifdef PYGIMLI
        if ( __SAVE_PYTHON_GIL__ ) PyEval_RestoreThread( save_ );
#endif
        saved_ = false; } }
private:
    bool saved_;
#ifdef PYGIMLI
    PyThreadState *save_;
#endif
};

#define ALLOW_PYTHON_THREADS PythonGILSave __pygil_t__;
#define RESTORE_PYTHON_THREADS __pygil_t__.restore();
#define SAVE_PYTHON_THREADS __pygil_t__.save();
#ifndef Py_BEGIN_ALLOW_THREADS
    #define Py_BEGIN_ALLOW_THREADS
    #define Py_END_ALLOW_THREADS
#endif

//DLLEXPORT std::string version();

inline std::string versionStr(){
  std::string vers( (std::string)(PACKAGE_NAME) + "-" + PACKAGE_VERSION );
  return vers;
}

DLLEXPORT std::string authors();

template < class T > inline T min( const T & a, const T & b ){ return std::min( a, b ); }
template < class T > inline T max( const T & a, const T & b ){ return std::max( a, b ); }

DLLEXPORT int openFile( const std::string & fname, std::fstream * file,
                        std::_Ios_Openmode farg, bool terminate );

inline int openInFileTerm( const std::string & fname, std::fstream * file ){
    return openFile( fname, file, std::ios::in, true);
}
inline int openInFile( const std::string & fname, std::fstream * file, bool terminate = true ){
    return openFile( fname, file, std::ios::in, terminate );
}
inline int openOutFile( const std::string & fname, std::fstream * file, bool terminate = true ){
    return openFile( fname, file, std::ios::out, terminate );
}

DLLEXPORT bool fileExist( const std::string & filename );
DLLEXPORT size_t countColumnsInFile( const std::string & fname, size_t & columnCount );
DLLEXPORT size_t countColumnsInFile( const std::string & fname );
DLLEXPORT size_t countRowsInFile( const std::string & fname );
DLLEXPORT size_t fileLength( std::fstream & file );
DLLEXPORT std::vector < std::string > getRowSubstrings( std::fstream & file, char comment = '#' );

inline std::vector < std::string > getRow( std::fstream & file, char comment = '#' ){
    return getRowSubstrings( file, comment );
}

DLLEXPORT std::vector < std::string > getNonEmptyRow( std::fstream & file, char comment = '#' );
DLLEXPORT std::vector < std::string > getCommentLine( std::fstream & file, char comment = '#' );

DLLEXPORT std::vector < std::string > getSubstrings( const std::string & str );
DLLEXPORT std::vector < std::string > split( const std::string & str, char delimiter = ':' );


DLLEXPORT std::map < float, float > loadFloatMap( const std::string & filename );
DLLEXPORT std::map < int, int > loadIntMap( const std::string & filename );

inline int       toInt( const std::string & str ){ return atoi( str.c_str() ); }
inline float   toFloat( const std::string & str ){ return strtof( str.c_str(), NULL ); }
inline double toDouble( const std::string & str ){ return strtod( str.c_str(), NULL ); }

//! General template for conversion to string, shoul supersede all sprintf etc.
template< typename T > inline std::string str( const T & value ){
  std::ostringstream streamOut;
  streamOut << value;
  return streamOut.str();
}

//! Deprecated! use str() instead, General template for conversion to string, shoul supersede all sprintf etc.
template< typename T > inline std::string toStr( const T & value ){
    return str( value );
}

inline std::string strReplaceBlankWithUnderscore( const std::string & str ) {
    std::string res( str );
    for ( uint i = 0; i < res.length(); i ++ ) if ( res[ i ] == ' ' ) res[ i ] = '_';
    return res;
}

inline std::string lower( const std::string & str ){
    std::string lo( str );
    std::transform( lo.begin(), lo.end(), lo.begin(), tolower );
    return lo;
}

template < typename T > inline void swapVal( T & a, T & m ){
    T tmp( a ); a = m; m = tmp;
}

template < typename T > void test( T a, T b, std::vector < bool > & success ){
    success.push_back( a == b );
    if ( !success.back() ){
        std::cout << "test " << success.size() << " ist: " << a << " soll: " << b << std::endl;
    }
}

template < typename T > bool test( T a, T b, bool verbose = false ){
    if ( verbose ){
        std::cout << "ist: " << a << " soll: " << b << std::endl;
    }
    return a == b;
}

/*! General template for deleting an object. This is not exception-safe unless you use some kind of smart pointer.\n
Example: Delete all objects in a container.
vector < ptr * > vP;
// ... // something that fills vP with the new operator.
for_each( vP.begin(), vP.end(), deletePtr() );
*/
struct DLLEXPORT deletePtr{
    template < typename T > void operator()( T * p ) { delete p; }
};
struct DLLEXPORT cerrPtr{
  template < typename T > void operator() ( const T * p ) const { std::cerr << p << " " << std::endl; }
};
struct DLLEXPORT cerrPtrObject{
  template < typename T > void operator() ( const T * p ) const { std::cerr << *p << " " << std::endl; }
};
struct DLLEXPORT coutPtr{
  template < typename T > void operator() ( const T * p ) const { std::cout << p << " " << std::endl; }
};
struct DLLEXPORT coutPtrObject{
  template < typename T > void operator() ( const T * p ) const { std::cout << *p << " " << std::endl; }
};

template < typename Set > inline void intersectionSet( Set & dest, const Set & a, const Set & b ){
    dest.clear();
    set_intersection( a.begin(), a.end(), b.begin(), b.end(), std::inserter( dest, dest.begin() ) );
}
template < typename Set > inline void intersectionSet( Set & dest, const Set & a, const Set & b,                                                                   
                                                       const Set & c ){
    dest.clear();
    set_intersection( a.begin(), a.end(), b.begin(), b.end(), std::inserter( dest, dest.begin() ) );
    Set tmp( dest );
    dest.clear();
    set_intersection( tmp.begin(), tmp.end(), c.begin(), c.end(), std::inserter( dest, dest.begin() ) );
}
template < typename Set > inline void intersectionSet( Set & dest, const std::vector < Set > & a ){
    if ( a.size() > 1 ) {
        intersectionSet( dest, a[ 0 ], a[ 1 ] );
        for ( uint i = 2; i < a.size(); i ++ ){
            Set tmp( dest );
            dest.clear();
            set_intersection( tmp.begin(), tmp.end(), a[ i ].begin(), a[ i ].end(),
                              std::inserter( dest, dest.begin() ) );
        }
    } else {
        dest.clear();
    }
}

template < class T > class IncrementSequence {
public:
    IncrementSequence( T initialValue = 0 ) : value_( initialValue ) {   }
    inline T operator() () { return value_++; }
private:
    T value_;
};

} // namespace GIMLI

#endif // _GIMLI_GIMLI__H
