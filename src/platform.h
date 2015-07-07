/***************************************************************************
 *   Copyright (C) 2006-2013 by the resistivity.net development team       *
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

#ifndef _GIMLI_PLATFORM__H
#define _GIMLI_PLATFORM__H

#if defined(WINDOWS) || defined(_WIN32) || defined(WIN32)
    //#error checkme
    
    #define PATHSEPARATOR "\\"

	#define WIN32_LEAN_AND_MEAN
    #define BOOST_USE_WINDOWS_H
    //#define BOOST_PLAT_WINDOWS_RUNTIME
    
#if defined(DLL_EXPORT) || defined(gimli_EXPORTS)
    //* We are building this library
    //#warning (We are building this library)
    #define DLLEXPORT __declspec(dllexport)
#else
    //* We are using this library
    //#warning (We are using this library)
    #define DLLEXPORT __declspec(dllimport)
#endif

    // Don't let win32api windef.h define min and max as macros
    // if included after c++config.h.
    #ifndef NOMINMAX
    #define NOMINMAX
#endif
    #include <stddef.h>
#ifndef PYGIMLI_GCCXML
    #include <windows.h>
#endif
    #undef near
    #undef far
#else /* ELSE NO WINDOWS */
    #define PATHSEPARATOR "/"
    #define DLLEXPORT
#endif /* NO WINDOWS */

/*!
 * Return the numbers of virtual CPUS on this system
 */

#include <cmath>
namespace GIMLI{
DLLEXPORT int numberOfCPU();

// Microsoft Visual C++ 10 does not provide some C99 functions
#if defined(_MSC_VER) 
template< typename T > T rint( T x ){ return std::floor(x + 0.5); }
template< typename T > inline bool isinf(T value){return value >= std::numeric_limits<T>::min() && 
														 value <= std::numeric_limits<T>::max();}
template< typename T > inline bool isnan(T value){return (value) != (value);}
#elif defined(_WIN32) // mingw
template< typename T > inline bool isinf(T value){return std::isinf(value);}
template< typename T > inline bool isnan(T value){return std::isnan(value);}
#else
template< typename T > inline bool isinf(T value){return std::isinf(value);}
template< typename T > inline bool isnan(T value){return std::isnan(value);}
#endif

}

#if defined ( __APPLE__ ) || ( defined (__SVR4) && defined (__sun) )
// #include <ieeefp.h>
// inline int isinf( double x ){ return !finite(x) && x==x; }
                              // /*
#include <unistd.h>
#include <sys/types.h>
                            
//   Apple (OS X) and Sun systems declare getopt in unistd.h,
//   other systems (Linux) use getopt.h
// */
// #if defined ( __APPLE__ ) || ( defined (__SVR4) && defined (__sun) )
// #include <unistd.h>
// struct option {
// # if defined __STDC__ && __STDC__
//   const char *name;
// # else
//   char *name;
// # endif
//   /* has_arg can't be an enum because some compilers complain about
//      type mismatches in all the code that assumes it is an int.  */
//   int has_arg;
//   int *flag;
//   int val;
// };
// #else

#endif

#endif // _GIMLI_PLATFORM__H
