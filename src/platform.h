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

#ifndef _GIMLI_PLATFORM__H
#define _GIMLI_PLATFORM__H

#if defined(WINDOWS) || defined(_WIN32)
#define PATHSEPARATOR "\\"

#ifdef DLL_EXPORT
#define DLLEXPORT __declspec(dllexport)
#else /* NO BUILDING_DLL */
#define DLLEXPORT __declspec(dllimport)
#endif /* NO BUILDING_DLL */

#include <windows.h>
#undef near
#undef far

inline BOOL APIENTRY DllMain (HINSTANCE hInst     /* Library instance handle. */ ,
                       DWORD reason        /* Reason this function is being called. */ ,
                       LPVOID reserved     /* Not used. */ ){
    switch ( reason )    {
      case DLL_PROCESS_ATTACH:
        break;
      case DLL_PROCESS_DETACH:
        break;
      case DLL_THREAD_ATTACH:
        break;
      case DLL_THREAD_DETACH:
        break;
    }
    /* Returns TRUE on success, FALSE on failure */
    return TRUE;
}

//inline bool isnan( double x ) { return x != x; }

#else /* ELSE NO WINDOWS */
#define PATHSEPARATOR "/"
#define DLLEXPORT
#endif /* NO WINDOWS */

/*!
 * Return the numbers of virtual CPUS on this system
 */

namespace GIMLI{
    int numberOfCPU();
}
//DLLEXPORT int numberOfCPU();

#if defined ( __APPLE__ ) || ( defined (__SVR4) && defined (__sun) )
// #include <ieeefp.h>
// #include <unistd.h>
// inline int isinf( double x ){ return !finite(x) && x==x; }
                              // /*
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
