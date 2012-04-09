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

#include "platform.h"

#ifdef WIN32_LEAN_AND_MEAN

#else
    #include <unistd.h>
#endif

#include <iostream>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace GIMLI{

int numberOfCPU(){

    long nprocs = -1;
    long nprocs_max = -1;
#if defined(WINDOWS) || defined(_WIN32)
    #ifndef _SC_NPROCESSORS_ONLN
        SYSTEM_INFO info;
        GetSystemInfo(&info);
    #define sysconf(a) info.dwNumberOfProcessors
    #define _SC_NPROCESSORS_ONLN
    #endif
#endif // no windows

#ifdef _SC_NPROCESSORS_ONLN
    nprocs = sysconf(_SC_NPROCESSORS_ONLN);

    if ( nprocs < 1 ) {
        std::cerr << "Could not determine number of CPUs online:" << std::strerror(errno) << std::endl;
    } else {
        nprocs = 1;
    }
    nprocs_max = sysconf(_SC_NPROCESSORS_CONF);

    if ( nprocs_max < 1 ) {
        std::cerr << "Could not determine number of CPUs configured:" << std::strerror(errno) << std::endl;
    }
#else
    std::cerr << "Could not determine number of CPUs (!_SC_NPROCESSORS_ONLN)" << std::endl;
    nprocs = 1;
#endif
    return nprocs_max;
}

} // namespace GIMLI{
