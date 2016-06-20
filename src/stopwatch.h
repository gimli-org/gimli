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

#ifndef _GIMLI_STOPWATCH__H
#define _GIMLI_STOPWATCH__H

#include "gimli.h"
#include <sys/timeb.h>

#if defined(__i386__)
static __inline__ size_t rdtsc__( void ){
    size_t x;
    __asm__ volatile (".byte 0x0f, 0x31" : "=A" ( x ));
    return x;
}
#elif defined(__x86_64__)
static __inline__ size_t rdtsc__( void ){
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ( (size_t )lo)|( ((size_t )hi)<<32 );
}
#else
static inline size_t rdtsc__( void ){
    return 0;
}
#endif

namespace GIMLI{

class DLLEXPORT CycleCounter{
public:
    CycleCounter() : var_( 0 ) {}

    ~CycleCounter(){}

    inline void tic( ){ var_ = rdtsc__(); }

    inline size_t toc( ) const { return ( rdtsc__() - var_ ); }

protected:

    size_t var_;
};

class DLLEXPORT Stopwatch {
public:
    Stopwatch( bool start = false );

    ~Stopwatch();

    void start();

    void stop( bool verbose=false );

    /*! Restart the stopwatch.*/
    void restart();

    /*! Reset the stopwatch, same like \ref restart.*/
    void reset();

    /*! Returns the current duration in seconds. Optional you can restart the stopwatch.*/
    double duration( bool restart=false );

    /*! Returns the cpu cycles. Optional you can restart the stopwatch.*/
    size_t cycles( bool restart=false );

    const CycleCounter & cycleCounter() const { return cCounter_; }

protected:
    timeb starttime, stoptime;
    enum watchstate {undefined,halted,running} state_;
    CycleCounter cCounter_;
};

#define TIC__ std::cout.precision( 12 ); GIMLI::Stopwatch __swatch__( true );
#define TOC__ std::cout << __swatch__.duration( true ) << std::endl;
#define toc__ __swatch__.duration()

} // namespace GIMLI

#endif // _GIMLI_STOPWATCH__H


