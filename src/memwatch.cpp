/***************************************************************************
 *   Copyright (C) 2012 by the resistivity.net development team            *
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

#include "memwatch.h"
#include "stopwatch.h"

#include <iostream>

#ifdef WIN32_LEAN_AND_MEAN
    #include <psapi.h>
#elif READPROC_FOUND || defined(HAVE_PROC_READPROC)
    #include <proc/readproc.h>
    #define USE_PROC_READPROC TRUE
#else
    #define USE_PROC_READPROC 0
#endif

#if USE_BOOST_THREAD
    #include <boost/thread.hpp>
    /*! Lock proc reading to be thread safe */
    static boost::mutex __readproc__mutex__;
#endif

namespace GIMLI {

// Global static pointer used to ensure a single instance of the class.
template < > DLLEXPORT MemWatch * Singleton < MemWatch>::pInstance_ = NULL;

MemWatch::MemWatch(){
    last_ = inUse();
    swatchAll_ = new Stopwatch(true);
    swatchDur_ = new Stopwatch(true);
}

MemWatch::~MemWatch(){
    delete swatchAll_; swatchAll_ = NULL;
    delete swatchDur_; swatchDur_ = NULL;
}

double MemWatch::current(){
    double ret = inUse() - last_;
    last_ = inUse();
    return ret;
}

double MemWatch::inUse() {
 #if USE_BOOST_THREAD
        boost::mutex::scoped_lock lock(__readproc__mutex__);
 #endif

 #ifdef WIN32_LEAN_AND_MEAN

    PROCESS_MEMORY_COUNTERS pmc;

    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))){
        double ret = mByte(pmc.WorkingSetSize);
        return ret;
    } else { return 0; }

#else

#if USE_PROC_READPROC
    struct proc_t usage;
    look_up_our_self(& usage);
    // vsize .. number of pages of virtual memory ...
     //__MS("vsize: " << usage.vsize/1024) // virtual memory
//      __MS("size: " << usage.size/1024)//   total # of pages of memory
//      __MS("resident: " << usage.resident/1024)//         number of resident set (non-swapped) pages (4k)
//      __MS("share: " << usage.share/1024)
//      __MS("rss: " << usage.rss/1024)
    double ret = mByte(usage.vsize);
    return ret;
    #else // no windows and no libproc

    #endif // no libproc
#endif // no windows
    return 0;
}

void MemWatch::info(const std::string & str){
    if (debug()){
#if defined(WIN32_LEAN_AND_MEAN) || USE_PROC_READPROC
        std::cout << "\t" << str << " Memory "
#if USE_BOOST_THREAD
                    << "(mt)"
#endif // HAVE_BOOST_THREAD_HPP
                    << " in use: abs: " << inUse() << " rel: "
                    << current() << " MByte. t = "
                    << swatchAll_->duration() << "/" << swatchDur_->duration(true) << " s " <<  std::endl;
    #else
        std::cout << "\t" << str << " Memory no info"  <<  std::endl;
    #endif
    }
}

} // namespace GIMLI
