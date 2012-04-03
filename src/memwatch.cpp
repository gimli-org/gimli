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

#ifdef HAVE_PROC_READPROC
    #include <proc/readproc.h>
#endif

namespace GIMLI { 

// Global static pointer used to ensure a single instance of the class.
template <typename MemWatch> MemWatch* Singleton < MemWatch>::pInstance_ = NULL;

MemWatch::MemWatch( ){ 
    last_ = inUse(); 
    swatchAll_ = new Stopwatch( true );
    swatchDur_ = new Stopwatch( true );
} 
     
MemWatch::~MemWatch( ){ 
    delete swatchAll_; swatchAll_ = NULL;
    delete swatchDur_; swatchDur_ = NULL;
} 
     
double MemWatch::current( ){
    double ret = inUse() - last_;
    last_ = inUse();
    return ret;
}
    
double MemWatch::inUse( ) {
#ifdef HAVE_PROC_READPROC
    #ifdef HAVE_BOOST_THREAD_HPP
    boost::mutex::scoped_lock lock( mutex_ ); // slows down alot
    #endif
    
    struct proc_t usage;
    look_up_our_self( & usage );
    double ret = mByte( usage.vsize );
    return ret;
#else
    return 0;
#endif
}
    
void MemWatch::info( const std::string & str ) {
    #ifdef HAVE_PROC_READPROC
    if ( __GIMLI_DEBUG__ ){
    
        std::cout << "\t" << str << " Memory in use: abs: " << inUse() << " rel: " 
                    << current() << " MByte. t = " 
                    << swatchAll_->duration() << "/" << swatchDur_->duration( true ) << " s " <<  std::endl;
    }
    #endif
}

} // namespace GIMLI