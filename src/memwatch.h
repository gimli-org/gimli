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

#ifndef _GIMLI_MEMWATCH__H
#define _GIMLI_MEMWATCH__H

#include "gimli.h"

#ifdef HAVE_BOOST_THREAD_HPP
    #include <boost/thread.hpp>
#endif

namespace GIMLI{

/*! Convert byte into KByte. */
inline double kByte( long byte ){ return double( byte / 1024.0 ); }

/*! Convert byte into MByte */
inline double mByte( long byte ){ return double( byte / ( 1024.0 * 1024.0 ) ); }

//! Memory watch.
/*! Class that might help debugging memory usage.
 * Informations are taken from /proc system, so only available for linux systems.
 * This is a singleton class to ensure a single instance.
 * To call it use e.g.: MemWatch::instance().info( WHERE );*/
class DLLEXPORT MemWatch  : public Singleton< MemWatch > {
public:
    friend class Singleton< MemWatch >;

    /*! Return the current memory usage of the process. Values are in MByte. */
    double inUse( );
        
    /*! Return the current memory usage relative to the last call of this method. Values are in MByte. */
    double current( );
            
    /*! Shows the current and the relative memory usage. */
    void info( const std::string & str = "" );
         
protected:
    double last_;
    
private:
    /*! Private so that it can not be called */
    MemWatch( ); 
    /*! Private so that it can not be called */
    virtual ~MemWatch( ); 
    /*! Copy constructor is private, so don't use it */
    MemWatch( const MemWatch & ){}; 
    /*! Assignment operator is private, so don't use it */
    void operator = ( const MemWatch & ){ };

    Stopwatch * swatchAll_;
    Stopwatch * swatchDur_;
    
    /*! Lock proc reading to be thread safe */
    #ifdef HAVE_BOOST_THREAD_HPP
    boost::mutex mutex_;
    #endif
};

   
} // namespace GIMLI

#endif