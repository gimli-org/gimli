/***************************************************************************
 *   Copyright (C) 2006-2016 by the resistivity.net development team       *
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

#include "stopwatch.h"

#include <iostream>

namespace GIMLI{

Stopwatch::Stopwatch( bool start ) {
    state_ = undefined;
    if ( start ) this->start();
}

Stopwatch::~Stopwatch() {
}

void Stopwatch::start(){
    ftime( & starttime );
    state_ = running;
    cCounter_.tic();
}

void Stopwatch::stop( bool verbose ){
    ftime( & stoptime );
    state_ = halted;
    if ( verbose ) std::cout << "time: " << duration() << "s" << std::endl;
}

void Stopwatch::restart(){
    stop();
    start();
}

void Stopwatch::reset(){
    restart();
}

double Stopwatch::duration( bool res ){
    if ( state_ == undefined ) std::cerr << "Stopwatch not started!" << std::endl;
    if ( state_ == running ) ftime( &stoptime );
    double t = ( stoptime.time - starttime.time ) + double( stoptime.millitm - starttime.millitm ) / 1000.0;
    if ( res ) restart();
    return t;
}

size_t Stopwatch::cycles( bool res ){
    if ( state_ == undefined ) std::cerr << "Stopwatch not started!" << std::endl;
    size_t t = 0;
    if ( state_ == running ) t = cCounter_.toc();
    if ( res ) restart();
    return t;
}

} // namespace GIMLI
