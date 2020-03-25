/******************************************************************************
 *   Copyright (C) 2006-2020 by the GIMLi development team                    *
 *   Carsten Rücker carsten@resistivity.net                                   *
 *                                                                            *
 *   Licensed under the Apache License, Version 2.0 (the "License");          *
 *   you may not use this file except in compliance with the License.         *
 *   You may obtain a copy of the License at                                  *
 *                                                                            *
 *       http://www.apache.org/licenses/LICENSE-2.0                           *
 *                                                                            *
 *   Unless required by applicable law or agreed to in writing, software      *
 *   distributed under the License is distributed on an "AS IS" BASIS,        *
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 *   See the License for the specific language governing permissions and      *
 *   limitations under the License.                                           *
 *                                                                            *
 ******************************************************************************/

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
