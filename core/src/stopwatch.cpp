/******************************************************************************
 *   Copyright (C) 2006-2024 by the GIMLi development team                    *
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
#include "vector.h"

#include <chrono>
#include <thread>

#include <iostream>

#include <ranges>

namespace GIMLI{

Stopwatch::Stopwatch(bool start) : _store(nullptr) {
    //newPtr(_store);
    _state = undefined;
    if (start) this->start();
    this->_store = new RVector();
}

Stopwatch::~Stopwatch() {
    // deletePtr(this->_store);
    delete this->_store;
}

void Stopwatch::start(){
    this->_start = std::chrono::high_resolution_clock::now();
    _state = running;
    _cCounter.tic();
}

void Stopwatch::stop(bool verbose){
    this->_stop = std::chrono::high_resolution_clock::now();
    _state = halted;
    if (verbose) std::cout << "time: " << duration() << "s" << std::endl;
}

void Stopwatch::restart(){
    stop();
    start();
}

void Stopwatch::reset(){
    restart();
    this->_store->clear();
}

double Stopwatch::duration(bool restart){
    std::chrono::time_point<std::chrono::high_resolution_clock> now;
    
    if (_state == undefined) {
        log(Error, "Stopwatch not started!");
    }
        
    if (_state == running) {
        //ftime(&stoptime);
        // double t = (stoptime.time - starttime.time) + 
        //     double(stoptime.millitm - starttime.millitm) / 1000.0;
            
        now = std::chrono::high_resolution_clock::now();
    } else {
        now = this->_stop;
    }
    
    std::chrono::duration< double > t = now - this->_start;

    if (restart) this->restart();
    return t.count();
}

void Stopwatch::store(bool restart){
    this->_store->push_back(this->duration(restart));
}

size_t Stopwatch::cycles(bool res){
    if (_state == undefined) {
        log(Error, "Stopwatch not started!");
    }
    size_t t = 0;
    if (_state == running) t = _cCounter.toc();
    if (res) restart();
    return t;
}

void waitms(Index ms, Index count){
    for (Index i = 0; i < count; i ++ ){
        std::this_thread::sleep_for(std::chrono::milliseconds(ms));
    }
    return;
}
void waitmsOMP(Index ms, Index count){

    #pragma omp parallel for
    for (Index i = 0; i < count; i ++ ){
        std::this_thread::sleep_for(std::chrono::milliseconds(ms));
    }
    return;
}

#include <unistd.h>
void waitus(Index ms, Index count){
    for (Index i = 0; i < count; i ++ ){
        usleep(ms);
        // std::this_thread::sleep_for(std::chrono::microseconds(ms));
    }
    return;
}
void waitusOMP(Index ms, Index count){

    #pragma omp parallel for
    for (Index i = 0; i < count; i ++ ){
        std::this_thread::sleep_for(std::chrono::microseconds(ms));
    }
    return;
}


template < > DLLEXPORT Swatches * Singleton < Swatches >::pInstance_ = NULL;

Swatches::Swatches(){
// __MS("Swatches()")
}

Swatches::~Swatches(){
// __MS("~Swatches()")
}

Stopwatch & Swatches::operator[](const std::string & key) {
    // __MS("[]", this, key)
    if (this->_sw.count(key) == 0){
        this->_sw[key] = new Stopwatch(true);
    }
    _trace = key;
    return *this->_sw[key];
}

std::vector < std::string > Swatches::keys(){
    // __MS(this)
    std::vector< std::string > keys; 
    keys.reserve(this->_sw.size()); 
    for (auto &kv: this->_sw){ 
        keys.push_back(kv.first); 
    } 

    // auto kv = std::views::keys(this->_sw); ## since c++20
    // std::vector<std::string> keys{ kv.begin(), kv.end() };
    // __M
    return keys;
}

std::vector < const Stopwatch * > Swatches::vals(){
    std::vector< const Stopwatch * > vals; 
    vals.reserve(this->_sw.size()); 
    for (auto &kv: this->_sw){ 
        vals.push_back(kv.second); 
    } 

    // auto kv = std::views::keys(this->_sw); ## since c++20
    // std::vector<std::string> keys{ kv.begin(), kv.end() };
    return vals;
}

std::vector < std::pair< std::string, Stopwatch * > > Swatches::items(){
    std::vector< std::pair< std::string, Stopwatch * > > items; 
    items.reserve(this->_sw.size()); 
    for (auto &kv: this->_sw){ 
        items.push_back(std::pair< std::string, Stopwatch * >(kv.first, kv.second)); 
    } 
    return items;
}
    
void Swatches::remove(const std::string & key, bool isRoot){
    if (isRoot == false){
        Stopwatch * s = this->_sw[key];
        this->_sw.erase(key);
        delete s;
    } else {
        //     for k in list(self._sw.keys()):
        //         if k.startswith(key):
        //             self._sw.pop(k, None)
        THROW_TO_IMPL
    }
}

TicToc::TicToc(const std::string & name, bool reset){
    
    this->_parentTrace = Swatches::instance().trace();
    std::string curTrace;

    if (this->_parentTrace.size() > 0){
        curTrace = this->_parentTrace + '/' + name;
    } else{
        curTrace = name;
    }

    if (reset == true){
        THROW_TO_IMPL
        Swatches::instance().remove(this->_parentTrace, true);
    }
    // __MS(this->_parentTrace, "start")
    
    this->_sw = & Swatches::instance()[curTrace];
    Swatches::instance().setTrace(curTrace);
    this->_sw->start();
}

TicToc::~TicToc(){
    // __MS(this->_parentTrace, "store", this->_sw->duration())
    this->_sw->store();
    Swatches::instance().setTrace(this->_parentTrace);
}


} // namespace GIMLI