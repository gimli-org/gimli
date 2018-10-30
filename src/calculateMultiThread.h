/******************************************************************************
 *   Copyright (C) 2005-2018 by the GIMLi development team                    *
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

#ifndef _GIMLI_CALCULATE_MULTI_THREAD__H
#define _GIMLI_CALCULATE_MULTI_THREAD__H

#include "gimli.h"

#ifdef USE_BOOST_THREAD
    #include <boost/thread.hpp>
#else
    #include <thread>
#endif

namespace GIMLI{

class BaseCalcMT{
public:
    BaseCalcMT(Index count=0, bool verbose=false)
        : verbose_(verbose), start_(0), end_(0), threadNumber_(count){
    }

    virtual ~BaseCalcMT(){ }

    void operator () () { calc(threadNumber_); }

    void setRange(Index start, Index end, Index threadNumber=0){
        start_ = start;
        end_ = end;
        if (threadNumber_ > 0){
            threadNumber_ = threadNumber;
        }
    }

    virtual void calc(Index tNr=0)=0;

protected:
    bool verbose_;
    Index start_;
    Index end_;
    Index threadNumber_;
};

template < class T > void distributeCalc(T calc, uint nCalcs, uint nThreads, bool verbose=false){
    if (nThreads == 1){
        calc.setRange(0, nCalcs);
        calc();
    } else {
        uint singleCalcCount = (uint)ceil((double)nCalcs / (double)nThreads);

        std::vector < T > calcObjs;
        for (uint i = 0; i < nThreads; i ++){
            calcObjs.push_back(calc);
            Index start = singleCalcCount * i;
            Index end   = singleCalcCount * (i + 1);
            if (i == nThreads -1) end = nCalcs;
            if (debug()) std::cout << "Threaded calculation: " << i << ": " << start <<" " << end << std::endl;
            calcObjs.back().setRange(start, end, i);
            if (end >= nCalcs) break;
        }
#if USE_BOOST_THREAD
        boost::thread_group threads;
        for (uint i = 0; i < calcObjs.size(); i++) {
            if (debug()) std::cout << "start boost::thread: " << i << std::endl;
            threads.create_thread(calcObjs[i]);
        }
        threads.join_all();
#else

        std::vector<std::thread> threads;

        for (uint i = 0; i < calcObjs.size(); i++) {
            if (debug()) std::cout << "start std::thread: " << i << std::endl;
            threads.emplace_back(calcObjs[i]);
        }

        for (auto & th : threads) if (th.joinable()) th.join();

//         std::vector boost::thread_group threads;
//         for (uint i = 0; i < nThreads; i++) calcObjs[i]();
#endif
    }
}

} // namespace GIMLI{

#endif //_GIMLI_IPC_CLIENT__H
