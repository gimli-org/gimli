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

#ifndef _GIMLI_CALCULATE_MULTI_THREAD__H
#define _GIMLI_CALCULATE_MULTI_THREAD__H

#include "gimli.h"

#ifdef USE_BOOST_THREAD
#include <boost/thread.hpp>
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
        }
#if USE_BOOST_THREAD
        boost::thread_group threads;
        for (uint i = 0; i < nThreads; i++) {
            if (debug()) std::cout << "start thread: " << i << std::endl;
            threads.create_thread(calcObjs[i]);
        }
        threads.join_all();
#else
        for (uint i = 0; i < nThreads; i++) calcObjs[i]();
#endif
    }
}

} // namespace GIMLI{

#endif //_GIMLI_IPC_CLIENT__H
