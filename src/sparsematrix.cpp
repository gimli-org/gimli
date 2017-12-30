/******************************************************************************
 *   Copyright (C) 2007-2018 by the GIMLi development team                    *
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

#include "gimli.h"
#include "sparsematrix.h"

namespace GIMLI{

template<>
void SparseMatrix< double >::copy_(const SparseMapMatrix< double, Index > & S){
    this->clear();
    Index col = 0, row = 0;
    double val;
    cols_ = S.cols();
    rows_ = S.rows();
    std::vector < std::map < Index, double > > idxMap(S.cols());

    for (typename SparseMapMatrix< double, Index>::const_iterator
        it = S.begin(); it != S.end(); it ++){
        row = S.idx1(it);
        col = S.idx2(it);
        val = S.val(it);
        idxMap[col].insert(std::pair< Index, double >(row, val));
    }
    colPtr_.resize(S.cols() + 1);
    rowIdx_.resize(S.nVals());
    vals_.resize(S.nVals());
    stype_  = S.stype();

    colPtr_[0] = 0;

    Index colCounter = 0, rowCounter = 0;
    for (typename std::vector < std::map < Index, double > >::iterator
         it = idxMap.begin(); it != idxMap.end(); it++){

        for (typename std::map < Index, double >::iterator
             itR = (*it).begin(); itR != (*it).end(); itR++){
            rowIdx_[rowCounter] = itR->first;
            vals_[rowCounter] = (double)itR->second;
            rowCounter ++;
        }
        colCounter ++;
        colPtr_[colCounter] = rowCounter;
    }
    valid_ = true;
}


} // namespace GIMLI{


