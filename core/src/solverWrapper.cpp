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

#include "solverWrapper.h"
#include "sparsematrix.h"

namespace GIMLI{

SolverWrapper::SolverWrapper( ){ dummy_ = true; }

SolverWrapper::SolverWrapper(const RSparseMatrix & S,
                             bool verbose)
    : verbose_(verbose){
    dummy_ = true;
    dim_ = S.size();
    nVals_ = S.nVals();
    tolerance_ = 1e-12;
    dropTol_ = 0.0;
    isComplex_ = false;
}

SolverWrapper::SolverWrapper(const CSparseMatrix & S,
                             bool verbose)
    : verbose_(verbose){
    dummy_ = true;
    dim_ = S.size();
    nVals_ = S.nVals();
    tolerance_ = 1e-12;
    dropTol_ = 0.0;
    isComplex_ = true;
}

SolverWrapper::~SolverWrapper(){ }


} //namespace GIMLI;

