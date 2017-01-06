/***************************************************************************
 *   Copyright (C) 2006-2014 by the GIMLi development team       *
 *   Carsten R�cker carsten@resistivity.net                                *
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

