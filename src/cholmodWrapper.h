/***************************************************************************
 *   Copyright (C) 2007-2017 by the GIMLi development team       *
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

#ifndef _GIMLI_CHOLMODWRAPPER__H
#define _GIMLI_CHOLMODWRAPPER__H

#include "gimli.h"
#include "solverWrapper.h"

namespace GIMLI{

class DLLEXPORT CHOLMODWrapper : public SolverWrapper {
public:
    CHOLMODWrapper(RSparseMatrix & S, bool verbose=false, int stype=-2);
    
    CHOLMODWrapper(CSparseMatrix & S, bool verbose=false, int stype=-2);
    
    virtual ~CHOLMODWrapper();

    static bool valid();

    int factorise();
    
    virtual int solve(const RVector & rhs, RVector & solution);
    
    virtual int solve(const CVector & rhs, CVector & solution);

protected:
    void init();
    
    int initializeMatrix_(RSparseMatrix & S);
    
    int initializeMatrix_(CSparseMatrix & S);

    template < class ValueType > 
    void init_(SparseMatrix < ValueType > & S, int stype);

    template < class ValueType > 
    int initMatrixChol_(SparseMatrix < ValueType > & S, int xType);
    
    template < class ValueType > 
    int solveCHOL_(const Vector < ValueType > & rhs, Vector < ValueType > & solution);
    
    
    template < class ValueType > 
    int solveUmf_(const Vector < ValueType > & rhs, Vector < ValueType > & solution);
    
    
    int stype_;
    
    void *c_;
    void *A_;
    void *L_;
    
    bool useUmfpack_;
    void *Numeric_;
    void *NumericD_;
    int * Ap_;
    int * Ai_;
    int * ApR_;
    int * AiR_;
        
    RVector *AxV_;
    RVector *AzV_;
        
};

} //namespace GIMLI;

#endif // _GIMLI_CHOLMODWRAPPER__H
