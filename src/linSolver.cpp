/***************************************************************************
 *   Copyright (C) 2007-2011 by the GIMLi development team       *
 *   Carsten Rücker carsten@resistivity.net                                *
 *   Thomas Günther thomas@resistivity.net                                 *
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

#include "linSolver.h"
#include "sparsematrix.h"
#include "ldlWrapper.h"
#include "cholmodWrapper.h"

namespace GIMLI{

LinSolver::LinSolver(bool verbose) 
    : verbose_(verbose){
    init_();
    setSolverType(AUTOMATIC);
}

LinSolver::LinSolver(RSparseMatrix & S, bool verbose) 
    : verbose_(verbose) {
    init_();
    setSolverType(AUTOMATIC);
    setMatrix(S);
}

LinSolver::LinSolver(RSparseMapMatrix & S, bool verbose)
    : verbose_(verbose) {
    init_();
    setSolverType(AUTOMATIC);
    cacheMatrix_ = new RSparseMatrix(S);
    
    setMatrix(dynamic_cast < RSparseMatrix &>( *cacheMatrix_));
}
  
LinSolver::LinSolver(RSparseMatrix & S, SolverType solverType, bool verbose) 
    : verbose_(verbose) {
    init_();
    setSolverType(solverType);
    setMatrix(S);
}
  
LinSolver::LinSolver(CSparseMatrix & S, bool verbose) 
    : verbose_(verbose) {
    init_();
    setSolverType(AUTOMATIC);
    setMatrix(S);
}
  
LinSolver::LinSolver(CSparseMatrix & S, SolverType solverType, bool verbose) 
    : verbose_(verbose) {
    init_();
    setSolverType(solverType);
    setMatrix(S);
}

void LinSolver::init_(){
    rows_ = 0;
    cols_ = 0;
    solver_ = 0;
    cacheMatrix_ = 0;
}

LinSolver::~LinSolver(){
    if (cacheMatrix_) {
        delete cacheMatrix_;
        cacheMatrix_ = 0;
    }
    if (solver_) {
        delete solver_;
        solver_ = 0;
    }
}

void LinSolver::setSolverType(SolverType solverType){
   solverType_ = solverType;
   if (solverType_ == AUTOMATIC){
        solverType_ = UNKNOWN;

        if (LDLWrapper::valid()){
            solverType_ = LDL;
        } 
        if (CHOLMODWrapper::valid()){
            solverType_ = CHOLMOD;
        }
    }
}

void LinSolver::setMatrix(RSparseMatrix & S, int stype){
    initialize_(S, stype);
}

void LinSolver::setMatrix(CSparseMatrix & S, int stype){
    initialize_(S, stype);
}

// template <> void LinSolver::solve(const RVector & rhs, RVector & solution);
// template <> void LinSolver::solve(const CVector & rhs, CVector & solution);
// template <> RVector LinSolver::solve(const RVector & rhs);
// template <> CVector LinSolver::solve(const CVector & rhs);

void LinSolver::solve(const RVector & rhs, RVector & solution){
    solution.resize(rows_);
    if (rhs.size() != cols_){
        std::cerr << WHERE_AM_I << " rhs size mismatch: " << cols_ << "  " << rhs.size() << std::endl;
    }
    if (solver_) solver_->solve(rhs, solution);
}

RVector LinSolver::solve(const RVector & rhs){
    RVector solution(rhs.size());
    if (solver_) solver_->solve(rhs, solution);
    return solution;
}

void LinSolver::solve(const CVector & rhs, CVector & solution){
    solution.resize(rows_);
    if (rhs.size() != cols_){
        std::cerr << WHERE_AM_I << " rhs size mismatch: " << cols_ << "  " << rhs.size() << std::endl;
    }
    if (solver_) solver_->solve(rhs, solution);
}

CVector LinSolver::solve(const CVector & rhs){
    CVector solution(rhs.size());
    if (solver_) solver_->solve(rhs, solution);
    return solution;
}

void LinSolver::initialize_(RSparseMatrix & S, int stype){
    rows_ = S.rows();
    cols_ = S.cols();
    setSolverType(solverType_);

    switch(solverType_){
        case LDL:     solver_ = new LDLWrapper(S, verbose_); break;
        case CHOLMOD: solver_ = new CHOLMODWrapper(S, verbose_, stype); break;
        case UNKNOWN: 
    default:
            std::cerr << WHERE_AM_I << " no valid solver found"  << std::endl;
    }
}

void LinSolver::initialize_(CSparseMatrix & S, int stype){
    rows_ = S.rows();
    cols_ = S.cols();
    setSolverType(solverType_);

    switch(solverType_){
        case LDL:     solver_ = new LDLWrapper(S, verbose_); break;
        case CHOLMOD: solver_ = new CHOLMODWrapper(S, verbose_, stype); break;
        case UNKNOWN: 
    default:
            std::cerr << WHERE_AM_I << " no valid solver found"  << std::endl;
    }
}

    
std::string LinSolver::solverName() const {
  switch(solverType_){
  case LDL:     return "LDL"; break;
  case CHOLMOD: return "CHOLMOD"; break;
  case UNKNOWN: 
  default: return " no valid solver installed";
  }
}

} // namespace GIMLI

