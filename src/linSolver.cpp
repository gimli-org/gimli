/***************************************************************************
 *   Copyright (C) 2007-2011 by the resistivity.net development team       *
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

LinSolver::LinSolver(  ) 
    : solver_( NULL ), verbose_( false ), rows_( 0 ), cols_( 0 ) {
    setSolverType( AUTOMATIC );
}

LinSolver::LinSolver( DSparseMatrix & S, bool verbose ) 
    : solver_( NULL ), verbose_( verbose ) {
    setSolverType( AUTOMATIC );
    setMatrix(S);
}
  
LinSolver::LinSolver( DSparseMatrix & S, SolverType solverType, bool verbose ) 
    : solver_( NULL ), verbose_( verbose ) {
    setSolverType( solverType );
    setMatrix(S);
}
  
LinSolver::~LinSolver(){
    delete solver_;
}

void LinSolver::setSolverType(SolverType solverType){
   solverType_ = solverType;
   if ( solverType_ == AUTOMATIC ){
        solverType_ = UNKNOWN;

        if ( LDLWrapper::valid() ){
            solverType_ = LDL;
        } 
        if ( CHOLMODWrapper::valid() ){
            solverType_ = CHOLMOD;
        }
    }
}

void LinSolver::setMatrix(DSparseMatrix & S, bool verbose){
    verbose_ = verbose;
    initialize_(S);
}

void LinSolver::initialize_(DSparseMatrix & S){
    rows_ = S.rows();
    cols_ = S.cols();
    setSolverType(solverType_);
    
    switch( solverType_ ){
    case LDL:     solver_ = new LDLWrapper( S, verbose_ ); break;
    case CHOLMOD: solver_ = new CHOLMODWrapper(S, verbose_); break;
    case UNKNOWN: 
    default:
        std::cerr << WHERE_AM_I << " no valid solver found"  << std::endl;
    }
}

void LinSolver::solve( const RVector & rhs, RVector & solution ){
    solution.resize( rows_ );
    if ( rhs.size() != cols_ ){
        std::cerr << WHERE_AM_I << " rhs size mismatch: " << cols_ << "  " << rhs.size() << std::endl;
    }
    if ( solver_ ) solver_->solve( rhs, solution );
}

std::string LinSolver::solverName() const {
  switch( solverType_ ){
  case LDL:     return "LDL"; break;
  case CHOLMOD: return "CHOLMOD"; break;
  case UNKNOWN: 
  default: return " no valid solver installed";
  }
}

} // namespace GIMLI

