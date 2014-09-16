/***************************************************************************
 *   Copyright (C) 2006-2014 by the resistivity.net development team       *
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

#ifndef _GIMLI_SOLVERWRAPPER__H
#define _GIMLI_SOLVERWRAPPER__H

#include "gimli.h"

namespace GIMLI{

class DLLEXPORT SolverWrapper{
public:
    SolverWrapper();

    SolverWrapper(const RSparseMatrix & S, bool verbose=false);

    SolverWrapper(const CSparseMatrix & S, bool verbose=false);
    
    virtual ~SolverWrapper();

    virtual int solve(const RVector & rhs, RVector & solution) = 0;
    
    virtual int solve(const CVector & rhs, CVector & solution){ THROW_TO_IMPL return 0;}

protected:
   
    bool dummy_;
    bool verbose_;
    bool isComplex_;
    uint dim_;
    long nVals_;
    double dropTol_;
    double tolerance_;
    double maxiter_;

};

} //namespace GIMLI;

#endif // _GIMLI_SOLVERWRAPPER__H
