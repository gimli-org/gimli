/***************************************************************************
 *   Copyright (C) 2007-2011 by the resistivity.net development team       *
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

#ifdef CHOLMOD_FOUND
    #if CHOLMOD_FOUND==TRUE
        #define USE_CHOLMOD 1
    #endif
#endif

#ifdef USE_CHOLMOD
    #define CHOLMOD_MAXMETHODS 9
    #define UF_long long

    #include <cholmod.h>
#else
struct cholmod_sparse;
struct cholmod_factor;
struct cholmod_common;
#endif

namespace GIMLI{

class DLLEXPORT CHOLMODWrapper : public SolverWrapper {
public:
  CHOLMODWrapper(DSparseMatrix & S, bool verbose = false);
  ~CHOLMODWrapper();

  static bool valid();

  int factorise();
  int solve( const RVector & rhs, RVector & solution );

protected:
  int initialize_( DSparseMatrix & S );

  cholmod_common *c_;
  cholmod_sparse *A_;
  cholmod_factor *L_;
};

} //namespace GIMLI;

#endif // _GIMLI_CHOLMODWRAPPER__H
