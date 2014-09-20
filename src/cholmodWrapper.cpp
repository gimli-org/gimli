/***************************************************************************
 *   Copyright (C) 2007-2014 by the resistivity.net development team       *
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

#include "cholmodWrapper.h"
#include "vector.h"
#include "sparsematrix.h"

#if CHOLMOD_FOUND
    #define CHOLMOD_MAXMETHODS 9
    #define UF_long long

    #include <cholmod.h>
    #define USE_CHOLMOD 1
#endif

namespace GIMLI{

#if USE_CHOLMOD
    bool CHOLMODWrapper::valid() { return true; }
#else
    bool CHOLMODWrapper::valid() { return false; }
#endif

CHOLMODWrapper::CHOLMODWrapper(RSparseMatrix & S, bool verbose, int stype)
    : SolverWrapper(S, verbose){
    
    if (stype == -2){
        stype_ = S.stype();
    } else {
        stype_ = stype;
    }
    
        
  c_ = NULL;
  A_ = NULL;
  L_ = NULL;
#if USE_CHOLMOD
  c_ = new cholmod_common;
//   cholmod_common *c_;
//   cholmod_sparse *A_;
//   cholmod_factor *L_;
  
  int ret =  cholmod_start((cholmod_common*)c_);
  if (ret) dummy_ = false;

  initialize_(S);

#else
  std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
}

CHOLMODWrapper::CHOLMODWrapper(CSparseMatrix & S, bool verbose, int stype)
    : SolverWrapper(S, verbose){
    c_ = NULL;
    A_ = NULL;
    L_ = NULL;
  
    if (stype == -2){
        stype_ = S.stype();
    } else {
        stype_ = stype;
    }
    
#if USE_CHOLMOD
  c_ = new cholmod_common;
//   cholmod_common *c_;
//   cholmod_sparse *A_;
//   cholmod_factor *L_;
  
  int ret =  cholmod_start((cholmod_common*)c_);
  if (ret) dummy_ = false;

  initialize_(S);

#else
  std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
}
CHOLMODWrapper::~CHOLMODWrapper(){
#if USE_CHOLMOD
  cholmod_free_factor((cholmod_factor**)(&L_), (cholmod_common*)c_);
  //** We did not allocate the matrix so we dont need to free it
  //  cholmod_free_sparse(&A_, c_);
  cholmod_finish((cholmod_common*)c_);

  if (A_) delete (cholmod_sparse*)A_;
  if (c_) delete (cholmod_common*)c_;

#else
  std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
}

int CHOLMODWrapper::initialize_(CSparseMatrix & S){
  if (!dummy_){
#if USE_CHOLMOD

  //** We do not allocate the matrix since we use the allocated space from DSparsemarix
//    A_ = cholmod_allocate_sparse(dim_, dim_, nVals_, true, true, 1, CHOLMOD_REAL, c_) ;

    A_ = new cholmod_sparse;
    ((cholmod_sparse*)A_)->nrow  = S.nRows();         /* number of rows */
    ((cholmod_sparse*)A_)->ncol  = S.nCols();           /* number of columns */
    ((cholmod_sparse*)A_)->nzmax = S.nVals();               /* maximum number of entries */
    ((cholmod_sparse*)A_)->p     = (void*)S.colPtr();   /* column pointers (size n+1) or col indices (size nzmax) */
    ((cholmod_sparse*)A_)->i     = (void*)S.rowIdx();   /* row indices, size nzmax */
    
    ((cholmod_sparse*)A_)->x     = S.vals();     /* numerical values, size nzmax */

    ((cholmod_sparse*)A_)->stype = stype_;

    ((cholmod_sparse*)A_)->itype = CHOLMOD_INT;
    ((cholmod_sparse*)A_)->xtype = CHOLMOD_COMPLEX;  // data type for the pattern (Real, complex, zcomplex)
    ((cholmod_sparse*)A_)->dtype = CHOLMOD_DOUBLE; // data type for complex or real (float/double)
    ((cholmod_sparse*)A_)->packed = true;
    ((cholmod_sparse*)A_)->sorted = true; // testen, scheint schneller, aber hab ich das immer?
    
    factorise();
    
    return 1;
#else
    std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
  }
  return 0;
}

int CHOLMODWrapper::initialize_(RSparseMatrix & S){
  if (!dummy_){
#if USE_CHOLMOD

  //** We do not allocate the matrix since we use the allocated space from DSparsemarix
//    A_ = cholmod_allocate_sparse(dim_, dim_, nVals_, true, true, 1, CHOLMOD_REAL, c_) ;

    A_ = new cholmod_sparse;
    ((cholmod_sparse*)A_)->nrow  = dim_;         /* number of rows */
    ((cholmod_sparse*)A_)->ncol  = dim_;           /* number of columns */
    ((cholmod_sparse*)A_)->nzmax = nVals_;       /* maximum number of entries */
    ((cholmod_sparse*)A_)->p     = (void*)S.colPtr();   /* column pointers (size n+1) or col indices (size nzmax) */
    ((cholmod_sparse*)A_)->i     = (void*)S.rowIdx();   /* row indices, size nzmax */
    ((cholmod_sparse*)A_)->x     = S.vals();     /* numerical values, size nzmax */

     //std::cout << "CHOLMODWrapper::initialize: " << nVals_ << std::endl;

    ((cholmod_sparse*)A_)->stype  = stype_;

    ((cholmod_sparse*)A_)->itype = CHOLMOD_INT;
    ((cholmod_sparse*)A_)->xtype = CHOLMOD_REAL;
    ((cholmod_sparse*)A_)->dtype = CHOLMOD_DOUBLE;
    ((cholmod_sparse*)A_)->packed = true;
    ((cholmod_sparse*)A_)->sorted = true; // testen, scheint schneller, aber hab ich das immer?
    
    factorise();
      
    return 1;
#else
    std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
  }
  return 0;
}

int CHOLMODWrapper::factorise(){
  if (!dummy_){
#if USE_CHOLMOD
    if (verbose_) cholmod_print_sparse((cholmod_sparse *)A_, "A", (cholmod_common*)c_);
    L_ = cholmod_analyze((cholmod_sparse*)A_, (cholmod_common*)c_);		    /* analyze */
    if (verbose_) std::cout << "Cholmod analyze .. preordering: " << ((cholmod_factor *)(L_))->ordering << std::endl;
    
    cholmod_factorize((cholmod_sparse*)A_, (cholmod_factor *)L_, (cholmod_common*)c_);		    /* factorize */
    
    //    L_ = cholmod_super_symbolic (A_, c_);	/* analyze */
    //    cholmod_super_numeric (A_, L_, c_);		/* factorize */
    
    if (verbose_) cholmod_print_factor((cholmod_factor *)L_, "L", (cholmod_common*)c_);
    return 1;
#else
    std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
  }
  return 0;
}

int CHOLMODWrapper::solve(const RVector & rhs, RVector & solution){
  if (!dummy_){
#if USE_CHOLMOD
    cholmod_dense * b = cholmod_ones(((cholmod_sparse*)A_)->nrow, 1, 
                                     ((cholmod_sparse*)A_)->xtype,
                                     (cholmod_common*)c_);
    cholmod_dense * r = cholmod_zeros(((cholmod_sparse*)A_)->nrow, 1,
                                      ((cholmod_sparse*)A_)->xtype,
                                      (cholmod_common*)c_);
    double * bx = (double*)b->x;
    for (uint i = 0; i < dim_; i++) bx[i] = rhs[i];

    cholmod_dense * x = cholmod_solve(CHOLMOD_A, 
                                      (cholmod_factor *)L_,
                                      b,
                                      (cholmod_common*)c_);	    /* solve Ax=b */
    
    if (((cholmod_sparse*)A_)->stype == 0){
        double al[2] = {0,0}, be[2] = {1,0} ;       /* basic scalars */
        cholmod_sdmult((cholmod_sparse*)A_, 0, be, al, x, r, (cholmod_common*)c_);       
        bx = (double*)r->x; /* ret = Ax */
    } else {
        bx = (double*)x->x; /* ret = x */
    }

    for (uint i = 0; i < dim_; i++) solution[i] = bx[i];
    
    cholmod_free_dense(&r, (cholmod_common*)c_);
    cholmod_free_dense(&x, (cholmod_common*)c_);
    cholmod_free_dense(&b, (cholmod_common*)c_);
    return 1;
#else
    std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
  }
  return 0;
}

int CHOLMODWrapper::solve(const CVector & rhs, CVector & solution){
  if (!dummy_){
#if USE_CHOLMOD
    cholmod_dense * b = cholmod_ones(((cholmod_sparse*)A_)->nrow, 1,
                                     ((cholmod_sparse*)A_)->xtype,
                                     (cholmod_common*)c_);
    cholmod_dense * r = cholmod_ones(((cholmod_sparse*)A_)->nrow, 1,
                                     ((cholmod_sparse*)A_)->xtype,
                                     (cholmod_common*)c_);
    Complex * bx = (Complex*)b->x;
    
    for (uint i = 0; i < dim_; i++) bx[i] = rhs[i];
    
    cholmod_dense * x = cholmod_solve(CHOLMOD_A,                // solve Ax=b
                                      (cholmod_factor *)L_,
                                      b,
                                      (cholmod_common *)c_);     

    if (((cholmod_sparse*)A_)->stype == 0){
        double al[2] = {0,0}, be[2] = {1,0} ;       /* basic scalars */
        cholmod_sdmult((cholmod_sparse*)A_, 0, be, al, x, r, (cholmod_common*)c_);       
        bx = (Complex*)r->x; /* ret = Ax */
    } else {
        bx = (Complex*)x->x; /* ret = x */
    }


    for (uint i = 0; i < dim_; i++) solution[i] = bx[i];
    
    cholmod_free_dense(&x, (cholmod_common*)c_);
    cholmod_free_dense(&b, (cholmod_common*)c_);
    cholmod_free_dense(&r, (cholmod_common*)c_);
    return 1;
#else
    std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
  }
  return 0;
}

    
} //namespace GIMLI
