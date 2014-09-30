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
    
    #if UMFPACK_FOUND
        #include <umfpack.h>
        #define USE_UMFPACK 1
    #endif
#endif

namespace GIMLI{

#if USE_CHOLMOD
    bool CHOLMODWrapper::valid() { return true; }
#else
    bool CHOLMODWrapper::valid() { return false; }
#endif

CHOLMODWrapper::CHOLMODWrapper(RSparseMatrix & S, bool verbose, int stype)
    : SolverWrapper(S, verbose){
    init_(S, stype);
}

CHOLMODWrapper::CHOLMODWrapper(CSparseMatrix & S, bool verbose, int stype)
    : SolverWrapper(S, verbose){
    init_(S, stype);
}

CHOLMODWrapper::~CHOLMODWrapper(){
#if USE_CHOLMOD

    if (L_) cholmod_free_factor((cholmod_factor**)(&L_), (cholmod_common*)c_);
//     ** We did not allocate the matrix so we dont need to free it
//      cholmod_free_sparse(&A_, c_);
    
    cholmod_finish((cholmod_common*)c_);

    if (A_) delete (cholmod_sparse*)A_;
    if (c_) delete (cholmod_common*)c_;
  
    
    if (Numeric_) umfpack_di_free_numeric (&Numeric_);
    if (AxV_) delete AxV_;
    if (AzV_) delete AzV_;
   
#else
    std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
}

template < class ValueType > 
void CHOLMODWrapper::init_(SparseMatrix < ValueType > & S, int stype){

    useUmfpack_ = false;
    Numeric_ = 0;
    c_ = NULL;
    A_ = NULL;
    L_ = NULL;
    AxV_ = NULL;
    AzV_ = NULL;
    
    if (stype == -2){
        stype_ = S.stype();
    } else {
        stype_ = stype;
    }
    
#if USE_CHOLMOD
    int ret = 0;
    c_ = new cholmod_common;
    ret =  cholmod_start((cholmod_common*)c_);
    if (ret) dummy_ = false;

    initializeMatrix_(S);
#else
    std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
}

template < class ValueType > 
int CHOLMODWrapper::initMatrixChol_(SparseMatrix < ValueType > & S, int xType){
    if (!dummy_){
#if USE_CHOLMOD
        A_ = new cholmod_sparse;
        ((cholmod_sparse*)A_)->nrow  = S.nRows();         /* number of rows */
        ((cholmod_sparse*)A_)->ncol  = S.nCols();           /* number of columns */
        ((cholmod_sparse*)A_)->nzmax = S.nVals();               /* maximum number of entries */
        ((cholmod_sparse*)A_)->p     = (void*)S.colPtr();   /* column pointers (size n+1) or col indices (size nzmax) */
        ((cholmod_sparse*)A_)->i     = (void*)S.rowIdx();   /* row indices, size nzmax */
    
        ((cholmod_sparse*)A_)->x     = S.vals();     /* numerical values, size nzmax */
        ((cholmod_sparse*)A_)->stype = stype_;

        ((cholmod_sparse*)A_)->itype = CHOLMOD_INT;
        ((cholmod_sparse*)A_)->xtype = xType;  // data type for the pattern (Real, complex, zcomplex)
        ((cholmod_sparse*)A_)->dtype = CHOLMOD_DOUBLE; // data type for complex or real (float/double)
        ((cholmod_sparse*)A_)->packed = true;
        ((cholmod_sparse*)A_)->sorted = true; // testen, scheint schneller, aber hab ich das immer?
    
        factorise();
#else
        std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
    } // if ! dummy
    return 0;
}
    
int CHOLMODWrapper::initializeMatrix_(CSparseMatrix & S){
    if (!dummy_){
        // check for non-hermetian 
        if (S.stype() == 0){ //  matrix is symmetric
            for (Index i = 0; i < S.size(); i++){
                for (int j = S.vecColPtr()[i]; j < S.vecColPtr()[i + 1]; j ++){
                    if (std::imag(S.vecVals()[j]) != 0.0){
                        if (S.vecVals()[j] == S.getVal(S.vecRowIdx()[j], i)){
                            // non-hermetian symmetric
                            if (verbose_) std::cout << "non-hermetian symmetric matrix found .. switching to umfpack." << std::endl;
                            useUmfpack_ = true;
                        }
                        i = S.size();
                        break;
                    }
                }
            }
        }
     
        if (useUmfpack_){
#if USE_UMFPACK
            Ap_ = (int*)S.colPtr();
            Ai_ = (int*)S.rowIdx();
            AxV_ = new RVector(real(S.vecVals()));
            AzV_ = new RVector(imag(S.vecVals()));
            double * Ax_ = &((*AxV_)[0]);
            double * Az_ = &((*AzV_)[0]);
        
            double *null = (double *) NULL;
            void *Symbolic;
        
            if (verbose_) std::cout << "Using umfpack .. " << std::endl;        
        
            umfpack_zi_symbolic (S.nRows(), S.nRows(), Ap_, Ai_, Ax_, Az_, &Symbolic, null, null) ;
            umfpack_zi_numeric (Ap_, Ai_, Ax_, Az_, Symbolic, &Numeric_, null, null) ;
            umfpack_zi_free_symbolic (&Symbolic);
            return 1;
#else
        std::cerr << WHERE_AM_I << " umfpack not installed" << std::endl;
#endif
        } else {
            return initMatrixChol_(S, CHOLMOD_COMPLEX);
        }
    } // if ! dummy
    return 0;
}

int CHOLMODWrapper::initializeMatrix_(RSparseMatrix & S){
    return initMatrixChol_(S, CHOLMOD_REAL);
}

int CHOLMODWrapper::factorise(){
    if (!dummy_){
#if USE_CHOLMOD
        if (verbose_) cholmod_print_sparse((cholmod_sparse *)A_, "A", (cholmod_common*)c_);

        L_ = cholmod_analyze((cholmod_sparse*)A_,
                         (cholmod_common*)c_);		    /* analyze */
        cholmod_factorize((cholmod_sparse*)A_,
                      (cholmod_factor*)L_,
                      (cholmod_common*)c_);		    /* factorize */
    
    if (verbose_) std::cout << "Cholmod analyze .. preordering: " << ((cholmod_factor *)(L_))->ordering << std::endl;
    if (verbose_) cholmod_print_factor((cholmod_factor *)L_, "L", (cholmod_common*)c_);
    return 1;
#else
        std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
    }
    return 0;
}

template < class ValueType > 
    int CHOLMODWrapper::solveCHOL_(const Vector < ValueType > & rhs, Vector < ValueType > & solution){
    if (!dummy_){
#if USE_CHOLMOD
    
        cholmod_dense * b = cholmod_ones(((cholmod_sparse*)A_)->nrow, 1, 
                                        ((cholmod_sparse*)A_)->xtype,
                                        (cholmod_common*)c_);
        cholmod_dense * r = cholmod_zeros(((cholmod_sparse*)A_)->nrow, 1,
                                        ((cholmod_sparse*)A_)->xtype,
                                        (cholmod_common*)c_);
        ValueType * bx = (ValueType*)b->x;
        for (uint i = 0; i < dim_; i++) bx[i] = rhs[i];

        cholmod_dense * x = cholmod_solve(CHOLMOD_A, 
                                        (cholmod_factor *)L_,
                                        b,
                                        (cholmod_common*)c_);       /* solve Ax=b */
    
        if (((cholmod_sparse*)A_)->stype == 0){
            double al[2] = {0,0}, be[2] = {1,0} ;       /* basic scalars */
            cholmod_sdmult((cholmod_sparse*)A_, 0, be, al, x, r, (cholmod_common*)c_);       
            bx = (ValueType *)r->x; /* ret = Ax */
            for (uint i = 0; i < dim_; i++) solution[i] = conj(bx[i]);
        } else {
            bx = (ValueType *)x->x; /* ret = x */
            for (uint i = 0; i < dim_; i++) solution[i] = bx[i];
        }
    
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
    
int CHOLMODWrapper::solve(const RVector & rhs, RVector & solution){
    return solveCHOL_(rhs, solution);
}

int CHOLMODWrapper::solve(const CVector & rhs, CVector & solution){
    if (!dummy_){

        if (useUmfpack_){
#if USE_UMFPACK
            double * Ax_ = &((*AxV_)[0]);
            double * Az_ = &((*AzV_)[0]);
        
            double *null = (double *) NULL ;
        
            RVector xx(rhs.size());
            RVector xz(rhs.size());
            RVector bx(real(rhs));
            RVector bz(imag(rhs));
        
            umfpack_zi_solve(UMFPACK_A, Ap_, Ai_, Ax_, Az_,
                             &xx[0], &xz[0], &bx[0], &bz[0],
                             Numeric_, null, null) ;
        
            solution = toComplex(xx, xz);
            return 1;
#endif
        } else {
            return solveCHOL_(rhs, solution);
        }
    }
    return 0;
}

    
} //namespace GIMLI
