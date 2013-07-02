/***************************************************************************
 *   Copyright (C) 2007-2013 by the resistivity.net development team       *
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

#ifdef CHOLMOD_FOUND
    #if CHOLMOD_FOUND==TRUE
        #define USE_CHOLMOD 1
    #endif
#endif

#ifdef USE_CHOLMOD
    #define CHOLMOD_MAXMETHODS 9
    #define UF_long long

    #include <cholmod.h>
    #define USE_CHOLMOD 1

#elif HAVE_LIBCHOLMOD

extern "C" {
    #define USE_CHOLMOD
    #define CHOLMOD_MAXMETHODS 9
    #define UF_long long

  struct cholmod_common{
    double dbound ;
    double grow0 ;
    double grow1 ;
    size_t grow2 ;
    size_t maxrank ;
    double supernodal_switch ;
    int supernodal ;
    int final_asis ;
    int final_super ;
    int final_ll ;
    int final_pack ;
    int final_monotonic ;
    int final_resymbol ;
    double zrelax [3] ;
    size_t nrelax [3] ;
    int prefer_zomplex ;
    int prefer_upper ;
    int quick_return_if_not_posdef ;
    int print ;
    int precise ;
    int (*print_function) (const char *, ...) ;
    int try_catch ;
    void (*error_handler) (int status, char *file, int line, char *message) ;
    int nmethods ;
    int current ;
    int selected ;
    struct cholmod_method_struct    {
	double lnz ;
	double fl ;
	double prune_dense ;
	double prune_dense2 ;
	double nd_oksep ;
	double other1 [4] ;
	size_t nd_small ;
	size_t other2 [4] ;
	int aggressive ;
	int order_for_lu ;
	int nd_compress ;
	int nd_camd ;
	int nd_components ;
	int ordering ;
	size_t other3 [4] ;
    } method [CHOLMOD_MAXMETHODS + 1] ;

    int postorder ;
    void *(*malloc_memory) (size_t) ;
    void *(*realloc_memory) (void *, size_t) ;
    void (*free_memory) (void *) ;
    void *(*calloc_memory) (size_t, size_t) ;
    int (*complex_divide) (double ax, double az, double bx, double bz,
	    double *cx, double *cz) ;
    double (*hypotenuse) (double x, double y) ;
    double metis_memory ;
    double metis_dswitch ;
    size_t metis_nswitch ;
    size_t nrow ;
    UF_long mark ;
    size_t iworksize ;
    size_t xworksize ;
    void *Flag ;
    void *Head ;
    void *Xwork ;
    void *Iwork ;
    int itype ;
    int dtype ;
    int no_workspace_reallocate ;
    int status ;
    double fl ;
    double lnz ;
    double anz ;
    double modfl ;
    size_t malloc_count ;
    size_t memory_usage ;
    size_t memory_inuse ;
    double nrealloc_col ;
    double nrealloc_factor ;
    double ndbounds_hit ;
    double rowfacfl ;
    double aatfl ;
    double  other1 [16] ;
    UF_long other2 [16] ;
    int     other3 [13] ;
    int prefer_binary ;
    int default_nesdis ;
    int called_nd ;
    size_t  other4 [16] ;
    void   *other5 [16] ;
};

  struct cholmod_dense  {
    size_t nrow ;
    size_t ncol ;
    size_t nzmax ;
    size_t d ;
    void *x ;
    void *z ;
    int xtype ;
    int dtype ;
  };

  struct cholmod_sparse {
    size_t nrow ;
    size_t ncol ;
    size_t nzmax ;
    void *p ;
    void *i ;
    void *nz ;
    void *x ;
    void *z ;
    int stype ;
    int itype ;
    int xtype ;
    int dtype ;
    int sorted ;
    int packed ;
  };

#define CHOLMOD_A    0    /* solve Ax=b */
#define CHOLMOD_INT  0    /* all integer arrays are int */
#define CHOLMOD_REAL 1    /* a real matrix */
#define CHOLMOD_DOUBLE 0  /* all numerical values are double */

  int cholmod_start( cholmod_common * Common );
  int cholmod_finish( cholmod_common * Common );
  cholmod_factor * cholmod_analyze( cholmod_sparse * A,	cholmod_common * Common );
  int cholmod_factorize( cholmod_sparse *A, cholmod_factor *L, cholmod_common *Common );
  cholmod_dense * cholmod_ones( size_t nrow, size_t ncol, int xtype, cholmod_common *Common );
  cholmod_dense * cholmod_solve( int, cholmod_factor * L, cholmod_dense *b, cholmod_common *Common );
  int cholmod_free_dense( cholmod_dense **b, cholmod_common *Common ) ;

  cholmod_sparse * cholmod_allocate_sparse (size_t, size_t, size_t, int, int, int, int, cholmod_common *) ;
  int cholmod_free_sparse( cholmod_sparse ** A, cholmod_common *Common ) ;
  int cholmod_free_factor( cholmod_factor ** L, cholmod_common *Common ) ;
}

#endif

namespace GIMLI{


#ifdef USE_CHOLMOD
  bool CHOLMODWrapper::valid() { return true; }
#else
  bool CHOLMODWrapper::valid() { return false; }
#endif

CHOLMODWrapper::CHOLMODWrapper( DSparseMatrix & S, bool verbose ) : SolverWrapper( S, verbose ){
  c_ = NULL;
  A_ = NULL;
  L_ = NULL;
#ifdef USE_CHOLMOD
  c_ = new cholmod_common;
//   cholmod_common *c_;
//   cholmod_sparse *A_;
//   cholmod_factor *L_;
  
  int ret =  cholmod_start( (cholmod_common*)c_ );
  if ( ret ) dummy_ = false;

  initialize_( S );
  //  dummy_ = true;
  factorise();
#else

  std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
}

CHOLMODWrapper::~CHOLMODWrapper(){
#ifdef USE_CHOLMOD
  cholmod_free_factor( (cholmod_factor**)(&L_), (cholmod_common*)c_ );
  //** We did not allocate the matrix so we dont need to free it
  //  cholmod_free_sparse( &A_, c_ );
  cholmod_finish( (cholmod_common*)c_ );

  if ( A_ ) delete (cholmod_sparse*)A_;
  if ( c_ ) delete (cholmod_common*)c_;

#else
  std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
}

int CHOLMODWrapper::factorise(){
  if ( !dummy_ ){
#ifdef USE_CHOLMOD
    L_ = cholmod_analyze( (cholmod_sparse*)A_, (cholmod_common*)c_  );		    /* analyze */
    cholmod_factorize( (cholmod_sparse*)A_, (cholmod_factor *)L_, (cholmod_common*)c_ );		    /* factorize */
    //    L_ = cholmod_super_symbolic (A_, c_ );	/* analyze */
    //    cholmod_super_numeric (A_, L_, c_);		/* factorize */
    return 1;
#else
    std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
  }
  return 0;
}

int CHOLMODWrapper::solve( const RVector & rhs, RVector & solution ){
  if ( !dummy_ ){
#ifdef USE_CHOLMOD
    cholmod_dense * b = cholmod_ones( ((cholmod_sparse*)A_)->nrow, 1, ((cholmod_sparse*)A_)->xtype, (cholmod_common*)c_ );
    double * bx = (double*)b->x;
    for ( uint i = 0; i < dim_; i++) bx[ i ] = rhs[ i ];

    cholmod_dense * x = cholmod_solve( CHOLMOD_A, (cholmod_factor *)L_, b, (cholmod_common*)c_ );	    /* solve Ax=b */
    bx = (double*)x->x;

    for ( uint i = 0; i < dim_; i++) solution[ i ] = bx[ i ];
    cholmod_free_dense( &x, (cholmod_common*)c_ ) ;
    cholmod_free_dense( &b, (cholmod_common*)c_ ) ;
    return 1;
#else
    std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
  }
  return 0;
}

int CHOLMODWrapper::initialize_( DSparseMatrix & S ){
  if ( !dummy_ ){
#ifdef USE_CHOLMOD

  //** We do not allocate the matrix since we use the allocated space from DSparsemarix
//    A_ = cholmod_allocate_sparse( dim_, dim_, nVals_, true, true, 1, CHOLMOD_REAL, c_ ) ;

    A_ = new cholmod_sparse;
    ((cholmod_sparse*)A_)->nrow  = dim_;         /* number of rows */
    ((cholmod_sparse*)A_)->ncol  = dim_;	       /* number of columns */
    ((cholmod_sparse*)A_)->nzmax = nVals_;       /* maximum number of entries */
    ((cholmod_sparse*)A_)->p     = (void*)S.colPtr();   /* column pointers (size n+1) or col indices (size nzmax) */
    ((cholmod_sparse*)A_)->i     = (void*)S.rowIdx();   /* row indices, size nzmax */
    ((cholmod_sparse*)A_)->x     = S.vals();     /* numerical values, size nzmax */

     //std::cout << "CHOLMODWrapper::initialize: " << nVals_ << std::endl;

    ((cholmod_sparse*)A_)->stype  = 1;

    ((cholmod_sparse*)A_)->itype = CHOLMOD_INT;
    ((cholmod_sparse*)A_)->xtype = CHOLMOD_REAL;
    ((cholmod_sparse*)A_)->dtype = CHOLMOD_DOUBLE;
    ((cholmod_sparse*)A_)->packed = true;
    ((cholmod_sparse*)A_)->sorted = true; // testen, scheint schneller, aber hab ich das immer?
    return 1;
#else
    std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
  }
  return 0;
}

} //namespace GIMLI
