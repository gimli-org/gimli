/***************************************************************************
 *   Copyright (C) 2006-2011 by the resistivity.net development team       *
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

#include "ldlWrapper.h"
#include "vector.h"
#include "sparsematrix.h"

#include <iostream>

#ifdef HAVE_LIBLDL

    #ifndef AMD_INFO
        #define AMD_INFO 20
    #endif

    #ifndef AMD_OK
        #define AMD_OK 0
    #endif

    #ifdef HAVE_LDL_H 
        extern "C"{
            #include <ldl.h>
        }
    #else
        extern "C"{
            extern void ldl_symbolic( int n, int Ap [ ], int Ai [ ], int Lp [ ],
                                      int Parent [ ], int Lnz [ ], int Flag [ ], int P [ ], int Pinv [ ] ) ;

            extern int ldl_numeric( int n, int Ap [ ], int Ai [ ], double Ax [ ],
                                    int Lp [ ], int Parent [ ], int Lnz [ ], int Li [ ], double Lx [ ],
                                    double D [ ], double Y [ ], int Pattern [ ], int Flag [ ],
                                    int P [ ], int Pinv [ ] ) ;

            extern void ldl_lsolve( int n, double X [ ], int Lp [ ], int Li [ ], double Lx [ ] ) ;
            extern void ldl_dsolve( int n, double X [ ], double D [ ]) ;
            extern void ldl_ltsolve( int n, double X [ ], int Lp [ ], int Li [ ],  double Lx [ ] ) ;
            extern void ldl_perm( int n, double X [ ], double B [ ], int P [ ] ) ;
            extern void ldl_permt( int n, double X [ ], double B [ ], int P [ ] ) ;
            extern int ldl_valid_perm( int n, int P [ ], int Flag [ ] ) ;
            extern int ldl_valid_matrix( int n, int Ap [ ], int Ai [ ] ) ;
        }
    #endif // else no HAVE_LDL_H 
    #ifdef HAVE_AMD_H 
        #include <amd.h>
    #else
        extern "C"{
            extern int amd_order( int n, const int Ap [ ], const int Ai [ ], int P [ ], double Control [ ], double Info [ ] ) ;
        }
    #endif // else no HAVE_AMD_H 
#endif  // HAVE_LIBLDL

namespace GIMLI{

#ifdef HAVE_LIBLDL
  bool LDLWrapper::valid(){ return true; }
#else
  bool LDLWrapper::valid(){ return false; }
#endif


LDLWrapper::LDLWrapper( DSparseMatrix & S, bool verbose) : SolverWrapper( S, verbose ) {
  preordering_ = true;
  initialize_( S );
  factorise( );
}

LDLWrapper::~LDLWrapper(){
#ifdef HAVE_LIBLDL
  delete [] Li_;
  delete [] Lp_;
  delete [] Lx_;
  delete [] D_;

  if ( preordering_ ){
    delete [] P_;
    delete [] Pinv_;
  }
  return;
#endif //HAVE_LIBLDL
  std::cerr << WHERE_AM_I << " LDL not installed" << std::endl;
}

int LDLWrapper::initialize_( DSparseMatrix & S ){

#ifdef HAVE_LIBLDL
    colPtr_ = reinterpret_cast< int * >(S.colPtr());
    rowIdx_ = reinterpret_cast< int * >(S.rowIdx());
      THROW_TO_IMPL
     //vals_   = S.vals();

    if ( preordering_ ){
        P_ = new int[ dim_ ];
        Pinv_ = new int[ dim_ ];

#ifdef HAVE_LIBAMD
        double Info[ AMD_INFO ];

//     int * RcolPtr = new int[ dim_ + 1 ];
//     int * RrowIdx = new int[ nVals_ ];
//     amd_preprocess (dim_, AcolPtr_, ArowIdx_, RcolPtr, RrowIdx);
//     AcolPtr_ = RcolPtr;
//     ArowIdx_ = RrowIdx;

        if ( amd_order( dim_, colPtr_, rowIdx_, P_, (double *) NULL, Info) != AMD_OK ){
            throwError( EXIT_AMD_INTERNAL, WHERE_AM_I + " AMD failed for unknown reasons" );
        }
#else
        for ( int i = 0; i < dim_; i ++ ) P_[ i ] = i;
#endif
    } else{
        P_ = (int *) NULL ;
        Pinv_ = (int *) NULL ;
    }

    return 1;
    
#endif //HAVE_LIBLDL
    std::cerr << WHERE_AM_I << " Warning! LDL not installed" << std::endl;
    return 0;
}

int LDLWrapper::factorise(){
#ifdef HAVE_LIBLDL
  Lp_ = new int[ dim_ + 1 ];
  int *Parent = new int[ dim_ ];
  int *Lnz = new int[ dim_ ];
  int *Flag = new int[ dim_ ];

  //  cout << "ldl_valid_matrix: " << ldl_valid_matrix ( dim_, colPtr_, rowIdx_) << std::endl;

  ldl_symbolic( dim_, colPtr_, rowIdx_, Lp_, Parent, Lnz, Flag, P_, Pinv_);
  int lnz = Lp_[ dim_ ];

  Li_ = new int[ lnz ];
  Lx_ = new double [ lnz ];
  D_ = new double[ dim_ ];

  double *Y = new double[ dim_ ];
  int *Pattern = new int[ dim_ ];
  //  for ( int i = 0; i < dim_; i ++ ) std::cout << vals_[ i ] << std::endl;

  uint d = ldl_numeric( dim_, colPtr_, rowIdx_, vals_, Lp_, Parent, Lnz, Li_, Lx_, D_, Y, Pattern, Flag, P_, Pinv_) ;

  if ( d != dim_ ) {
    std::cout << WHERE_AM_I << " Warning! ldl_numeric failed, D (" << d << "," << d << ") is zero\n";
  }

  delete [] Y;
  delete [] Parent;
  delete [] Lnz;
  delete [] Flag;
  delete [] Pattern;

  return 1;
#endif //HAVE_LIBLDL
  std::cerr << WHERE_AM_I << " Warning! LDL not installed" << std::endl;
  return 0;
}

int LDLWrapper::solve( const RVector & rhs, RVector & solution ){
#ifdef HAVE_LIBLDL
  if ( solution.size() != dim_ ) solution.resize( dim_ );

  if ( rhs.size() == dim_ ){
    RVector b( rhs );

    THROW_TO_IMPL
//     if ( preordering_ ){
//       RVector bP( dim_ );
//       ldl_perm ( dim_, &bP[ 0 ], &b[ 0 ], P_ );			/* y = Pb */
//       ldl_lsolve( dim_, &bP[ 0 ], Lp_, Li_, Lx_ );
//       ldl_dsolve( dim_, &bP[ 0 ], D_ ) ;
//       ldl_ltsolve( dim_, &bP[ 0 ], Lp_, Li_, Lx_ );
//       ldl_permt ( dim_, &solution[ 0 ], &bP[ 0 ], P_ );		/* x = P'y */
//     } else {
//       ldl_lsolve( dim_, &b[ 0 ], Lp_, Li_, Lx_ ) ;
//       ldl_dsolve( dim_, &b[ 0 ], D_ ) ;
//       ldl_ltsolve( dim_, &solution[ 0 ], Lp_, Li_, Lx_ ) ;
//     }
  } else {
    std::cerr << WHERE_AM_I << " rhs is not defined" << std::endl;
  }

  return 1;
#endif //HAVE_LIBLDL
  std::cerr << WHERE_AM_I << " Warning! LDL not installed" << std::endl;
  return 0;
}

} //namespace GIMLI
