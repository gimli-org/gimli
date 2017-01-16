/***************************************************************************
 *   Copyright (C) 2007-2014 by the GIMLi development team       *
 *   Carsten Ruecker carsten@resistivity.net                               *
 *   Thomas Guenther thomas@resistivity.net                                *
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

#ifndef _GIMLI_LINSOLVER__H
#define _GIMLI_LINSOLVER__H

#include "gimli.h"
#include <cmath>

namespace GIMLI{

class SolverWrapper;

enum SolverType{AUTOMATIC,LDL,CHOLMOD,UNKNOWN};

class DLLEXPORT LinSolver{
public:
    LinSolver(bool verbose=false);
  
    LinSolver(RSparseMatrix & S, bool verbose=false);
    
    LinSolver(RSparseMapMatrix & S, bool verbose=false);
    
    LinSolver(CSparseMatrix & S, bool verbose=false);

    LinSolver(RSparseMatrix & S, SolverType solverType, bool verbose=false);

    LinSolver(CSparseMatrix & S, SolverType solverType, bool verbose=false);
    
    ~LinSolver();

    void solve(const RVector & rhs, RVector & solution);
    void solve(const CVector & rhs, CVector & solution);
    RVector solve(const RVector & rhs);
    CVector solve(const CVector & rhs);

    void setSolverType(SolverType solverType = AUTOMATIC);
  
    /*! Forwarded to the wrapper to overwrite settings within S. stype =-2 -> use S.stype()*/
    void setMatrix(RSparseMatrix & S, int stype=-2);
    
    /*! Verbose level = -1, use Linsolver.verbose(). */
    void setMatrix(CSparseMatrix & S, int stype=-2);
    
    SolverType solverType() const { return solverType_; }

    std::string solverName() const;
    
protected:
    void init_();
        
    void initialize_(RSparseMatrix & S, int stype);
    void initialize_(CSparseMatrix & S, int stype);
        
    MatrixBase * cacheMatrix_;
    SolverType      solverType_;
    SolverWrapper * solver_;
    bool            verbose_;
    uint rows_;
    uint cols_;
};

template < class Mat, class Vec > int solveLU(const Mat & A, Vec & x, const Vec & b){

	//** from TETGEN

    Mat lu(A);

    int N = 0;
    uint n = b.size();
    int ps[n];

    double scales[n];
    double pivot, biggest, mult, tempf;
    uint pivotindex = 0, tmpIdx = 0;

    for (uint i = N; i < n + N; i++) {
        // Find the largest element in each row for row equilibration
        biggest = 0.0;
        for (uint j = N; j < n + N; j++) {
            if (biggest < (tempf = std::fabs(lu[i][j]))) biggest  = tempf;
        }

        if (biggest != 0.0) scales[i] = 1.0 / biggest;
        else {
            scales[i] = 0.0;
            std::cerr << WHERE_AM_I << " Zero row: singular matrix" << std::endl;
            return false;                      // Zero row: singular matrix.
        }
        ps[i] = i;                                 // Initialize pivot sequence.
    }

    for (uint k = N; k < n + N - 1; k++) {      // For each column.
        // Find the largest element in each column to pivot around.
        biggest = 0.0;
        for (uint i = k; i < n + N; i++) {
            if (biggest < (tempf = std::fabs(lu[ps[i]][k]) * scales[ps[i]])) {
                biggest = tempf;
                pivotindex = i;
            }
        }
        if (biggest == 0.0) {
            std::cerr << WHERE_AM_I << " Zero column: singular matrix" << std::endl;
            return false;                         // Zero column: singular matrix.
        }
        
        if (pivotindex != k) {                         // Update pivot sequence.
            tmpIdx = ps[k];
            ps[k] = ps[pivotindex];
            ps[pivotindex] = tmpIdx;
            //      *d = -(*d);               // ...and change the parity of d.
        }

        // Pivot, eliminating an extra variable  each time
        pivot = lu[ps[k]][k];
    
        for (uint i = k + 1; i < n + N; i++) {
            lu[ps[i]][k] = mult = lu[ps[i]][k] / pivot;
            if (mult != 0.0) {
                for (uint j = k + 1; j < n + N; j++) lu[ps[i]][j] -= mult * lu[ps[k]][j];
            } 
        }
    }

    double dot = 0.0;

    // Vector reduction using U triangular matrix.
    for (uint i = N; i < n + N; i++) {
        dot = 0.0;
        for (uint j = N; j < i + N; j++) dot += lu[ps[i]][j] * x[j];
        x[i] = b[ps[i]] - dot;
    }

    // Back substitution, in L triangular matrix.
    for (int i = n + N - 1; i >= N; i--) {
        dot = 0.0;
        for (uint j = i + 1; j < n + N; j++) dot += lu[ps[i]][j] * x[j];
        x[i] = (x[i] - dot) / lu[ps[i]][i];
    }

    if (rms(Vec(A * x - b)) > 1e-9){
        std::cerr << "rms(A * x -b) " << rms((const Vec)(A * x - b)) << std::endl;
    //     std::cout << Vector(A * x - b) << std::endl;
    //     std::cout << A << b << x << std::endl;
        return -1;
    }
    return 1;
}


} // namespace GIMLI{

#endif //_GIMLI_LINSOLVER__H
