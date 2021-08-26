/******************************************************************************
 *   Copyright (C) 2006-2021 by the GIMLi development team                    *
 *   Carsten Rücker carsten@resistivity.net                                   *
 *                                                                            *
 *   Licensed under the Apache License, Version 2.0 (the "License");          *
 *   you may not use this file except in compliance with the License.         *
 *   You may obtain a copy of the License at                                  *
 *                                                                            *
 *       http://www.apache.org/licenses/LICENSE-2.0                           *
 *                                                                            *
 *   Unless required by applicable law or agreed to in writing, software      *
 *   distributed under the License is distributed on an "AS IS" BASIS,        *
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 *   See the License for the specific language governing permissions and      *
 *   limitations under the License.                                           *
 *                                                                            *
 ******************************************************************************/

#pragma once

#include "gimli.h"
#include "elementmatrix.h"
#include "vector.h"
#include "matrix.h"

namespace GIMLI{

class DLLEXPORT ElementMatrixMap {
public:
    
    void push_back(const ElementMatrix < double > & Ai);
    
    void resize(Index size);

    #define DEFINE_INTEGRATOR(A_TYPE) \
        /*! R = \int_mesh this * f \d mesh and R = RVector(dof) and \
            f = A_TYPE \
        */ \
        void integrate(const A_TYPE & f, RVector & R, bool neg=false) const; \
        /*! Integrate into bilinear form R = \int_mesh this * f * R \d mesh an\
        R = RSparseMapMatrix(dof, dof) and \
            f = A_TYPE \
        */ \
        RVector integrate(const A_TYPE & f, bool neg=false) const; \
        /*! R = this * f */ \
        void mult(const A_TYPE & f, ElementMatrixMap & ret) const; 

    DEFINE_INTEGRATOR(double)   // const scalar for all cells
    DEFINE_INTEGRATOR(RMatrix)  // const Matrix for all cells
    DEFINE_INTEGRATOR(RVector)      // const scalar for each cells
    DEFINE_INTEGRATOR(Pos)      // const vector for all cells //!calling order!
    DEFINE_INTEGRATOR(PosVector)    // const vector for each cells
    DEFINE_INTEGRATOR(std::vector< RMatrix >)// const matrix for each cells
    DEFINE_INTEGRATOR(std::vector< RVector >)// scalar for quadr. on each cells
    DEFINE_INTEGRATOR(std::vector< PosVector >)// vector for quadr. on each cells
    DEFINE_INTEGRATOR(std::vector< std::vector< RMatrix > >)// mat for quadr. on each cells

    #undef DEFINE_INTEGRATOR

    #define DEFINE_INTEGRATOR(A_TYPE) \
        /*! Integrate into bilinear form R = \int_mesh this * f * R \d mesh an\
        R = RSparseMapMatrix(dof, dof) and \
            f = A_TYPE \
        */ \
        void integrate(const ElementMatrixMap & R, const A_TYPE & f, \
                       SparseMatrixBase & A, bool neg=false) const; \
        RSparseMapMatrix integrate(const ElementMatrixMap & R, \
                                   const A_TYPE & f, bool neg=false) const; \

    DEFINE_INTEGRATOR(double)   // const scalar for all cells
    DEFINE_INTEGRATOR(RMatrix)  // const Matrix for all cells
    DEFINE_INTEGRATOR(RVector)      // const scalar for each cells
    DEFINE_INTEGRATOR(std::vector< RMatrix >)// const matrix for each cells
    DEFINE_INTEGRATOR(std::vector< RVector >)// scalar for quadr. on each cells
    DEFINE_INTEGRATOR(std::vector< std::vector< RMatrix > >)// mat for quadr. on each cells

    #undef DEFINE_INTEGRATOR

    const std::vector< ElementMatrix < double > > & mats() const;
    
    ElementMatrix < double > * pMat(Index i){ return & mats_[i]; }

    const std::vector < PosVector > & quadraturePoints() const;

    /*!Calculate copy of this + B, depending on requested dim. */
    void add(const ElementMatrixMap & B, ElementMatrixMap & ret, Index dim=1, double b=1.0) const;
    
    void add(Index row, const ElementMatrix < double > & Ai);

    //TODO .. check if its the same like mult(a-b, m-n))
    RVector mult(const RVector & a, const RVector & b,
                 const RVector & m, const RVector & n) const;

    /*! Return (S_i * a) * b for all i*/
    RVector mult(const RVector & a, const RVector & b) const;

    Index size() const { return mats_.size();}

    Index rows() const { return rows_; }

    Index cols() const { return cols_; }

protected:
    std::vector< ElementMatrix < double > > mats_;
    mutable std::vector < PosVector > quadrPnts_; // cache Q for mats_

    std::vector< SmallMatrix > mat_;
    std::vector< IndexArray > _ids;
    std::vector< Index > row_;

    Index rows_;
    Index cols_;
};


DLLEXPORT void
createUMap(const Mesh & mesh, Index order,
           ElementMatrixMap & ret,
           Index nCoeff=1, Index dofOffset=0);

DLLEXPORT ElementMatrixMap 
createUMap(const Mesh & mesh, Index order,
           Index nCoeff=1, Index dofOffset=0);

DLLEXPORT void 
createdUMap(const Mesh & mesh, Index order,
            ElementMatrixMap & ret,
            bool elastic, bool div, bool kelvin,
            Index nCoeff=1, Index dofOffset=0);

DLLEXPORT ElementMatrixMap
createdUMap(const Mesh & mesh, Index order,
            bool elastic, bool div, bool kelvin,
            Index nCoeff=1, Index dofOffset=0);

DLLEXPORT void 
createIdentityMap(const Mesh & mesh, Index order, 
                  ElementMatrixMap & ret,
                  Index nCoeff=1, Index dofOffset=0);

DLLEXPORT ElementMatrixMap 
createIdentityMap(const Mesh & mesh, Index order,
                  Index nCoeff, Index dofOffset);

DLLEXPORT ElementMatrixMap sym(const ElementMatrixMap & A);
DLLEXPORT void sym(const ElementMatrixMap & A, ElementMatrixMap & ret);

DLLEXPORT RVector tr(const ElementMatrixMap & A);
DLLEXPORT void tr(const ElementMatrixMap & A, RVector & ret);


} // namespace GIMLI{
