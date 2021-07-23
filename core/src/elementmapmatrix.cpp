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

#include "elementmapmatrix.h"
#include "shape.h"
#include "meshentities.h"
#include "node.h"
#include "pos.h"
#include "sparsematrix.h"

#include "integration.h"

namespace GIMLI{

void ElementMatrixMap::resize(Index size){
    mats_.resize(size);
    rows_ = this->mats_.size();
}

void ElementMatrixMap::push_back(const ElementMatrix < double > & Ai){
    mats_.push_back(Ai);
    rows_ = this->mats_.size();
}

template < class ValueType >
void integrateLConstT_(const ElementMatrixMap * self, 
                     const ValueType & f, RVector & R, bool neg){
    RVector rt;
    for (auto &m : self->mats()){
        m.integrate(f, rt);
        
        if (neg) {
            R.subVal(rt, m.rowIDs());
        } else {
            R.addVal(rt, m.rowIDs());
        }
    }
}

template < class ValueType >
void integrateLPerCellT_(const ElementMatrixMap * self, 
                        const ValueType & f, RVector & R, bool neg){
    ASSERT_EQUAL_SIZE(self->mats(), f)
    
    RVector rt;
    for (auto &m : self->mats()){
        m.integrate(f[m.entity()->id()], rt);
        
        if (neg) {
            R.subVal(rt, m.rowIDs());
        } else {
            R.addVal(rt, m.rowIDs());
        }
    }
}


// //** R*f space f=RVector at quadrature points per cell
// void ElementMatrixMap::integrate(const std::vector< RVector > & f,
//                                  RVector & R, bool neg) const{
//     ASSERT_EQUAL_SIZE(this->mats_, f)
//     for (auto &m : this->mats_){
//         const RVector &w(*m.w());
//         const RVector &fi(f[m.entity()->id()]);

//         Index nRules(w.size());

//         // __MS(w)
//         // __MS(b)
//         ASSERT_VEC_SIZE(fi, nRules)
//         ASSERT_VEC_SIZE(m.matX(), nRules)

//         RVector rt(m.rows(), 0.0);
// #if USE_EIGEN3
//                     __MS("EIGENCHECK")
// #endif
//         for (Index r = 0; r < nRules; r++){
//             // __MS(m.matX()[r], b[r])
//             rt += m.matX()[r](0) * fi[r] * w[r];
//         }
//         rt *= m.entity()->size();
//         if (neg) R.subVal(rt, m.rowIDs());
//         else R.addVal(rt, m.rowIDs());
//     }
// }

// //** R*f Vector space f=PosVector at quadrature points per cell
// void ElementMatrixMap::integrate(const std::vector< PosVector > & f,
//                                  RVector & R, bool neg) const{
//     ASSERT_EQUAL_SIZE(this->mats_, f)
//     for (auto &m : this->mats_){

//         const RVector &w(*m.w());
//         const PosVector &b(f[m.entity()->id()]);

//         Index nRules(w.size());

//         // __MS(w)
//         // __MS(b)
//         ASSERT_VEC_SIZE(b, nRules)
//         ASSERT_VEC_SIZE(m.matX(), nRules)

//         RVector rt(m.rows(), 0.0);

//         for (Index r = 0; r < nRules; r++){
//             const SmallMatrix &mr(m.matX()[r]);

//             for (Index k = 0; k < mr.rows(); k ++){
//                 // __MS(r << " " << k << " " << b[r][k])
//                 #if USE_EIGEN3
//                     __MS("EIGENCHECK")
//                 #endif
//                 rt += mr(k) * b[r][k]* w[r];
//             }
//             //rt *= w[r];
//         }

//         rt *= m.entity()->size();
//         if (neg) R.subVal(rt, m.rowIDs());
//         else R.addVal(rt, m.rowIDs());
//     }
// }

#define DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL(A_TYPE) \
void ElementMatrixMap::integrate(const A_TYPE & f, RVector & R, bool neg) const { integrateLConstT_(this, f, R, neg); }\

DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL(double)
DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL(Pos)
DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL(RMatrix)
#undef DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL

#define DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL(A_TYPE) \
void ElementMatrixMap::integrate(const A_TYPE & f, RVector & R, bool neg) const { integrateLPerCellT_(this, f, R, neg); }\

DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL(RVector)
DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL(PosVector)
DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL(std::vector< RMatrix >)
DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL(std::vector< RVector >) 
DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL(std::vector< PosVector >)
DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL(std::vector< std::vector< RMatrix > >)
#undef DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL

#define DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL_RET(A_TYPE) \
RVector ElementMatrixMap::integrate(const A_TYPE & f, bool neg) const { \
    Index maxR = 0; \
    for (auto &m : this->mats_){ \
        maxR = max(maxR, max(m.rowIDs()));\
    } \
    RVector R(maxR+1); \
    integrate(f, R, neg); \
    return R; \
}\
RSparseMapMatrix ElementMatrixMap::integrate(const ElementMatrixMap & R, \
                                             const A_TYPE & f, bool neg) const { \
    RSparseMapMatrix A(0,0);\
    integrate(R, f, A, neg);\
    return A;\
}\

DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL_RET(double)
DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL_RET(Pos)
DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL_RET(RMatrix)
DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL_RET(RVector)
DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL_RET(PosVector)
DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL_RET(std::vector< RMatrix >)
DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL_RET(std::vector< RVector >)
DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL_RET(std::vector< PosVector >)
DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL_RET(std::vector< std::vector< RMatrix > >)
#undef DEFINE_INTEGRATE_ELEMENTMAP_R_IMPL_RET


// M * scalar * M
void ElementMatrixMap::integrate(const ElementMatrixMap & R,
                                 const double & f,
                                 SparseMatrixBase & A, bool neg) const {
    ASSERT_EQUAL_SIZE(this->mats_, R.mats())

    Index i = 0;
    ElementMatrix < double > dAB;
    double scale = 1.0;
    if (neg) scale = -1.0;

    for (auto &m : this->mats_){
        //dot(m, R.mats()[i], f, A, scale);
        dot(m, R.mats()[i], f, dAB);
        A.add(dAB, scale);
        i++;
    }
}
// M * Vector * M
void ElementMatrixMap::integrate(const ElementMatrixMap & R,
                                 const Pos & f,
                                 SparseMatrixBase & A, bool neg) const {
    ASSERT_EQUAL_SIZE(this->mats_, R.mats())

    Index i = 0;
    ElementMatrix < double > dAB;
    ElementMatrix < double > mf;
    double scale = 1.0;
    if (neg) scale = -1.0;

    for (auto &m : this->mats_){
        // refactor without mult
        GIMLI::mult(m, f, mf);

        //dot(m, R.mats()[i], f, A, scale);
        dot(mf, R.mats()[i], 1, dAB);
        A.add(dAB, scale);
        i++;
    }
}
// M * scalar(cell) * M
void ElementMatrixMap::integrate(const ElementMatrixMap & R,
                                 const RVector & f,
                                 SparseMatrixBase & A, bool neg) const {
    
    
    ASSERT_EQUAL_SIZE(this->mats_, f)
    ASSERT_EQUAL_SIZE(this->mats_, R.mats())

    Index i = 0;
    double scale = 1.0;
    if (neg) scale = -1;
    ElementMatrix < double > dAB;
    
    for (auto &m : this->mats_){
        dot(m, R.mats()[i], f[m.entity()->id()], dAB);
        
        A.add(dAB, scale);
        i++;
    }
}

void ElementMatrixMap::integrate(const ElementMatrixMap & R,
                                 const std::vector< RVector > & f,
                                 SparseMatrixBase & A, bool neg) const {
    ASSERT_EQUAL_SIZE(this->mats_, f)
    ASSERT_EQUAL_SIZE(this->mats_, R.mats())

    ElementMatrix < double > dAB;
    ElementMatrix < double > mf;
    
    Index i = 0;
    double scale = 1.0;
    if (neg) scale = -1;

    for (auto &m : this->mats_){
        const RVector &fi(f[m.entity()->id()]);
        // refactor without mult
        GIMLI::mult(m, fi, mf);

        // void mult(const ElementMatrix < double > & A, const RVector & b,
        //           ElementMatrix < double > & C){

        dot(mf, R.mats()[i], 1.0, dAB);
        
        A.add(dAB, scale);
        i++;
    }
}

// bilinear forms R*f*R
#define DEFINE_INTEGRATE_ELEMENTMAP_A_IMPL(A_TYPE) \
void ElementMatrixMap::integrate(const ElementMatrixMap & R, const A_TYPE & f, \
                                 SparseMatrixBase & A, bool neg) const {\
    THROW_TO_IMPL \
} \

//DEFINE_INTEGRATE_ELEMENTMAP_A_IMPL(double)
//DEFINE_INTEGRATE_ELEMENTMAP_A_IMPL(Pos)
DEFINE_INTEGRATE_ELEMENTMAP_A_IMPL(RMatrix)
// DEFINE_INTEGRATE_ELEMENTMAP_A_IMPL(RVector)
DEFINE_INTEGRATE_ELEMENTMAP_A_IMPL(PosVector)
DEFINE_INTEGRATE_ELEMENTMAP_A_IMPL(std::vector< RMatrix >)
//DEFINE_INTEGRATE_ELEMENTMAP_A_IMPL(std::vector< RVector >)
DEFINE_INTEGRATE_ELEMENTMAP_A_IMPL(std::vector< PosVector >)
DEFINE_INTEGRATE_ELEMENTMAP_A_IMPL(std::vector< std::vector< RMatrix > >)
#undef DEFINE_INTEGRATE_ELEMENTMAP_A_IMPL

const std::vector< ElementMatrix < double > > & ElementMatrixMap::mats() const{
    return mats_;
}

const std::vector < PosVector > & ElementMatrixMap::quadraturePoints() const {

    if (this->quadrPnts_.size() != this->mats_.size()){
        this->quadrPnts_.clear();

        for (auto &m: this->mats_){
            const PosVector &x(*m.x());
            this->quadrPnts_.push_back(PosVector(x.size()));
            for (Index i = 0; i < x.size(); i ++){
                this->quadrPnts_.back()[i] = m.entity()->shape().xyz(x[i]);
            }
        }
    }
    return this->quadrPnts_;
}

void ElementMatrixMap::add(Index row, const ElementMatrix < double > & Ai){

    rows_ = max(row + 1, rows_);
    cols_ = max(max(Ai.ids()) + 1, cols_);

    mat_.push_back(Ai.mat());
    _ids.push_back(Ai.ids());
    row_.push_back(row);
}

RVector ElementMatrixMap::mult(const RVector & a, const RVector & b,
                               const RVector & m, const RVector & n) const{
    RVector ret(rows_);

/// refactor me with ElementMatrix.mult
    for (Index r = 0; r < row_.size(); r ++ ){
        double s = 0.0;
        const SmallMatrix & mat = mat_[r];
        const IndexArray & idx = _ids[r];
        for (Index i = 0; i < mat.rows(); i ++) {
            double t = 0;
            for (Index j = 0; j < mat.cols(); j ++) {
                t += mat(i,j) * (a[idx[j]]-b[idx[j]]);
            }
            s += t * (m[idx[i]] - n[idx[i]]);
        }

        ret[row_[r]] += s;
    }

    return ret;
}

RVector ElementMatrixMap::mult(const RVector & a, const RVector & b) const{
    RVector ret(rows_);

    for (Index r = 0; r < row_.size(); r ++ ){
        double s = 0.0;
        const SmallMatrix & mat = mat_[r];
        const IndexArray & idx = _ids[r];
        for (Index i = 0; i < mat.rows(); i ++) {
            double t = 0;
            for (Index j = 0; j < mat.cols(); j ++) {
                t += mat(i,j) * (a[idx[j]]);
            }
            s += t * (b[idx[i]]);
        }

        ret[row_[r]] += s;
    }

    return ret;
}


void createUMap(const Mesh & mesh, Index order, ElementMatrixMap & ret,
                Index nCoeff, Index dofOffset){
    // don't use cache here // check!
    ret.resize(mesh.cellCount());

    for (auto &cell: mesh.cells()){
        ret.pMat(cell->id())->pot(*cell, order, true,
                                   nCoeff, mesh.nodeCount(), dofOffset);
    }
}

ElementMatrixMap createUMap(const Mesh & mesh, Index order,
                            Index nCoeff, Index dofOffset){
    ElementMatrixMap ret;
    createUMap(mesh, order, ret, nCoeff, dofOffset);
    return ret;
}

void createdUMap(const Mesh & mesh, Index order,
                 ElementMatrixMap & ret,
                 bool elastic, bool div, bool kelvin,
                 Index nCoeff, Index dofOffset){
        // don't use cache here // check!
    
    ret.resize(mesh.cellCount());

    for (auto &cell: mesh.cells()){
        // grad(const MeshEntity & ent, Index order,
        //      bool elastic, bool sum, bool div,
        //      Index nCoeff, Index dof, Index dofOffset,
        //      bool kelvin)
        ret.pMat(cell->id())->grad(*cell, order,
                                   elastic, false, div,
                                    nCoeff, mesh.nodeCount(), dofOffset,
                                    kelvin);
    }
}

ElementMatrixMap createdUMap(const Mesh & mesh, Index order,
                             bool elastic, bool div, bool kelvin,
                             Index nCoeff, Index dofOffset){
    ElementMatrixMap ret;
    createdUMap(mesh, order, ret, 
               elastic, div, kelvin, 
               nCoeff, dofOffset);
    return ret;
}

} // namespace GIMLI
