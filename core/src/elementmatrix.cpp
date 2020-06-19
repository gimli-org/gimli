/******************************************************************************
 *   Copyright (C) 2006-2020 by the GIMLi development team                    *
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

#include "elementmatrix.h"
#include "shape.h"
#include "meshentities.h"
#include "node.h"
#include "pos.h"

#include "integration.h"

namespace GIMLI{

template < >
std::ostream & operator << (std::ostream & str,
                            const ElementMatrix< double > & e){
    for (Index i = 0; i < e.colIDs().size(); i ++) str << e.colIDs()[i] << " " ;

    str << std::endl;
    for (Index i = 0; i < e.size(); i ++){
        str << e.rowIDs()[i] << "\t: ";
        for (Index j = 0; j < e.colIDs().size(); j ++){
            str << e.getVal(i, j) << " ";
        }
        str << std::endl;
    }
    return str;
}


template < > DLLEXPORT ElementMatrix < double > &
ElementMatrix < double >::u(const MeshEntity & ent,
                            const RVector & w,
                            const R3Vector & x,
                            bool verbose){

    uint nVerts = ent.nodeCount();
    std::map< uint, RVector >::const_iterator it = uCache_.find(ent.rtti());

    if (it == uCache_.end()) {
        uint nRules = w.size();

        RVector u(nVerts);
        RMatrix N(nVerts, nRules);

        RVector tmp;
        for (uint i = 0; i < nRules; i ++){
            tmp = ent.N(x[i]);
            N.setCol(i, tmp);
        }
        for (uint i = 0; i < nVerts; i ++){
            u[i] = sum(w * N[i]);
        }
        uCache_[ent.rtti()] = u;
        it = uCache_.find(ent.rtti());
    }

    double A = ent.shape().domainSize();
    // double J = det(ent.shape().createJacobian());
    // __MS(A << " " << J)
    for (Index i = 0; i < nVerts; i ++){
        // __MS(i << " " << A << " " << it->second[i])
        mat_[0][i] = A * it->second[i];
        if (this->_nDof > 0){
            if (ent.dim() == 2){
                mat_[nVerts].setVal(mat_[0][i], nVerts + i);
            }
            if (ent.dim() == 3){
                mat_[2*nVerts].setVal(mat_[0][i], 2*nVerts + i);
            }
        }
    }

    if (verbose) std::cout << "int u " << *this << std::endl;
    return *this;
}

template < >
DLLEXPORT ElementMatrix < double > & ElementMatrix < double >::u2(const MeshEntity & ent,
                                                        const RVector & w,
                                                        const R3Vector & x,
                                                        bool verbose){

    uint nVerts = ent.nodeCount();
    std::map< uint, RMatrix>::const_iterator it = u2Cache_.find(ent.rtti());

    if (it == u2Cache_.end()) {
        uint nRules = w.size();

        RMatrix u2(nVerts, nVerts);
        RMatrix N(nVerts, nRules);

        RVector tmp;
        for (uint i = 0; i < nRules; i ++){
            tmp = ent.N(x[i]);
            N.setCol(i, tmp);
        }
        for (uint i = 0; i < nVerts; i ++){
            for (uint j = i; j < nVerts; j ++){
                u2[i][j] = sum(w * N[j] * N[i]);
                u2[j][i] = u2[i][j];
            }
        }
        u2Cache_[ent.rtti()] = u2;
        it = u2Cache_.find(ent.rtti());
    }

    double A = ent.shape().domainSize();
//** very slow yet, pimp this with expressions (matrix::operator = (matExpression &))
//    mat_ = it->second * A;
    for (uint i = 0; i < nVerts; i ++){
        for (uint j = 0; j < nVerts; j ++){
            mat_[i][j] = A * it->second[i][j];
        }
    }

    if (verbose) std::cout << "int u2 " << *this << std::endl;

    return *this;
}

template <> DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::dudi(
                                                        const MeshEntity & ent,
                                                        const RVector & w,
                                                        const R3Vector & x,
                                                        Index dim,
                                                        bool verbose){

    this->fillIds(ent);
    Index nVerts = ent.nodeCount();
    Index nRules = w.size();

    if (dNdr_.rows() != nVerts){
        dNdr_.resize(nVerts, nRules);
        for (Index i = 0; i < nRules; i ++){
            dNdr_.setCol(i, ent.dNdL(x[i], 0));
            if (ent.dim() > 1){
                dNds_.resize(nVerts, nRules);
                dNds_.setCol(i, ent.dNdL(x[i], 1));
            }
            if (ent.dim() > 2){
                dNdt_.resize(nVerts, nRules);
                dNdt_.setCol(i, ent.dNdL(x[i], 2));
            }
        }
        dNdx_.resize(nVerts, nRules);
    }

    double drdi = ent.shape().drstdxyz(0, dim);
    double dsdi = ent.shape().drstdxyz(1, dim);
    double dtdi = ent.shape().drstdxyz(2, dim);

    double A = ent.shape().domainSize();

    for (Index i = 0; i < nVerts; i ++){
        switch (ent.dim()){
            case 1: dNdx_[i].assign(drdi * dNdr_[i]); break;
            case 2: dNdx_[i].assign(drdi * dNdr_[i] + dsdi * dNds_[i]); break;
            case 3: dNdx_[i].assign(drdi * dNdr_[i] + dsdi * dNds_[i] + dtdi * dNdt_[i]); break;
        }

        mat_[i][i] = sum(w * dNdx_[i]);
    }
    if (verbose) std::cout << "int dudx " << *this << std::endl;
    return *this;
}


template < >
DLLEXPORT ElementMatrix < double > & ElementMatrix < double >::ux2(const MeshEntity & ent,
                                                         const RVector & w,
                                                         const R3Vector & x,
                                                         bool verbose){

    uint nVerts = ent.nodeCount();
    uint nRules = w.size();

    if (dNdr_.rows() != nVerts){
        dNdr_.resize(nVerts, nRules);
        for (uint i = 0; i < nRules; i ++){
            dNdr_.setCol(i, ent.dNdL(x[i], 0));
        }
        dNdx_.resize(nVerts, nRules);
    }

    double drdx = ent.shape().drstdxyz(0, 0);
//     double drdx = ent.shape().invJacobian()[0];

    double A = ent.shape().domainSize();

    for (uint i = 0; i < nVerts; i ++){
        dNdx_[i].assign(drdx * dNdr_[i]);
    }
    for (uint i = 0; i < nVerts; i ++){
        for (uint j = i; j < nVerts; j ++){
            mat_[i][j] = A * sum(w * (dNdx_[i] * dNdx_[j]));
            mat_[j][i] = mat_[i][j];
        }
    }
    if (verbose) std::cout << "int ux2uy2 " << *this << std::endl;
    return *this;
}

template < > ElementMatrix < double > &
DLLEXPORT ElementMatrix < double >::ux2uy2(const MeshEntity & ent,
                                 const RVector & w,
                                 const R3Vector & x,
                                 bool verbose){

//     __MS(w)

    Index nVerts = ent.nodeCount();
    Index nRules = w.size();

    if (dNdr_.rows() != nVerts){
        dNdr_.resize(nVerts, nRules);
        dNds_.resize(nVerts, nRules);

        for (uint i = 0; i < nRules; i ++){
            dNdr_.setCol(i, ent.dNdL(x[i], 0));
            dNds_.setCol(i, ent.dNdL(x[i], 1));

//             __MS(i << " " << dNdr_[i])
        }

        dNdx_.resize(nVerts, nRules);
        dNdy_.resize(nVerts, nRules);
    }

//     double drdx = ent.shape().invJacobian()[0]; //ent.shape().drstdxyz(0,0)
//     double drdy = ent.shape().invJacobian()[1];
//     double dsdx = ent.shape().invJacobian()[2];
//     double dsdy = ent.shape().invJacobian()[3];

    double drdx = ent.shape().drstdxyz(0, 0);
    double drdy = ent.shape().drstdxyz(0, 1);
    double dsdx = ent.shape().drstdxyz(1, 0);
    double dsdy = ent.shape().drstdxyz(1, 1);

    double A = ent.shape().domainSize();
    for (Index i = 0; i < nVerts; i ++){
        dNdx_[i].assign(drdx * dNdr_[i] + dsdx * dNds_[i]);
        dNdy_[i].assign(drdy * dNdr_[i] + dsdy * dNds_[i]);
    }

    for (Index i = 0; i < nVerts; i ++){
        // dNidx_.assign(drdx * dNdr_[i] + dsdx * dNds_[i]);
        // dNidy_.assign(drdy * dNdr_[i] + dsdy * dNds_[i]);

        for (Index j = i; j < nVerts; j ++){
            mat_[i][j] = A * sum(w * (dNdx_[i] * dNdx_[j] + dNdy_[i] * dNdy_[j]));

//             mat_[i][j] = A * sum(w * (dNidx_ * (drdx * dNdr_[j] + dsdx * dNds_[j]) +
//                                       dNidy_ * (drdy * dNdr_[j] + dsdy * dNds_[j])));
//             mat_[i][j] = A * sum(w * ((drdx * dNdr_[i] + dsdx * dNds_[i]) * (drdx * dNdr_[j] + dsdx * dNds_[j]) +
//                                       (drdy * dNdr_[i] + dsdy * dNds_[i]) * (drdy * dNdr_[j] + dsdy * dNds_[j])));
            mat_[j][i] = mat_[i][j];
        }
    }

    if (verbose) std::cout << "int ux2uy2 " << *this << std::endl;
    return *this;
}

template < > DLLEXPORT ElementMatrix < double > &
ElementMatrix < double >::ux2uy2uz2(const MeshEntity & ent,
                                    const RVector & w,
                                    const R3Vector & x,
                                    bool verbose){

    Index nVerts = ent.nodeCount();
    Index nRules = w.size();

    if (dNdr_.rows() != nVerts){
        dNdr_.resize(nVerts, nRules);
        dNds_.resize(nVerts, nRules);
        dNdt_.resize(nVerts, nRules);

        for (Index i = 0; i < nRules; i ++){
            dNdr_.setCol(i, ent.dNdL(x[i], 0));
            dNds_.setCol(i, ent.dNdL(x[i], 1));
            dNdt_.setCol(i, ent.dNdL(x[i], 2));
        }

        dNdx_.resize(nVerts, nRules);
        dNdy_.resize(nVerts, nRules);
        dNdz_.resize(nVerts, nRules);
    }

    double drdx = ent.shape().drstdxyz(0, 0);
    double drdy = ent.shape().drstdxyz(0, 1);
    double drdz = ent.shape().drstdxyz(0, 2);

    double dsdx = ent.shape().drstdxyz(1, 0);
    double dsdy = ent.shape().drstdxyz(1, 1);
    double dsdz = ent.shape().drstdxyz(1, 2);

    double dtdx = ent.shape().drstdxyz(2, 0);
    double dtdy = ent.shape().drstdxyz(2, 1);
    double dtdz = ent.shape().drstdxyz(2, 2);

    // double drdx = ent.shape().invJacobian()[0];
    // double drdy = ent.shape().invJacobian()[1];
    // double drdz = ent.shape().invJacobian()[2];
    //
    // double dsdx = ent.shape().invJacobian()[3];
    // double dsdy = ent.shape().invJacobian()[4];
    // double dsdz = ent.shape().invJacobian()[5];
    //
    // double dtdx = ent.shape().invJacobian()[6];
    // double dtdy = ent.shape().invJacobian()[7];
    // double dtdz = ent.shape().invJacobian()[8];

    double A = ent.shape().domainSize();
    for (Index i = 0; i < nVerts; i ++){
        dNdx_[i].assign(drdx * dNdr_[i] + dsdx * dNds_[i] + dtdx * dNdt_[i]);
        dNdy_[i].assign(drdy * dNdr_[i] + dsdy * dNds_[i] + dtdy * dNdt_[i]);
        dNdz_[i].assign(drdz * dNdr_[i] + dsdz * dNds_[i] + dtdz * dNdt_[i]);
    }

    for (Index i = 0; i < nVerts; i ++){
        for (Index j = i; j < nVerts; j ++){
            mat_[i][j] = A * sum(w * (dNdx_[i] * dNdx_[j] + dNdy_[i] * dNdy_[j] + dNdz_[i] * dNdz_[j]));

//             mat_[i][j] = A * sum(w * ((drdx * dNdr_[i] + dsdx * dNds_[i] + dtdx * dNdt_[i]) *
//                                       (drdx * dNdr_[j] + dsdx * dNds_[j] + dtdx * dNdt_[j]) +
//                                       (drdy * dNdr_[i] + dsdy * dNds_[i] + dtdy * dNdt_[i]) *
//                                       (drdy * dNdr_[j] + dsdy * dNds_[j] + dtdy * dNdt_[j]) +
//                                       (drdz * dNdr_[i] + dsdz * dNds_[i] + dtdz * dNdt_[i]) *
//                                       (drdz * dNdr_[j] + dsdz * dNds_[j] + dtdz * dNdt_[j])));

            mat_[j][i] = mat_[i][j];
        }
    }
    if (verbose) std::cout << "int ux2uy2uz2 " << *this << std::endl;
    return *this;
}

template < class ValueType >
void ElementMatrix < ValueType >::getWeightsAndPoints(const MeshEntity & ent,
const RVector * &w, const R3Vector * &x, int order){
    switch (ent.rtti()) {
        case MESH_EDGE_CELL_RTTI:
        case MESH_EDGE3_CELL_RTTI: {
            w = &IntegrationRules::instance().edgWeights(2);
            x = &IntegrationRules::instance().edgAbscissa(2);
        } break;
        case MESH_TRIANGLE_RTTI: {
            w = &IntegrationRules::instance().triWeights(1);
            x = &IntegrationRules::instance().triAbscissa(1);
        } break;
        case MESH_TRIANGLE6_RTTI: {
            w = &IntegrationRules::instance().triWeights(2);
            x = &IntegrationRules::instance().triAbscissa(2);
        } break;
        case MESH_QUADRANGLE_RTTI: {
            w = &IntegrationRules::instance().quaWeights(2);
            x = &IntegrationRules::instance().quaAbscissa(2);
        } break;
        case MESH_QUADRANGLE8_RTTI: {
            w = &IntegrationRules::instance().quaWeights(3);
            x = &IntegrationRules::instance().quaAbscissa(3);
        } break;
        case MESH_TETRAHEDRON_RTTI: {
            w = &IntegrationRules::instance().tetWeights(1);
            x = &IntegrationRules::instance().tetAbscissa(1);
        } break;
        case MESH_TETRAHEDRON10_RTTI: {
            w = &IntegrationRules::instance().tetWeights(2);
            x = &IntegrationRules::instance().tetAbscissa(2);
        } break;
        case MESH_HEXAHEDRON_RTTI: {
            w = &IntegrationRules::instance().hexWeights(2);
            x = &IntegrationRules::instance().hexAbscissa(2);
        } break;
        case MESH_HEXAHEDRON20_RTTI: {
            w = &IntegrationRules::instance().hexWeights(4);
            x = &IntegrationRules::instance().hexAbscissa(4);
        } break;
        case MESH_TRIPRISM_RTTI: {
            w = &IntegrationRules::instance().priWeights(2);
            x = &IntegrationRules::instance().priAbscissa(2);
        } break;
        case MESH_TRIPRISM15_RTTI: {
            w = &IntegrationRules::instance().priWeights(4);
            x = &IntegrationRules::instance().priAbscissa(4);
        } break;
        default:
            std::cerr << ent.rtti() << std::endl;
            THROW_TO_IMPL
            break;
    }
}
template void ElementMatrix < double >::getWeightsAndPoints(
                                        const MeshEntity & ent,
                                        const RVector * &w,
                                        const R3Vector * &x, int order);

template < > DLLEXPORT
void ElementMatrix < double >::fillGradientBase(
                                    const MeshEntity & ent,
                                    const RVector & w,
                                    const R3Vector & x,
                                    Index nC,
                                    bool voigtNotation){

    Index nRules = x.size();
    Index nDof = this->mat_.cols();
    Index nVerts = ent.nodeCount();
    if (_B.size() != nRules){
        _B.resize(nRules);
        for (Index i = 0; i < nRules; i ++ ){
            _B[i].resize(nC, nDof);
        }
    }
    if (dNdr_.rows() != nRules){
        if (ent.dim() > 0) dNdr_.resize(nRules, nVerts);
        if (ent.dim() > 1) dNds_.resize(nRules, nVerts);
        if (ent.dim() > 2) dNdt_.resize(nRules, nVerts);

        for (Index i = 0; i < nRules; i ++){
            if (ent.dim() > 0) dNdr_[i] = ent.dNdL(x[i], 0);
            if (ent.dim() > 1) dNds_[i] = ent.dNdL(x[i], 1);
            if (ent.dim() > 2) dNdt_[i] = ent.dNdL(x[i], 2);
        }

        if (ent.dim() > 0) dNdx_.resize(nRules, nVerts);
        if (ent.dim() > 1) dNdy_.resize(nRules, nVerts);
        if (ent.dim() > 2) dNdz_.resize(nRules, nVerts);
    }
    double drdx = ent.shape().drstdxyz(0, 0);
    double drdy = ent.shape().drstdxyz(0, 1);
    double drdz = ent.shape().drstdxyz(0, 2);
    double dsdx = ent.shape().drstdxyz(1, 0);
    double dsdy = ent.shape().drstdxyz(1, 1);
    double dsdz = ent.shape().drstdxyz(1, 2);
    double dtdx = ent.shape().drstdxyz(2, 0);
    double dtdy = ent.shape().drstdxyz(2, 1);
    double dtdz = ent.shape().drstdxyz(2, 2);

    for (Index i = 0; i < nRules; i ++){
        if (ent.dim() == 1){
            dNdx_[i].assign(dNdr_[i] * drdx);
        } else if (ent.dim() == 2){
            dNdx_[i].assign(dNdr_[i] * drdx + dNds_[i] * dsdx);
            dNdy_[i].assign(dNdr_[i] * drdy + dNds_[i] * dsdy);
        } else if (ent.dim() == 3){
            dNdx_[i].assign(dNdr_[i] * drdx + dNds_[i] * dsdx + dNdt_[i] * dtdx);
            dNdy_[i].assign(dNdr_[i] * drdy + dNds_[i] * dsdy + dNdt_[i] * dtdy);
            dNdz_[i].assign(dNdr_[i] * drdz + dNds_[i] * dsdz + dNdt_[i] * dtdz);
        }
    }

    double a = std::sqrt(2.);
    if (voigtNotation){
        a = 1.0;
    }

    for (Index i = 0; i < nRules; i ++){
        // __MS(i << " " << x[i])
        //dNdx_[i] = (dNdr_[i] * drdx);
        //dNdx_.setCol(i, (dNdr_[i] * drdx));

        if (this->_nDof == 0){
            // scalar field
            if (ent.dim() > 0){
                _B[i][0].setVal(dNdx_[i], 0, nVerts);
            }
            if (ent.dim() > 1){
                _B[i][1].setVal(dNdy_[i], 0, nVerts);
            }
            if (ent.dim() > 2){
                _B[i][2].setVal(dNdz_[i], 0, nVerts);
            }
        } else {
            // vector field
            if (ent.dim() == 1){
                _B[i][0].setVal(dNdx_[i], 0 * nVerts, 1 * nVerts);
            } else if (ent.dim() == 2){
                _B[i][0].setVal(dNdx_[i], 0 * nVerts, 1 * nVerts);
                _B[i][1].setVal(dNdy_[i], 1 * nVerts, 2 * nVerts);

                if (nC > ent.dim()){ // elastic material
                    _B[i][2].setVal(dNdy_[i] * a, 0 * nVerts, 1 * nVerts);
                    _B[i][2].setVal(dNdx_[i] * a, 1 * nVerts, 2 * nVerts);
                }
            } else if (ent.dim() == 3){
                _B[i][0].setVal(dNdx_[i], 0 * nVerts, 1 * nVerts);
                _B[i][1].setVal(dNdy_[i], 1 * nVerts, 2 * nVerts);
                _B[i][2].setVal(dNdz_[i], 2 * nVerts, 3 * nVerts);

                if (nC > ent.dim()){ // elastic material
                    _B[i][3].setVal(dNdy_[i] * a, 0 * nVerts, 1 * nVerts);
                    _B[i][3].setVal(dNdx_[i] * a, 1 * nVerts, 2 * nVerts);

                    _B[i][4].setVal(dNdz_[i] * a, 1 * nVerts, 2 * nVerts);
                    _B[i][4].setVal(dNdy_[i] * a, 2 * nVerts, 3 * nVerts);

                    _B[i][5].setVal(dNdz_[i] * a, 0 * nVerts, 1 * nVerts);
                    _B[i][5].setVal(dNdx_[i] * a, 2 * nVerts, 3 * nVerts);
                }
            }
        }
    }
}

template < > DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::gradU(const MeshEntity & ent,
                                                           const RVector & w,
                                                           const R3Vector & x,
                                                           Index nC,
                                                           bool voigtNotation){
    this->fillIds(ent, nC); // also cleans
    this->fillGradientBase(ent, w, x, nC, voigtNotation);

    for (Index i = 0; i < w.size(); i ++ ){
        for (Index j = 0; j < nC; j ++) {
            // __MS(*this)
            // __MS(_B[i])
            // __MS(_B[i][j])
            mat_[j * ent.nodeCount()] += _B[i][j] * w[i] * ent.size();
        }
        // check performance if this works
        // iterator over weights in fill Gradient
        // this *= ent.size();
    }
    //mat_ *= ent.size();

    return * this;
}

template < > DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::gradU(const Cell & cell,
                                 Index nC,
                                 bool voigtNotation){
    const RVector * w = 0;
    const R3Vector * x = 0;
    this->getWeightsAndPoints(cell, w, x, 1);
    return this->gradU(cell, *w, *x, nC, voigtNotation);
}

template < > DLLEXPORT
RVector ElementMatrix < double >::stress(const MeshEntity & ent,
                                         const RMatrix & C,
                                         const RVector & u, bool voigtNotation){
    const RVector * w = 0;
    const R3Vector * x = 0;

    this->getWeightsAndPoints(ent, w, x, 1);
    this->fillIds(ent, C.size()); // also cleans
    this->fillGradientBase(ent, *w, *x,
                           max(C.size(), ent.dim()),
                           voigtNotation);

    RVector ret(C.rows(), 0.0);
    for (Index i = 0; i < w->size(); i ++ ){
        // C * B * u
        ret += C * (_B[i] * u[_ids]) * (*w)[i];
    }

    // performance optimization using _grad from fill GradientBase
    // _grad = sum(_B*w) // sum in GradientBase
    //ret =  C * (_grad * u[_ids])

    return ret;
}

template < > DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::gradU2(const MeshEntity & ent,
                                                            const Matrix< double > & C,
                                                            const RVector & w,
                                                            const R3Vector & x,
                                                            bool voigtNotation){
    this->fillIds(ent, C.size()); // also cleans
    this->fillGradientBase(ent, w, x,
                           max(C.size(), ent.dim()),
                           voigtNotation);
    if (C.size() == 1){
        for (Index i = 0; i < w.size(); i ++ ){
            // B.T * B
            matTransMult(_B[i], _B[i], mat_, w[i] * ent.size() * C[0][0]);
        }
    } else {
        for (Index i = 0; i < w.size(); i ++ ){
            // B.T * C * B
            matMultABA(_B[i], C, mat_, _abaTmp, w[i] * ent.size());
        }
        // check performance if this works
        // iterator over weights in fill Gradient
        // matMultABA(_B, C, mat_, _abaTmp, 1.0);
        // this *= ent.size();
    }
    // M += w[i] * B.T @ C @ B
    //mat_ *= ent.size();

    return * this;
}


template < > DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::gradU2(const Cell & cell,
                                                    const Matrix< double > & C,
                                                    bool voigtNotation){

    const RVector * w = 0;
    const R3Vector * x = 0;
    this->getWeightsAndPoints(cell, w, x, 1);
    return this->gradU2(cell, C, *w, *x, voigtNotation);
}


template < > DLLEXPORT ElementMatrix < double > &
ElementMatrix < double >::u(const MeshEntity & ent){
    this->fillIds(ent);

    switch(ent.rtti()){
        case MESH_BOUNDARY_NODE_RTTI:
            mat_[0][0] = 1.0;
        break;
        case MESH_EDGE_CELL_RTTI:
        case MESH_EDGE_RTTI:
        case MESH_EDGE3_CELL_RTTI:
        case MESH_EDGE3_RTTI:
            return u(ent, IntegrationRules::instance().edgWeights(2),
                     IntegrationRules::instance().edgAbscissa(2), false); //ch
        case MESH_TRIANGLE_RTTI:
        case MESH_TRIANGLEFACE_RTTI:
            return u(ent, IntegrationRules::instance().triWeights(2),
                     IntegrationRules::instance().triAbscissa(2), false); //ch
        case MESH_TRIANGLE6_RTTI:
        case MESH_TRIANGLEFACE6_RTTI:
            return u(ent, IntegrationRules::instance().triWeights(2),
                     IntegrationRules::instance().triAbscissa(2), false); //ch
        case MESH_QUADRANGLE_RTTI:
        case MESH_QUADRANGLE8_RTTI:
            return u(ent, IntegrationRules::instance().quaWeights(2),
                     IntegrationRules::instance().quaAbscissa(2), false); //ch
        case MESH_QUADRANGLEFACE_RTTI:
        case MESH_QUADRANGLEFACE8_RTTI:
            return u(ent, IntegrationRules::instance().quaWeights(2),
                     IntegrationRules::instance().quaAbscissa(2), false); //ch
        case MESH_TETRAHEDRON_RTTI:
        case MESH_TETRAHEDRON10_RTTI:
            return u(ent, IntegrationRules::instance().tetWeights(2),
                     IntegrationRules::instance().tetAbscissa(2), false);
        case MESH_HEXAHEDRON_RTTI:
        case MESH_HEXAHEDRON20_RTTI:
            return u(ent, IntegrationRules::instance().hexWeights(2),
                     IntegrationRules::instance().hexAbscissa(2), false);
        case MESH_TRIPRISM_RTTI:
        case MESH_TRIPRISM15_RTTI:
            return u(ent, IntegrationRules::instance().priWeights(2),
                     IntegrationRules::instance().priAbscissa(2), false);

        default: std::cerr << WHERE_AM_I << " celltype not specified " << ent.rtti() << std::endl;
    }
    return *this;
}

template < > DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::u2(const MeshEntity & ent){
    this->fillIds(ent);

    switch(ent.rtti()){
    case MESH_BOUNDARY_NODE_RTTI:
        mat_[0][0] = 1.0;
    break;
    case MESH_EDGE_CELL_RTTI:
    case MESH_EDGE_RTTI:
    //{
//         double J = ent.shape().jacobianDeterminant();
//         if (J < 0) {
//             std::cerr << WHERE_AM_I << "JacobianDeterminant < 0 (" << J << ")" << std::endl;
//             std::cerr << ent.shape() << std::endl;
//             J = std::fabs(J);
//         }
//         mat_[0][0] = J / 3.0;
//         mat_[1][0] = J / 6.0;
//         mat_[0][1] = J / 6.0;
//         mat_[1][1] = J / 3.0;
//         std::cout << "2 " << *this << std::endl;
/*}
    break;*/
        u2(ent, IntegrationRules::instance().edgWeights(2), IntegrationRules::instance().edgAbscissa(2), false);
        break;
    case MESH_EDGE3_CELL_RTTI:
    case MESH_EDGE3_RTTI:
    //{
//         double J = ent.shape().jacobianDeterminant();
//         if (J < 0) {
//             std::cerr << WHERE_AM_I << "JacobianDeterminant < 0 (" << J << ")" << std::endl;
//             std::cerr << ent.shape() << std::endl;
//             J = std::fabs(J);
//         }
//         mat_[0][0] =   J / 30.0 * 4.0;
//         mat_[1][0] = - J / 30.0;
//         mat_[2][0] =   J / 15.0;
//         mat_[0][1] = - J / 30.0;
//         mat_[1][1] =   J / 30.0 * 4.0;
//         mat_[2][1] =   J / 15.0;
//         mat_[0][2] =   J / 15.0;
//         mat_[1][2] =   J / 15.0;
//         mat_[2][2] =   J / 30.0 * 16.0;
//          std::cout << "2 " << *this << std::endl;
    //} break;
        u2(ent, IntegrationRules::instance().edgWeights(3), IntegrationRules::instance().edgAbscissa(3), false);
        break;
    case MESH_TRIANGLE_RTTI:
    case MESH_TRIANGLEFACE_RTTI:
    //{
//         double J = ent.shape().jacobianDeterminant();
//         if (J < 0) {
//             std::cerr << WHERE_AM_I << "JacobianDeterminant < 0 (" << J << ")" << std::endl;
//             std::cerr << ent.shape() << std::endl;
//             J = std::fabs(J);
//         }
//         double Jl = J / 24.0;
//         mat_[0][0] = Jl * 2.0;
//         mat_[1][0] = Jl;
//         mat_[2][0] = Jl;
//         mat_[0][1] = Jl;
//         mat_[1][1] = Jl * 2.0;
//         mat_[2][1] = Jl;
//         mat_[0][2] = Jl;
//         mat_[1][2] = Jl;
//         mat_[2][2] = Jl * 2.0;
//         std::cout << "2 " << *this << std::endl;
//}break;
        u2(ent, IntegrationRules::instance().triWeights(2), IntegrationRules::instance().triAbscissa(2), false);
        break;
    case MESH_QUADRANGLE_RTTI:
    case MESH_QUADRANGLEFACE_RTTI:
        u2(ent, IntegrationRules::instance().quaWeights(2), IntegrationRules::instance().quaAbscissa(2), false);
        break;
    case MESH_QUADRANGLE8_RTTI:
    case MESH_QUADRANGLEFACE8_RTTI:
        u2(ent, IntegrationRules::instance().quaWeights(3), IntegrationRules::instance().quaAbscissa(3), false);
        break;
    case MESH_TRIANGLE6_RTTI:
    case MESH_TRIANGLEFACE6_RTTI:
    //{
//         double J = ent.shape().jacobianDeterminant();
//         if (J < 0) {
//             std::cerr << WHERE_AM_I << "JacobianDeterminant < 0 (" << J << ")" << std::endl;
//             std::cerr << ent.shape() << std::endl;
//             J = std::fabs(J);
//         }
//         for (uint i = 0; i < nVerts; i++){
//             for (uint j = 0; j < nVerts; j++){
//                 mat_[i][j] = J * (*Tri6_u2)[i][j];//compound[i][j];
//             }
//         }
//         std::cout << "2 " << *this << std::endl;
    //} break;
        return u2(ent, IntegrationRules::instance().triWeights(4), IntegrationRules::instance().triAbscissa(4), false); //ch
    case MESH_TETRAHEDRON_RTTI:
        return u2(ent, IntegrationRules::instance().tetWeights(2), IntegrationRules::instance().tetAbscissa(2), false);
    case MESH_TETRAHEDRON10_RTTI:
        return u2(ent, IntegrationRules::instance().tetWeights(4), IntegrationRules::instance().tetAbscissa(4), false);
    case MESH_HEXAHEDRON_RTTI:
        return u2(ent, IntegrationRules::instance().hexWeights(2), IntegrationRules::instance().hexAbscissa(2), false);
    case MESH_HEXAHEDRON20_RTTI:
        return u2(ent, IntegrationRules::instance().hexWeights(4), IntegrationRules::instance().hexAbscissa(4), false);
    case MESH_TRIPRISM_RTTI:
        return u2(ent, IntegrationRules::instance().priWeights(2), IntegrationRules::instance().priAbscissa(2), false);
    case MESH_TRIPRISM15_RTTI:
        return u2(ent, IntegrationRules::instance().priWeights(4), IntegrationRules::instance().priAbscissa(4), false);
    default:
      std::cerr << ent.rtti() << std::endl;
      THROW_TO_IMPL
    }

    return *this;
}

template < > DLLEXPORT ElementMatrix < double > &
ElementMatrix < double >::ux2uy2uz2(const Cell & cell, bool useCache){

    fillIds(cell);

    if (cell.uxCache().rows() > 0 && useCache){
        mat_ = cell.uxCache();
        return *this;
    }

//     double J = cell.jacobianDeterminant();
//     if (J <= 0) std::cerr << WHERE_AM_I << " JacobianDeterminant < 0 (" << J << ") " << cell << std::endl;
//      std::cout << J << std::endl;

    switch (cell.rtti()) {
    case MESH_EDGE_CELL_RTTI:
    case MESH_EDGE3_CELL_RTTI:
        ux2(cell, IntegrationRules::instance().edgWeights(2),
            IntegrationRules::instance().edgAbscissa(2), false); break;
    case MESH_TRIANGLE_RTTI: {
        double J = cell.size()*2.;
    ////////////////////////////////////////////////////////////////////
//         double dN1dx = cell.shape().deriveCoordinates(0, 0);
//         double dN2dx = cell.shape().deriveCoordinates(1, 0);
//         double dN3dx = cell.shape().deriveCoordinates(2, 0);
//         double dN1dy = cell.shape().deriveCoordinates(0, 1);
//         double dN2dy = cell.shape().deriveCoordinates(1, 1);
//         double dN3dy = cell.shape().deriveCoordinates(2, 1);
//         mat_[0][0] = J / 2.0 * (dN1dx * dN1dx + dN1dy * dN1dy);
//         mat_[1][0] = J / 2.0 * (dN2dx * dN1dx + dN2dy * dN1dy);
//         mat_[2][0] = J / 2.0 * (dN3dx * dN1dx + dN3dy * dN1dy);
//         mat_[1][1] = J / 2.0 * (dN2dx * dN2dx + dN2dy * dN2dy);
//         mat_[2][1] = J / 2.0 * (dN3dx * dN2dx + dN3dy * dN2dy);
//         mat_[2][2] = J / 2.0 * (dN3dx * dN3dx + dN3dy * dN3dy);
//         mat_[0][1] = mat_[1][0];
//         mat_[0][2] = mat_[2][0];
//         mat_[1][2] = mat_[2][1];
////////////////////////////////////////////////////////////////////
        // this is much faster than numerical integration

//         double x21 = cell.node(1).x()-cell.node(0).x();
//         double x31 = cell.node(2).x()-cell.node(0).x();
//         double y21 = cell.node(1).y()-cell.node(0).y();
//         double y31 = cell.node(2).y()-cell.node(0).y();
//
//         double a =   ((x31) * (x31) + (y31) * (y31)) / J;
//         double b = - ((x31) * (x21) + (y31) * (y21)) / J;
//         double c =   ((x21) * (x21) + (y21) * (y21)) / J;

        double x1 = cell.node(0).x();
        double x2 = cell.node(1).x();
        double x3 = cell.node(2).x();
        double y1 = cell.node(0).y();
        double y2 = cell.node(1).y();
        double y3 = cell.node(2).y();

        double a =   ((x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1)) / J;
        double b = - ((x3 - x1) * (x2 - x1) + (y3 - y1) * (y2 - y1)) / J;
        double c =   ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)) / J;

        mat_[0][0] = a *  0.5 + b        + c *  0.5 ;
        mat_[1][0] = a * -0.5 + b * -0.5            ;
        mat_[2][0] =            b * -0.5 + c * -0.5 ;
        mat_[1][1] = a *  0.5                       ;
        mat_[2][1] =            b *  0.5            ;
        mat_[2][2] =                       c *  0.5 ;

        mat_[0][1] = mat_[1][0];
        mat_[0][2] = mat_[2][0];
        mat_[1][2] = mat_[2][1];

        //std::cout << "2" << *this << std::endl;*/
        //ux2uy2(cell, IntegrationRules::instance().triWeights(1), IntegrationRules::instance().triAbscissa(1), false);
    } break;
    case MESH_TRIANGLE6_RTTI: {
///////////////////////////////////////////////////////////////////////////////////
//         double x1 = cell.node(0).x();
//         double x2 = cell.node(1).x();
//         double x3 = cell.node(2).x();
//         double y1 = cell.node(0).y();
//         double y2 = cell.node(1).y();
//         double y3 = cell.node(2).y();
//
//         double b1 = y2-y3;
//         double b2 = y3-y1;
//         double b3 = y1-y2;
//         double c1 = x3-x2;
//         double c2 = x1-x3;
//         double c3 = x2-x1;
//         double a = (c3*b2-c2*b3)/2.0;
//         double b1sqr=b1*b1;
//         double b2sqr=b2*b2;
//         double b3sqr=b3*b3;
//         double c1sqr=c1*c1;
//         double c2sqr=c2*c2;
//         double c3sqr=c3*c3;
//         double b1b2=b1*b2;
//         double b1b3=b1*b3;
//         double b2b3=b2*b3;
//         double c1c2=c1*c2;
//         double c1c3=c1*c3;
//         double c2c3=c2*c3;
//
//         mat_[0][0]=(b1sqr+c1sqr) / (4.0*a);
//         mat_[0][1]=(-b1b2-c1c2)/(12.0*a);
//         mat_[0][2]=(-b1b3-c1c3)/(12.0*a);
//         mat_[0][3]=(b1b2+c1c2)/(3.0*a);
//         mat_[0][4]=0.0;
//         mat_[0][5]=(b1b3+c1c3)/(3.0*a);
//         mat_[1][1]=(b2sqr+c2sqr)/(4.0*a);
//         mat_[1][2]=(-b2b3-c2c3)/(12.0*a);
//         mat_[1][3]=(b1b2+c1c2)/(3.0*a);
//         mat_[1][4]=(b2b3+c2c3)/(3.0*a);
//         mat_[1][5]=0.0;
//         mat_[2][2]=(b3sqr+c3sqr)/(4.0*a);
//         mat_[2][3]=0.0;
//         mat_[2][4]=(b2b3+c2c3)/(3.0*a);
//         mat_[2][5]=(b1b3+c1c3)/(3.0*a);
//         mat_[3][3]=2.0*((b1sqr+b1b2+b2sqr)+(c1sqr+c1c2+c2sqr))/(3.0*a);
//         mat_[3][4]=((b1b2+b2sqr+2.0*b1b3+b2b3)+(c1c2+c2sqr+2.0*c1c3+c2c3))/(3.0*a);
//         mat_[3][5]=((b1sqr+b1b3+b1b2+2.0*b2b3)+(c1sqr+c1c3+c1c2+2.0*c2c3))/(3.0*a);
//         mat_[4][4]=2.0*((b2sqr+b2b3+b3sqr)+(c2sqr+c2c3+c3sqr))/(3.0*a);
//         mat_[4][5]=((2.0*b1b2+b2b3+b1b3+b3sqr)+(2.0*c1c2+c2c3+c1c3+c3sqr))/(3.0*a);
//         mat_[5][5]=2.0*((b1sqr+b1b3+b3sqr)+(c1sqr+c1c3+c3sqr))/(3.0*a);
//
//         for (int i = 1, imax = 6; i < imax; i ++){
//              for (int j = 0, jmax = i; j < jmax; j ++){
//                 mat_[i][j]=mat_[j][i];
//              }
//         }

//working version
//         double x1 = cell.node(0).x();
//         double x2 = cell.node(1).x();
//         double x3 = cell.node(2).x();
//         double y1 = cell.node(0).y();
//         double y2 = cell.node(1).y();
//         double y3 = cell.node(2).y();
//
//         double a = ((x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1)) / J / 6.0;
//         double b = - ((x3 - x1) * (x2 - x1) + (y3 - y1) * (y2 - y1)) / J / 6.0;
//         double c = ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)) / J / 6.0;
//         for (int i = 0, imax = 6; i < imax; i ++){
//             for (int j = 0, jmax = 6; j < jmax; j ++){
//                 mat_[i][j] = Triangle6_S1[i][j] * a +
//                                  Triangle6_S2[i][j] * b +
//                                  Triangle6_S3[i][j] * c;
//             }
//         }
//         std::cout << "2" << *this << std::endl;
        ux2uy2(cell,
               IntegrationRules::instance().triWeights(2),
               IntegrationRules::instance().triAbscissa(2), false); //ch
    } break;
    case MESH_QUADRANGLE_RTTI:
        ux2uy2(cell,
               IntegrationRules::instance().quaWeights(2),
               IntegrationRules::instance().quaAbscissa(2), false); break;
    case MESH_QUADRANGLE8_RTTI:
        ux2uy2(cell,
               IntegrationRules::instance().quaWeights(3),
               IntegrationRules::instance().quaAbscissa(3), false); break;
    case MESH_TETRAHEDRON_RTTI:
    //{
//         double x_xi = cell.shape().partDerivationRealToUnity(0, 1);
//         double x_eta = cell.shape().partDerivationRealToUnity(0, 2);
//         double x_zeta = cell.shape().partDerivationRealToUnity(0, 3);
//
//         double y_xi = cell.shape().partDerivationRealToUnity(1, 1);
//         double y_eta = cell.shape().partDerivationRealToUnity(1, 2);
//         double y_zeta = cell.shape().partDerivationRealToUnity(1, 3);
//
//         double z_xi = cell.shape().partDerivationRealToUnity(2, 1);
//         double z_eta = cell.shape().partDerivationRealToUnity(2, 2);
//         double z_zeta = cell.shape().partDerivationRealToUnity(2, 3);
//
//         double xi_x =  1.0 / J * det(y_eta, z_eta, y_zeta, z_zeta);
//         double xi_y = -1.0 / J * det(x_eta, z_eta, x_zeta, z_zeta);
//         double xi_z =  1.0 / J * det(x_eta, y_eta, x_zeta, y_zeta);
//
//         double eta_x = -1.0 / J * det(y_xi, z_xi, y_zeta, z_zeta);
//         double eta_y =  1.0 / J * det(x_xi, z_xi, x_zeta, z_zeta);
//         double eta_z = -1.0 / J * det(x_xi, y_xi, x_zeta, y_zeta);
//
//         double zeta_x =  1.0 / J * det(y_xi, z_xi, y_eta, z_eta);
//         double zeta_y = -1.0 / J * det(x_xi, z_xi, x_eta, z_eta);
//         double zeta_z =  1.0 / J * det(x_xi, y_xi, x_eta, y_eta);
//
//         double a = J / 6.0 * (xi_x * xi_x + xi_y * xi_y + xi_z * xi_z);
//         double b = J / 6.0 * (eta_x * eta_x + eta_y * eta_y + eta_z * eta_z);
//         double c = J / 6.0 * (zeta_x * zeta_x + zeta_y * zeta_y + zeta_z * zeta_z);
//
//         double d = J / 6.0 * (xi_x * eta_x + xi_y * eta_y + xi_z * eta_z);
//         double e = J / 6.0 * (xi_x * zeta_x + xi_y * zeta_y + xi_z * zeta_z);
//         double f = J / 6.0 * (eta_x * zeta_x + eta_y * zeta_y + eta_z * zeta_z);
//
//     double u_xi2[4][4] = {
//       {   a,  -a, 0.0, 0.0 },
//       {  -a,   a, 0.0, 0.0 },
//       { 0.0, 0.0, 0.0, 0.0 },
//       { 0.0, 0.0, 0.0, 0.0 }
//     };
//     double u_eta2[4][4] = {
//       {   b, 0.0,  -b, 0.0 },
//       { 0.0, 0.0, 0.0, 0.0 },
//       {  -b, 0.0,   b, 0.0 },
//       { 0.0, 0.0, 0.0, 0.0 }
//     };
//     double u_zeta2[4][4] = {
//       {   c, 0.0, 0.0,  -c },
//       { 0.0, 0.0, 0.0, 0.0 },
//       { 0.0, 0.0, 0.0, 0.0 },
//       {  -c, 0.0, 0.0,   c }
//     };
//     double u_xi__u_eta[4][4] = {
//       { 2.0 * d,  -d,  -d, 0.0 },
//       {      -d, 0.0,   d, 0.0 },
//       {      -d,   d, 0.0, 0.0 },
//       {     0.0, 0.0, 0.0, 0.0 }
//     };
//     double u_xi__u_zeta[4][4] = {
//       { 2.0 * e,  -e, 0.0,  -e },
//       {      -e, 0.0, 0.0,   e },
//       {     0.0, 0.0, 0.0, 0.0 },
//       {      -e,   e, 0.0, 0.0 }
//     };
//     double u_eta__u_zeta[4][4] = {
//       { 2.0 * f, 0.0,  -f,  -f },
//       {     0.0, 0.0, 0.0, 0.0 },
//       {      -f, 0.0, 0.0,   f },
//       {      -f, 0.0,   f, 0.0 }
//     };
//
//     for (uint i = 0; i < dim; i++){
//       for (uint j = 0; j < dim; j++){
//  mat_[i][j] = u_xi2[i][j] + u_eta2[i][j] +  u_zeta2[i][j] +
//    u_xi__u_eta[i][j] + u_xi__u_zeta[i][j] + u_eta__u_zeta[i][j];
//       }
//     }
//     std::cout << "0 " << *this << std::endl;
//} break;
        ux2uy2uz2(cell,
                  IntegrationRules::instance().tetWeights(1),
                  IntegrationRules::instance().tetAbscissa(1), false); //ch
        break;
    case MESH_TETRAHEDRON10_RTTI:
//     {
// ////////////////////////////////////////////////////////////////////
//     double x_xi = cell.shape().partDerivationRealToUnity(0, 1);
//     double x_eta = cell.shape().partDerivationRealToUnity(0, 2);
//     double x_zeta = cell.shape().partDerivationRealToUnity(0, 3);
//
//     double y_xi = cell.shape().partDerivationRealToUnity(1, 1);
//     double y_eta = cell.shape().partDerivationRealToUnity(1, 2);
//     double y_zeta = cell.shape().partDerivationRealToUnity(1, 3);
//
//     double z_xi = cell.shape().partDerivationRealToUnity(2, 1);
//     double z_eta = cell.shape().partDerivationRealToUnity(2, 2);
//     double z_zeta = cell.shape().partDerivationRealToUnity(2, 3);
//
//     double xi_x =  1.0 / J * det(y_eta, z_eta, y_zeta, z_zeta);
//     double xi_y = -1.0 / J * det(x_eta, z_eta, x_zeta, z_zeta);
//     double xi_z =  1.0 / J * det(x_eta, y_eta, x_zeta, y_zeta);
//
//     double eta_x = -1.0 / J * det(y_xi, z_xi, y_zeta, z_zeta);
//     double eta_y =  1.0 / J * det(x_xi, z_xi, x_zeta, z_zeta);
//     double eta_z = -1.0 / J * det(x_xi, y_xi, x_zeta, y_zeta);
//
//     double zeta_x =  1.0 / J * det(y_xi, z_xi, y_eta, z_eta);
//     double zeta_y = -1.0 / J * det(x_xi, z_xi, x_eta, z_eta);
//     double zeta_z =  1.0 / J * det(x_xi, y_xi, x_eta, y_eta);
//
//     double a = J / 6.0 * (xi_x * xi_x + xi_y * xi_y + xi_z * xi_z);
//     double b = J / 6.0 * (eta_x * eta_x + eta_y * eta_y + eta_z * eta_z);
//     double c = J / 6.0 * (zeta_x * zeta_x + zeta_y * zeta_y + zeta_z * zeta_z);
//
//     double d = J / 6.0 * (xi_x * eta_x + xi_y * eta_y + xi_z * eta_z);
//     double e = J / 6.0 * (xi_x * zeta_x + xi_y * zeta_y + xi_z * zeta_z);
//     double f = J / 6.0 * (eta_x * zeta_x + eta_y * zeta_y + eta_z * zeta_z);
//
//     RSTLMatrix compound(10, 10);
//     compound = a * (*Tet10_u_xi2) + b * (*Tet10_u_eta2) + c * (*Tet10_u_zeta2)
//         + (2.0*d) * (*Tet10_u_xi_u_eta)
//         + (2.0*e) * (*Tet10_u_xi_u_zeta)
//         + (2.0*f) * (*Tet10_u_eta_u_zeta);
//     for (uint i = 0; i < dim; i++){
//       for (uint j = 0; j < dim; j++){
//         //** * 6.0 because a b c d e f / 6
//         mat_[i][j] = compound[i][j] * 6.0;
//       }
//     }
// //    std::cout << " 2 " << *this << std::endl;
//
//     break;
//   }
        ux2uy2uz2(cell,
                  IntegrationRules::instance().tetWeights(2),
                  IntegrationRules::instance().tetAbscissa(2), false);   break;
    case MESH_HEXAHEDRON_RTTI:
        ux2uy2uz2(cell,
                  IntegrationRules::instance().hexWeights(2),
                  IntegrationRules::instance().hexAbscissa(2), false);   break;
    case MESH_HEXAHEDRON20_RTTI:
        ux2uy2uz2(cell,
                  IntegrationRules::instance().hexWeights(4),
                  IntegrationRules::instance().hexAbscissa(4), false);   break;
    case MESH_TRIPRISM_RTTI:
        ux2uy2uz2(cell,
                  IntegrationRules::instance().priWeights(2),
                  IntegrationRules::instance().priAbscissa(2), false);   break;
    case MESH_TRIPRISM15_RTTI:
        ux2uy2uz2(cell,
                  IntegrationRules::instance().priWeights(4),
                  IntegrationRules::instance().priAbscissa(4), false);   break;

    default:
        std::cerr << cell.rtti() << std::endl;
        THROW_TO_IMPL
        break;
    }

    if (useCache) cell.setUxCache(mat_);

    return *this;
}

template < class ValueType >
void ElementMatrix < ValueType >::fillIds(const MeshEntity & ent, Index nC){
    Index nDims = 1;
    Index nNodes = ent.nodeCount();

    // __MS(nC)
    if (this->_nDof > 0){
        nDims = ent.dim();
        if (nC > ent.dim()){
            //THROW_TO_IMPL
        }

        if (size() != nNodes * nDims) resize(nNodes * nDims);

        for (Index dim = 0; dim < nDims; dim++){
            for (Index i = 0; i < nNodes; i ++) {
                _ids[dim * nNodes + i]  = dim * this->_nDof + ent.node(i).id();
                _idsC[dim * nNodes + i] = _ids[i + dim * nNodes];
                _idsR[dim * nNodes + i] = _ids[i + dim * nNodes];
            }
        }
    } else {
        nDims = 1; //nC;
        this->resize(nNodes*nC, nNodes);
        // mat_.resize(nNodes*nC, nNodes);
        // _ids.resize(nNodes*nC);

        for (Index dim = 0; dim < nDims; dim++){
            for (Index i = 0; i < nNodes; i ++) {
                _ids[dim * nNodes + i] = ent.node(i).id();
                _idsC[dim * nNodes + i] = _ids[dim * nNodes + i];
                _idsR[i] = _ids[dim * nNodes + i];
            }
        }
    }

    // __MS(_ids)
    // __MS(_idsR)
    // __MS(_idsC)

    *this *= 0.0;
}

template void ElementMatrix < double >::fillIds(const MeshEntity & ent, Index Crows);


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

    for (Index r = 0; r < row_.size(); r ++ ){
        double s = 0.0;
        const RMatrix & mat = mat_[r];
        const IndexArray & idx = _ids[r];
        for (Index i = 0; i < mat.rows(); i ++) {
            double t = 0;
            for (Index j = 0; j < mat.cols(); j ++) {
                t += mat[i][j] * (a[idx[j]]-b[idx[j]]);
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
        const RMatrix & mat = mat_[r];
        const IndexArray & idx = _ids[r];
        for (Index i = 0; i < mat.rows(); i ++) {
            double t = 0;
            for (Index j = 0; j < mat.cols(); j ++) {
                t += mat[i][j] * (a[idx[j]]);
            }
            s += t * (b[idx[i]]);
        }

        ret[row_[r]] += s;
    }

    return ret;
}


} // namespace GIMLI
