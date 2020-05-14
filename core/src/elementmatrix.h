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

#ifndef _GIMLI_ELEMENTMATRIX__H
#define _GIMLI_ELEMENTMATRIX__H

#include "gimli.h"
#include "vector.h"
#include "matrix.h"

namespace GIMLI{

template < class ValueType > class DLLEXPORT ElementMatrix {
public:
    /*! If dof != 0 then scalar field approximation is to be supposed.
    For vector field solution give a dof, means be the number of nodes of the current mesh. */
    ElementMatrix(Index dof=0) {
        this->_nDof = dof;
    }

    ~ElementMatrix() {}

    inline const Vector< ValueType > & operator[](Index row) const {
        return mat_[row]; }

    void resize(Index rows, Index cols=0) {
        if (cols == 0) cols = rows;
        _idsR.resize(rows);
        _idsC.resize(cols);
        _ids.resize(rows);
        mat_.resize(rows, cols);
    }

    ElementMatrix < ValueType > & operator += (const ElementMatrix < ValueType > & E){
        for (uint i = 0; i < size(); i ++){ mat_[i] += E.row(i); }
        return *this;
    }

    #define DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(OP)                   \
        ElementMatrix < ValueType > & operator OP##= (ValueType val) { \
            for (Index i = 0; i < size(); i ++) mat_[i] OP##= val; \
            return *this;\
        } \

        DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(+)
        DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(-)
        DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(/)
        DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(*)

    #undef DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__

    inline Index size() const { return mat_.rows(); }
    inline Index rows() const { return mat_.rows(); }
    inline Index cols() const { return mat_.cols(); }

    inline const ValueType & getVal(Index i, Index j) const {
        return mat_[i][j]; }
    inline void setVal(Index i, Index j, const ValueType & v) {
        mat_[i][j] = v; }
    inline void addVal(Index i, Index j, const ValueType & v) {
        mat_[i][j] += v; }

    /*! Set data matrix. */
    inline void setMat(const Matrix < ValueType > & m) { mat_ = m; }
    /*! Return data matrix. */
    inline const Matrix < ValueType > & mat() const { return mat_; }
    /*! Return data for row i. */
    inline const Vector < ValueType > & row(Index i) const { return mat_[i]; }
    /*! Return data for col i. */
    inline Vector < ValueType > col(Index i) const {
        Vector < ValueType > ret(this->rows());
        for (Index j = 0; j < ret.size(); j ++){ ret[j] = mat_[j][i];}
        return ret;
    }

    /*! Fill the node ids with a number of coefficents.
    For vector field approximation give field dimension 2 or 3.
    Please note that you need to give the number of nodes to the ElementMatrix constructor.
    */
    void fillIds(const MeshEntity & ent, Index nC=1);

    /*! Set all node indices for row and columns. Can be unsymmetric.*/
    inline void setIds(const IndexArray & idsR, const IndexArray & idsC) {
        _idsR = idsR; _idsC = idsC; _ids = idsR;
    }

    /*! Set all node indices.*/
    inline void setIds(const IndexArray & ids) {
        _idsR = ids; _idsC = ids; _ids = ids;
    }
    /*! Return all node indices.*/
    inline const IndexArray & ids() const { return _ids; }
    /*! Return all row node indices.*/
    inline const IndexArray & rowIDs() const { return _idsR; }
    /*! Return all column node indices.*/
    inline const IndexArray & colIDs() const { return _idsC; }
    /*! Return the node index for node i.*/
    inline const Index idx(Index i) const {
        return _ids[i]; }

    /*! Fill this element matrix with int_boundary C * u */
    ElementMatrix < ValueType > & u(const MeshEntity & ent
                                    // const Matrix< ValueType > & C
                                    );

    /*! Fill this element matrix with int_domain C * u * u */
    ElementMatrix < ValueType > & u2(const MeshEntity & ent
                                    //  const Matrix< ValueType > & C
                                     );

    /*! Get integration weights and points for the entity. */
    void getWeightsAndPoints(const MeshEntity & ent,
                const RVector * &w, const R3Vector * &x, int order);


    /*! Return the stress matrix for this entity.*/
    Vector < ValueType > stress(const MeshEntity & ent,
                                const Matrix< ValueType > & C,
                                const RVector & u, bool voigtNotation=false);

    /*! Return gradient base for the last entity, i.e., integrate over gradient u over the cell. */
    const Matrix < ValueType > & gradientBase() const { return this->_grad;}

    /*! Fill Element Gradients Matrix for all integration points. */
    void fillGradientBase(const MeshEntity & ent,
                          const RVector & w,
                          const R3Vector & x,
                          Index nC,
                          bool voigtNotation);

    /*! Fill this element matrix with int_domain C * grad u*/
    ElementMatrix < ValueType > & gradU(const Cell & cell,
                                        Index nC,
                                        bool voigtNotation=false);

    /*! Fill this element matrix with int_domain C * grad u */
    ElementMatrix < ValueType > & gradU(const MeshEntity & ent,
                                        const RVector & w,
                                        const R3Vector & x,
                                        Index nC,
                                        bool voigtNotation=false);


    /*! Fill this element matrix with int_domain C * grad u * grad u.
    For scalar field approximation define C.size() = (1x1) isotropic or anisotropic
    C.size() = (cell.dim() x cell.dim()) parameter. For vector field approximation
    create the ElementMatrix with appropriate dof and C can be of size = (1x1)
    for isotropic, (cell.dim() x cell.dim()) for anisotropic, or (3x3)
    for 2D elastic and (6x6) 3D elastic parameter.
    Notation for elastic parameters Kelvin notation as default and can be Voigt's notation if needed.
     */
    ElementMatrix < ValueType > & gradU2(const Cell & cell,
                                         const Matrix< ValueType > & C,
                                         bool voigtNotation=false);

    /*! Fill this element matrix with int_domain C * grad u * grad u.*/
    ElementMatrix < ValueType > & gradU2(const Cell & cell, ValueType c){
        Matrix < ValueType > C(1, 1);
        C[0][0] = c;
        return this->gradU2(cell, C, false);
    }

    /*! Fill this element matrix with int_domain C * grad u * grad u. */
    ElementMatrix < ValueType > & gradU2(const MeshEntity & ent,
                                         const Matrix< ValueType > & C,
                                         const RVector & w,
                                         const R3Vector & x,
                                         bool voigtNotation=false);

    ElementMatrix < ValueType > & ux2uy2uz2(const Cell & cell,
                                            bool useCache=false);

    ElementMatrix < ValueType > & u(const MeshEntity & ent,
                                    const RVector & w,
                                    const R3Vector & x,
                                    bool verbose=false);
    ElementMatrix < ValueType > & u2(const MeshEntity & ent,
                                     const RVector & w,
                                     const R3Vector & x,
                                     bool verbose=false);
    ElementMatrix < ValueType > & ux2(const MeshEntity & ent,
                                      const RVector & w,
                                      const R3Vector & x,
                                      bool verbose=false);
    ElementMatrix < ValueType > & ux2uy2(const MeshEntity & ent,
                                         const RVector & w,
                                         const R3Vector & x,
                                         bool verbose=false);
    ElementMatrix < ValueType > & ux2uy2uz2(const MeshEntity & ent,
                                            const RVector & w,
                                            const R3Vector & x,
                                            bool verbose=false);

    ElementMatrix < double > & dudi(const MeshEntity & ent,
                                  const RVector & w,
                                  const R3Vector & x,
                                  Index i, bool verbose=false);

    ElementMatrix < double > & ux(const MeshEntity & ent,
                                  const RVector & w,
                                  const R3Vector & x,
                                  bool verbose=false){
        return dudi(ent, w, x, 0, verbose);
    }
    ElementMatrix < double > & uy(const MeshEntity & ent,
                                  const RVector & w,
                                  const R3Vector & x,
                                  bool verbose=false){
        return dudi(ent, w, x, 1, verbose);
    }
    ElementMatrix < double > & uz(const MeshEntity & ent,
                                  const RVector & w,
                                  const R3Vector & x,
                                  bool verbose=false){
        return dudi(ent, w, x, 2, verbose);
    }

    Vector < ValueType > mult(const Vector < ValueType > & v){
        Vector < ValueType > ret(this->size());
        this->mult(v, ret);
        return ret;
    }

    /*! Return S * a */
    void mult(const Vector < ValueType > & a, Vector < ValueType > & ret){
        ASSERT_EQUAL(size(), ret.size())
        for (Index i = 0; i < size(); i ++) {
            for (Index j = 0; j < size(); j ++) {
                ret[i] += mat_[i][j] * a[_ids[j]];
            }
        }
    }
    /*! Return (S * a) * b */
    ValueType mult(const Vector < ValueType > & a, const Vector < ValueType > & b){
        ValueType ret = 0;
        for (Index i = 0; i < size(); i ++) {
            ValueType t = 0;
            for (Index j = 0; j < size(); j ++) {
                t += mat_[i][j] * a[_ids[j]];
            }
            ret += t * b[_ids[i]];
        }
        return ret;
    }

    /*! Return (S * (a-b)) * (m-n) TODO Check if its the same like mult(a-b, m-n)*/
    template < class Val > Val mult_(const Vector < Val > & a,
                                     const Vector < Val > & b,
                                     const Vector < Val > & m,
                                     const Vector < Val > & n){
        Val ret = 0;
        for (Index i = 0; i < size(); i ++) {
            Val t = 0;
            for (Index j = 0; j < size(); j ++) {
                t += mat_[i][j] * (a[_ids[j]]-b[_ids[j]]);
            }
            ret += t * (m[_ids[i]]-n[_ids[i]]);
        }
        return ret;
    }

    double mult(const RVector & a, const RVector & b,
                const RVector & m, const RVector & n){
        return mult_(a, b, m, n);
    }
    Complex mult(const CVector & a, const CVector & b,
                const CVector & m, const CVector & n){
        return mult_(a, b, m, n);
    }


protected:
    Matrix < ValueType > mat_;
    IndexArray _ids;
    IndexArray _idsC;
    IndexArray _idsR;

    std::map< uint, RVector > uCache_;
    std::map< uint, Matrix < ValueType > > u2Cache_;

    std::vector< Matrix < ValueType > > _B;
    Matrix < ValueType > _grad;

    // number of single dof
    Index _nDof;
    // number of coefficients: 1, 2, 3 for scalar(dim), 1, 3, 6 for vector(dim)
    // Index _nC;

    RMatrix dNdr_;
    RMatrix dNds_;
    RMatrix dNdt_;

    RMatrix dNdx_; // (nRules, nVerts)
    RMatrix dNdy_; // (nRules, nVerts)
    RMatrix dNdz_; // (nRules, nVerts)

    RMatrix _abaTmp; // temp workspace

private:
    /*! No copy operator. */
    ElementMatrix(const ElementMatrix < ValueType > & E) {
        std::cout << "ElementMatrix(const ElementMatrix < ValueType > & E) " << std::endl;
        THROW_TO_IMPL
//     this->resize(E.size());
//         for (uint i = 0; i < E.size(); i ++) mat_[i] = E.mat(i);
//         _ids = E.idx();
//         initBaseMatricies();
    }
    /*! No assignment operator. */
    ElementMatrix < ValueType > & operator = (const ElementMatrix < ValueType > & E) {
        std::cout << "ElementMatrix::operator = (" << std::endl;
        THROW_TO_IMPL
        if (this != & E){
//             this->resize(E.size());
//             for (uint i = 0; i < E.size(); i ++) mat_[i] = E.row(i);
//             _ids = E.idx();
        } return *this;
    }
};

class DLLEXPORT ElementMatrixMap {
public:

    void add(Index row, const ElementMatrix < double > & Ai);

    //TODO .. check if its the same like mult(a-b, m-n))
    RVector mult(const RVector & a, const RVector & b,
                 const RVector & m, const RVector & n) const;

    /*! Return (S_i * a) * b for all i*/
    RVector mult(const RVector & a, const RVector & b) const;

    Index rows() const { return rows_; }

    Index cols() const { return cols_; }


protected:
    std::vector< RMatrix > mat_;
    std::vector< IndexArray > _ids;
    std::vector< Index > row_;

    Index rows_;
    Index cols_;
};

template < class ValueType > std::ostream & operator << (std::ostream & str,
                                                         const ElementMatrix< ValueType > & e);

template < > DLLEXPORT std::ostream & operator << (std::ostream & str,
                                                    const ElementMatrix< double > & e);


} // namespace GIMLI{

#endif // _GIMLI_ELEMENTMATRIX__H
