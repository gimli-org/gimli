/***************************************************************************
 *   Copyright (C) 2006-2015 by the resistivity.net development team       *
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

#ifndef _GIMLI_ELEMENTMATRIX__H
#define _GIMLI_ELEMENTMATRIX__H

#include "gimli.h"
#include "vector.h"
#include "matrix.h"

namespace GIMLI{

template < class ValueType > class DLLEXPORT ElementMatrix {
public:
    ElementMatrix() { }

    ElementMatrix(const ElementMatrix < ValueType > & E) {
        std::cout << "ElementMatrix(const ElementMatrix < ValueType > & E) " << std::endl;
        THROW_TO_IMPL
//     this->resize(E.size());
//         for (uint i = 0; i < E.size(); i ++) mat_[i] = E.mat(i);
//         idx_ = E.idx();
//         initBaseMatricies();
    }

    ElementMatrix < ValueType > & operator = (const ElementMatrix < ValueType > & E) {
        std::cout << "ElementMatrix::operator = (" << std::endl;
        THROW_TO_IMPL
        if (this != & E){
//             this->resize(E.size());
//             for (uint i = 0; i < E.size(); i ++) mat_[i] = E.row(i);
//             idx_ = E.idx();
        } return *this;
    }

    ~ElementMatrix() {}

    inline const Vector< ValueType > & operator[](uint row) const { return mat_[row]; }
    
    void resize(uint newSize) {
        idx_.resize(newSize);
        mat_.resize(newSize, newSize);
    }

    ElementMatrix < ValueType > & operator += (const ElementMatrix < ValueType > & E){
        for (uint i = 0; i < size(); i ++){ mat_[i] += E.row(i); } return *this;
    }
  
    #define DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(OP)                   \
        ElementMatrix < ValueType > & operator OP##= (ValueType val) { \
            for (register size_t i = 0; i < size(); i ++) mat_[i] OP##= val; \
            return *this;\
        } \

        DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(+)
        DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(-)
        DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(/)
        DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(*)

    #undef DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__

    inline const Index idx(Index i) const { return idx_[i]; }
    inline Index size() const { return mat_.rows(); }
    inline const ValueType & getVal(Index i, Index j) const { return mat_[i][j]; }

    inline const Vector < ValueType > & row(Index i) const { return mat_[i]; }
    inline const Matrix < ValueType > & mat() const { return mat_; }
    inline const IndexArray & idx() const { return idx_; }

    ElementMatrix < ValueType > & u(const MeshEntity & ent);

    ElementMatrix < ValueType > & u2(const MeshEntity & ent);

    ElementMatrix < ValueType > & ux2uy2uz2(const Cell & cell, bool useCache=false);

    ElementMatrix < ValueType > & u(const MeshEntity & ent, 
                            const RVector & w,
                            const std::vector < RVector3 > & integrationPnts,
                            bool verbose = false);
    ElementMatrix < ValueType > & u2(const MeshEntity & ent, const RVector & w, 
                                     const std::vector < RVector3 > & integrationPnts,
                                     bool verbose = false);
    ElementMatrix < ValueType > & ux2(const MeshEntity & ent, const RVector & w, 
                               const std::vector < RVector3 > & integrationPnts,
                              bool verbose = false);
    ElementMatrix < ValueType > & ux2uy2(const MeshEntity & ent, const RVector & w, 
                                  const std::vector < RVector3 > & integrationPnts,
                                 bool verbose = false);
    ElementMatrix < ValueType > & ux2uy2uz2(const MeshEntity & ent, const RVector & w, 
                                    const std::vector < RVector3 > & integrationPnts, 
                                    bool verbose = false);

protected:
    //RMatrix mat_;
    Matrix < ValueType > mat_;
    IndexArray idx_;

    RMatrix functx_;
    RMatrix functy_;
    RMatrix functz_;

    std::map< uint, RVector > uCache_;
    std::map< uint, Matrix < ValueType > > u2Cache_;
    
    RMatrix dNdr_;
    RMatrix dNds_;
    RMatrix dNdt_;

};

template < class ValueType > std::ostream & operator << (std::ostream & str, const ElementMatrix< ValueType > & e){
    for (uint i = 0; i < e.idx().size(); i ++) str << e.idx(i) << " " ;

    str << std::endl;
    for (uint i = 0; i < e.size(); i ++){
        str << e.idx(i) << "\t: ";
        for (uint j = 0; j < e.size(); j ++){
            str << e.getVal(i , j) << " ";
        }
        str << std::endl;
    }
    return str;
}

    
} // namespace GIMLI{

#endif // _GIMLI_ELEMENTMATRIX__H
