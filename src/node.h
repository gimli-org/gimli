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

#ifndef _GIMLI_NODE__H
#define _GIMLI_NODE__H

#include "gimli.h"

#include "baseentity.h"
#include "pos.h"

#include <set>

namespace GIMLI{

//! 3D Node
/*!
 * Node is a basic entity of a mesh at a 3D position x/y/z (a vertex),
 * which owns a marker, an id, and information about connected boundarys and cells.
 * For a 2D node ignore $z$ or let $z$-koord = 0.0.
 */
class DLLEXPORT Node : public BaseEntity {

public:
    /*! Construct node with non valid position and marker = 0 */
    Node();
    
    /*! Construct node from koordinates $x,y,[z]$ with marker = 0 */
    Node(double x, double y, double z=0.0);

    /*! Construct node from RVector3 */
    Node(const RVector3 & pos);

    /*! Construct node from RVector3 with marker and optional id */
    Node(const RVector3 & pos, int marker, int id=-1);

    /*! Copy constructor */
    Node(const Node & node);

    /*! Assignement operator */
    Node & operator = (const Node & node);

    /*! Destruct the node and all containing informations */
    ~Node();
    
    /*! Unchecked index operator to pos */
    inline double & operator [] (const Index i) { return this->at(i); }

    /*! Unchecked index operator to pos */
    inline const double & operator [] (const Index i) const { return this->at(i); }
    
    inline uint rtti() const { return MESH_NODE_RTTI; }

    inline void setPos(const RVector3 & pos) { pos_ = pos; }

    inline const RVector3 & pos() const { return pos_; }

    inline RVector3 & pos() { return pos_; }

    inline void insertBoundary(Boundary & bound){ boundSet_.insert(&bound); }

    inline void eraseBoundary(Boundary & bound){ boundSet_.erase(&bound); }

    inline void insertCell(Cell & cell){ cellSet_.insert(&cell); }

    inline void eraseCell(Cell & cell){ cellSet_.erase(&cell); }

    inline const std::set < Boundary * > & boundSet() const { return boundSet_; }

    inline std::set < Boundary * > & boundSet() { return boundSet_; }

    inline const std::set < Cell * > & cellSet() const { return cellSet_; }

    inline std::set < Cell * > & cellSet() { return cellSet_; }

    inline void scale(const RVector3 & s) { changed_(); pos_.scale(s); }

    inline void translate(const RVector3 & t) { changed_(); pos_.translate(t); }

    inline void rotate(const RVector3 & r) { changed_(); pos_.rotate(r); }

    inline const double & x() const { return pos_[0]; }

    inline const double & y() const { return pos_[1]; }

    inline const double & z() const { return pos_[2]; }

    inline double & at(Index i) { return pos_[i]; }
    
    inline const double & at(Index i) const { return pos_[i]; }
    
    inline double dist(const Node & n) const { return pos_.dist(n.pos()); }

    /*!*/
    void smooth(uint function);

protected:

    void copy_(const Node & node);
    
    void changed_();

    void init_();

    RVector3 pos_;

    std::set < Boundary * > boundSet_;
    std::set < Cell * >     cellSet_;

}; // class Node

DLLEXPORT std::ostream & operator << (std::ostream & str, const GIMLI::Node & node);

inline bool operator == (const Node & n1, const Node & n2) { return n1.pos() == n2.pos(); }

} // namespace GIMLI

#endif // _GIMLI_NODE__H
