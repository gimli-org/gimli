/***************************************************************************
 *   Copyright (C) 2006-2011 by the resistivity.net development team       *
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

#include "node.h"

#include "meshentities.h"
#include "shape.h"

namespace GIMLI{

std::ostream & operator << (std::ostream & str, const GIMLI::Node & n){
    str << "Node: "<< &n << " id: " << n.id() << "\t" << n.pos();
    str << " marker: " << n.marker();
    return str;
}

Node::Node(){
    init_();
    marker_ = 0;
}

Node::Node(double x, double y, double z)
: pos_(RVector3(x, y, z)) {
    init_();
}

Node::Node(const RVector3 & pos)
: pos_(pos) {
    init_();
}

Node::Node(const RVector3 & pos, int marker, int id)
: pos_(pos) {
    init_();
    setMarker(marker);
    setId(id);
}

Node::Node(const Node & node){
    copy_(node);
}

Node & Node::operator = (const Node & node){
    if (this != &node){
        copy_(node);
    }
    return *this;
}

Node::~Node(){
    //std::cout << " delete Node " << pos_ << " " << id_ << " at " << this << std::endl;
}

void Node::changed_(){
    for (std::set < Boundary * >::iterator it = boundSet_.begin(); it!= boundSet_.end(); it ++){
        (*it)->shape().changed();        
    }
    for (std::set < Cell * >::iterator it = cellSet_.begin(); it!= cellSet_.end(); it ++){
        (*it)->shape().changed();        
    }
}

void Node::copy_(const Node & node){
    //std::cout << "copy node from "  << & node << " into " << this << std::endl;
    pos_    = node.pos();
    marker_ = node.marker();
    setId(node.id());
}

void Node::init_(){
    setId(-1);
}

void Node::smooth(uint function){
    std::set< Node * > common(commonNodes(this->boundSet()));
    //** Achtung konkave gebiete koennen entstehen wenn zu festen knoten benachbarte gesmooth werden
    //** aufzeichen -> pruefen -> fixen.
    //** was passiert bei node at interface or boundary
    RVector3 c(0.0, 0.0, 0.0);
    for (std::set< Node * >::iterator it = common.begin(); it != common.end(); it ++){
        c += (*it)->pos();
    }
    this->setPos(c / (double)common.size());
}

} // namespace GIMLI{

