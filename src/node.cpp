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

namespace GIMLI{

std::ostream & operator << ( std::ostream & str, const GIMLI::Node & n ){
    str << "Node: "<< &n << " id: " << n.id() << "\t" << n.pos();
    str << " marker: " << n.marker();
    return str;
}

void Node::smooth( uint function ){
    std::set< Node * > common( commonNodes( this->boundSet() ) );
    //** Achtung konkave gebiete koennen entstehen wenn zu festen knoten benachbarte gesmooth werden
    //** aufzeichen -> pruefen -> fixen.
    //** was passiert bei node at interface or boundary
    RVector3 c( 0.0, 0.0, 0.0 );
    for ( std::set< Node * >::iterator it = common.begin(); it != common.end(); it ++){
        c += (*it)->pos();
    }
    this->setPos( c / (double)common.size() );
}

} // namespace GIMLI{

