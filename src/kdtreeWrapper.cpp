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

#include "kdtreeWrapper.h"

#include "node.h"
   
namespace GIMLI{
    
inline double tac( Node * n, size_t i ) { return n->pos()[ i ]; }
    
KDTreeWrapper::KDTreeWrapper(){
    tree_ = new NodeKDTree( std::ptr_fun( tac ) );
}

KDTreeWrapper::~KDTreeWrapper(){
    if ( tree_ ) delete tree_;
}

void KDTreeWrapper::insert( Node * node ){ 
    tree_->insert( node ); 
}

Node * KDTreeWrapper::nearest( const RVector3 & pos ){
    Node testNode( pos );
    return *tree_->find_nearest( &testNode ).first;
}

uint KDTreeWrapper::size( ) const{
    return tree_->size();
}

NodeKDTree * KDTreeWrapper::tree() { 
    return tree_; 
}

} // namespace GIMLI
