/******************************************************************************
 *   Copyright (C) 2006-2022 by the GIMLi development team                    *
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
