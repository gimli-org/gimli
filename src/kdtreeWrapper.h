/***************************************************************************
 *   Copyright (C) 2006-2014 by the GIMLi development team       *
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

#ifndef _GIMLI_KDTREEWRAPPER__H
#define _GIMLI_KDTREEWRAPPER__H

#include "gimli.h"

#include <iostream>

#ifdef _MSC_VER
	// warning from kdtree
	#pragma warning( disable: 4396)
#endif

#define KDTREE_DEFINE_OSTREAM_OPERATORS
#include <kdtree++/kdtree.hpp>

typedef KDTree::KDTree< 3, GIMLI::Node *, std::pointer_to_binary_function< GIMLI::Node *, size_t, double > > NodeKDTree;

namespace GIMLI{

//! Interface class for a kd-search tree. We use it for fast nearest neighbor point search in three dimensions.
/*! Interface class for a kd-search tree. Currently we use libkdtree++ from: http://libkdtreeplus-pplus-p.sourcearchive.com/
We use it for fast nearest neighbor point search in three dimensions. The tree is designed to cooperate with \ref Mesh thus it has to be feeded by pointers of \ref Node. */
class KDTreeWrapper{
public:
    /*! Standard constructor */
    KDTreeWrapper();
    
    /*! Standard destructor */
    ~KDTreeWrapper();

    /*! Insert new node to the tree */
    void insert( Node * node );

    /*! Find the nearest \ref Node to the coordinates pos. */
    Node * nearest( const RVector3 & pos );

    /*! Return the amount of nodes inside the tree. */
    uint size( ) const;

    /*! Return a pointer to the base libkdetree++  */
    NodeKDTree * tree();

protected:
    NodeKDTree * tree_;
};

} // namespace GIMLI

#endif // _GIMLI_KDTREEWRAPPER__H
