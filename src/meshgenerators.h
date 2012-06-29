/***************************************************************************
 *   Copyright (C) 2008-2011 by the resistivity.net development team       *
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

#ifndef GIMLI_MESHGENERATORS__H
#define GIMLI_MESHGENERATORS__H

#include "gimli.h"

namespace GIMLI{

/*! Generate simple one dimensional mesh with nCells cells with length = 1.0, and nCells + 1 nodes. */
DLLEXPORT Mesh createMesh1D( uint nCells, uint nClones = 1 );

/*! Generate simple one dimensional mesh with nodes at position in RVector pos. */
DLLEXPORT Mesh createMesh1D( const std::vector < double > & x );
DLLEXPORT Mesh createMesh1D( const RVector & x );

/*! Generate 1D block model of thicknesses and properties */
DLLEXPORT Mesh createMesh1DBlock( uint nLayers, uint nProperties = 1 );

/*! Generate simple two dimensional mesh with nRows x nCols cells with each length = 1.0 */
DLLEXPORT Mesh createMesh2D( uint xDim, uint yDim, int markerType = 0 );

/*! Generate simple two dimensional mesh with nodes at position in RVector x and y. */
DLLEXPORT Mesh createMesh2D( const RVector & x, const RVector & y, int markerType = 0 );

/*! Generate simple three dimensional mesh with nx x nx x nz cells with each length = 1.0 */
DLLEXPORT Mesh createMesh3D( uint xDim, uint yDim, uint zDim, int markerType = 0 );

/*! Generate simple three dimensional mesh with nodes at position in RVector x and y. */
DLLEXPORT Mesh createMesh3D( const RVector & x, const RVector & y, const RVector & z, int markerType = 0 );

/*! Add triangle boundary to the mesh. Return false on failors. */
DLLEXPORT bool addTriangleBoundary( Mesh & mesh, 
                                    double xBoundary, double yBoundary, int cellMarker, bool save = false );

/*! Generate a simple three dimensional mesh by extruding a two dimensional triangle mesh into RVector z using triangle prism */
DLLEXPORT Mesh createMesh3D( const Mesh mesh, const RVector & z, int markerType = 0 );

    
} // namespace GIMLI

#endif // GIMLI_MESHGENERATORS__H
