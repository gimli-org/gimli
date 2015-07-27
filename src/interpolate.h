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

#ifndef _GIMLI_INTERPOLATE__H
#define _GIMLI_INTERPOLATE__H

#include "gimli.h"
#include "matrix.h"
#include <vector>

namespace GIMLI{

/*!Interpolate a given input data regarding the mesh srcMesh to a set of positions 
 * and write the interpolated data to outMat. outMat will resized if necessary.
 * Each data vector in inMat have to correspond to mesh.nodeCount().
 * If data length is mesh.cellCount() \ref cellDataToPointData will performed. 
 * The interpolation rule depend on the shape functions of mesh cells. 
 * Several utility or shortcut functions are defined. */
DLLEXPORT void interpolate(const Mesh & srcMesh, const RMatrix & inMat,
                           const R3Vector & destPos, RMatrix & outMat,
                           bool verbose=false);

/*! Utility function for interpolation. */
DLLEXPORT void interpolate(const Mesh & srcMesh, const RVector & inVec,
                           const R3Vector & destPos, RVector & outVec,
                           bool verbose=false);

/*! Utility function for interpolation. */
DLLEXPORT RVector interpolate(const Mesh & srcMesh, const RVector & inVec,
                              const R3Vector & destPos,
                              bool verbose=false);

/*! Utility function for interpolation. */
DLLEXPORT void interpolate(const Mesh & srcMesh, const RVector & inVec,
                           const Mesh & destMesh, RVector & outVec,
                           bool verbose=false);

/*! Utility function for interpolation. Read in data from fileName and add the 
 interpolated data into the destination mesh. */
DLLEXPORT void interpolate(const Mesh & srcMesh, const std::string & fileName,
                           Mesh & destMesh, bool verbose=false);

/*! Utility function for interpolation.*/
DLLEXPORT RVector interpolate(const Mesh & srcMesh, const RVector & inVec,
                              const RVector & x, const RVector & y, 
                              const RVector & z, bool verbose=false);

/*! Utility function for interpolation.*/
DLLEXPORT RVector interpolate(const Mesh & srcMesh, const RVector & inVec,
                              const RVector & x, const RVector & y, 
                              bool verbose=false);

/*! Utility function for interpolation.*/
DLLEXPORT RVector interpolate(const Mesh & srcMesh, const RVector & inVec,
                              const RVector & x, 
                              bool verbose=false);

/*! Utility function for interpolation. 
 * Interpolate all export data from srcMesh to the destination mesh. 
 * Point and Cell based. */
DLLEXPORT void interpolate(const Mesh & srcMesh, Mesh & destMesh,
                           bool verbose=false);

/*! Utility function for interpolation. 
 * Interpolate the z-coordinate from mesh to the z-coordinate of the 
 * destination mesh.*/
DLLEXPORT void interpolateSurface(const Mesh & srcMesh, Mesh & destMesh, 
                                  bool verbose=false);

/*! Utility function. Convert cell data to point data with the 
 * corresponding the cell interpolation function */
DLLEXPORT RVector cellDataToPointData(const Mesh & mesh,
                                      const RVector & cellData);

DLLEXPORT void triangleMesh_(const Mesh & mesh, Mesh & tmpMesh);
      
} // namespace GIMLI


#endif // _GIMLI_INTERPOLATE__H
