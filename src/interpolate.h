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

#ifndef _GIMLI_INTERPOLATE__H
#define _GIMLI_INTERPOLATE__H

#include "gimli.h"
#include "matrix.h"
#include <vector>

namespace GIMLI{

/*! Interpolate a given data vector on a mesh to a set of positions and write the interpolated data to iData. iData will resized if nessecary. Amount of data have to correspond to mesh.nodeCount() if data correspond to mesh.cellCount() cellDataToPointData will performed. The interpolation rule depend on the shape functions of mesh cells. Several utility or shortcut functions are defined.*/
DLLEXPORT void interpolate(const Mesh & mesh, const RMatrix & data,
                           const std::vector< RVector3 > & pos, RMatrix & iData,
                           bool verbose = false);

/*! Utility function for interpolation. */
DLLEXPORT void interpolate(const Mesh & mesh, const RVector & data,
                           const std::vector< RVector3 > & pos, RVector & iData,
                           bool verbose = false);

/*! Utility function for interpolation. */
DLLEXPORT RVector interpolate(const Mesh & mesh, const RVector & data,
                              const std::vector< RVector3 > & pos,
                              bool verbose = false);

/*! Utility function for interpolation. */
DLLEXPORT void interpolate(const Mesh & mesh, const RVector & data,
                           const Mesh & pos, RVector & iData,
                           bool verbose = false);

DLLEXPORT void interpolate(const Mesh & mesh, const std::string & data,
                           Mesh & pos, bool verbose = false);

/*! Utility function for interpolation.*/
DLLEXPORT RVector interpolate(const Mesh & mesh, const RVector & data,
                              const RVector & x, const RVector & y, 
                              const RVector & z, bool verbose = false);

/*! Utility function for interpolation. Interpolate all export data from mesh to the query mesh. Point and Cell based. */
DLLEXPORT void interpolate(const Mesh & mesh, Mesh & qmesh,
                           bool verbose=false);

/*! Utility function for interpolation. Interpolate the z-coordinate from mesh to the z-coordinate of the query mesh qmesh.*/
DLLEXPORT void interpolateSurface(const Mesh & mesh, Mesh & qmesh, 
                                  bool verbose=false);

/*! Utility function. Convert cell data to point data with the corresponding the cell interpolation function */
DLLEXPORT RVector cellDataToPointData(const Mesh & mesh,
                                      const RVector & cellData);

DLLEXPORT void triangleMesh_(const Mesh & mesh, Mesh & tmpMesh);
      
//double interpolate(const RVector3 & queryPos, const MeshEntity & entity, const RVector & sol);

} // namespace GIMLI


#endif // _GIMLI_INTERPOLATE__H
