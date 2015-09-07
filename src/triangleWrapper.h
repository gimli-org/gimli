/***************************************************************************
 *   Copyright (C) 2005-2014 by the resistivity.net development team       *
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

#ifndef _GIMLI_TRIANGLEWRAPPER__H
#define _GIMLI_TRIANGLEWRAPPER__H

#include "gimli.h"

#include <vector>
#include <string>

struct triangulateio;

namespace GIMLI{

class DLLEXPORT TriangleWrapper{
public:
    /*! Construct a PLC from nodes and edges in inMesh and initialize
     * the trianglewrapper. Add region and hole marker later.*/
    TriangleWrapper(const Mesh & inMesh);

    /*! Constructor initialize with input mesh
     * (nodes and edges where used for input and also region marker,
     * if added to the mesh) and generate the 2d mesh with corresponding
     * triSwitches */
    TriangleWrapper(const Mesh & inMesh, Mesh & outMesh,
                    const std::string & triSwitches);

    /*! Default destructur */
    virtual ~TriangleWrapper();

    /*! Set the triangle commandline switches */
    void setSwitches(const std::string & s);

    /*! Return the triangle switches. */
    inline const std::string & switches() const { return switches_; }

    /*! Generate the mesh and store in mesh. */
    void generate(Mesh & mesh);

    /*! Generate and return the new mesh. */
    Mesh generate();

protected:
    /*! For internal use only. */
    void init_();

    /*! For internal use only. */
    void transformTriangleToMesh_(const triangulateio & trimesh, Mesh & mesh);

    /*! For internal use only. Only Edges and nodes(x,y) from mesh are used. */
    void transformMeshToTriangle_(const Mesh & mesh, triangulateio & trimesh);

    /*! For internal use only. */
    void allocateOutMemory_();

    /*! For internal use only. */
    void freeMemory_();

    struct triangulateio * mesh_input_;
    struct triangulateio * mesh_output_;
    struct triangulateio * mesh_voronoi_output_;

    std::string switches_;

    const Mesh * inMesh_;

};

}  // namespace GIMLI

#endif // _GIMLI_TRIANGLEWRAPPER__H
