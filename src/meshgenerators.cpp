/***************************************************************************
 *   Copyright (C) 2008-2015 by the resistivity.net development team       *
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

#include "meshgenerators.h"
#include "mesh.h"
#include "meshentities.h"
#include "triangleWrapper.h"

namespace GIMLI{

Mesh createGrid(const RVector & x){
    Mesh mesh(1);
    mesh.createGrid(x);
    return mesh;
}

Mesh createGrid(const RVector & x, const RVector & y){
    Mesh mesh(2);
    mesh.createGrid(x, y);
    return mesh;
}

Mesh createGrid(const RVector & x, const RVector & y, const RVector & z){
    Mesh mesh(3);
    mesh.createGrid(x, y, z);
    return mesh;
}

Mesh createMesh1D(Index nCells, Index nClones){
    RVector x(nCells * nClones + 1);

    std::generate(x.begin(), x.end(), IncrementSequence< double >(0.0));
    Mesh mesh(createMesh1D(x));
    for (Index i = 0; i < nClones; i ++) {
        for (Index j = 0; j < nCells; j ++) {
            mesh.cell((i * nCells) + j).setMarker(i);
        }
    }

  return mesh;
}

Mesh createMesh1D(const RVector & x){
    Mesh mesh(1);
    mesh.create1DGrid(x);
    return mesh;
}

Mesh createMesh1DBlock(Index nLayers, Index nProperties){
    Index nPar = nLayers * (nProperties + 1);
    RVector x(nPar);
    std::generate(x.begin(), x.end(), IncrementSequence< double>(0.0));
    Mesh mesh(createMesh1D(x));
    /*! Thicknesses have marker 0 */
    for (Index i = 0; i < nLayers - 1; i++) mesh.cell(i).setMarker(0);
    /*! Properties have markers 1,2,... */
    for (Index i = 0; i < nProperties; i ++) {
        for (Index j = 0; j < nLayers; j ++) {
            mesh.cell((i + 1) * nLayers + j -1).setMarker(i + 1);
        }
    }
    return mesh;
}

Mesh createMesh2D(Index xDim, Index yDim, int markerType){
    RVector x(xDim + 1); std::generate(x.begin(), x.end(), IncrementSequence< double >(0.0));
    RVector y(yDim + 1); std::generate(y.begin(), y.end(), IncrementSequence< double >(0.0));
    return createMesh2D(x, y, markerType);
}

Mesh createMesh2D(const RVector & x, const RVector & y, int markerType){
    Mesh mesh(2);
    mesh.create2DGrid(x, y, markerType);
    for (Index i = 0; i < mesh.boundaryCount(); i ++){
        if (!mesh.boundary(i).leftCell() || !mesh.boundary(i).rightCell()){
            mesh.boundary(i).setMarker(1);
        }
    }
    return mesh;
}

Mesh createMesh3D(Index xDim, Index yDim, Index zDim, int markerType){
    RVector x(xDim + 1); std::generate(x.begin(), x.end(), IncrementSequence< double >(0.0));
    RVector y(yDim + 1); std::generate(y.begin(), y.end(), IncrementSequence< double >(0.0));
    RVector z(zDim + 1); std::generate(z.begin(), z.end(), IncrementSequence< double >(0.0));
    return createMesh3D(x, y, z, markerType);
}

Mesh createMesh3D(const RVector & x, const RVector & y, const RVector & z, int markerType){
    Mesh mesh(3);
    mesh.create3DGrid(x, y, z, markerType);
    int marker = 1;
    for (Index i = 0; i < mesh.boundaryCount(); i ++){
        if (!mesh.boundary(i).leftCell() || !mesh.boundary(i).rightCell()){
            mesh.boundary(i).setMarker(marker);
        }
    }
    return mesh;
}

Mesh createMesh2D(const Mesh & mesh, const RVector & y, 
                  int front, int back, bool adjustBack){
    Mesh mesh2(2);
    
    bool first = true;
    for (Index iy = 0; iy < y.size(); iy ++){
        for (Index in = 0; in < mesh.nodeCount(); in ++){
            int marker = 0;
            if (first) marker = mesh.node(in).marker();
            
            mesh2.createNode(mesh.node(in).pos() + RVector3(0.0, y[iy]), marker);
        }
        first = false;
    }
    std::vector < Node * > nodes;
    for (Index iy = 1; iy < y.size(); iy ++){
        first = true;
        for (Index in = 0; in < mesh.boundaryCount(); in ++){
            Index nn = mesh.boundary(in).nodeCount();
            nodes.resize(nn * 2) ;
            
            nodes[0] = & mesh2.node((iy - 1) * mesh.nodeCount() + mesh.boundary(in).node(0).id());
            nodes[1] = & mesh2.node((iy - 1) * mesh.nodeCount() + mesh.boundary(in).node(1).id());
            nodes[2] = & mesh2.node((iy * mesh.nodeCount() + mesh.boundary(in).node(1).id()));
            nodes[3] = & mesh2.node((iy * mesh.nodeCount() + mesh.boundary(in).node(0).id()));
            
            mesh2.createCell(nodes, mesh.boundary(in).marker());
            
            if (iy == 1){
                // create top layer boundaries // in revers direction so the outer normal shows downward into the mesh
            }
            if (iy == y.size()-1){
                // create bottom layer boundaries
            }
        }
        first = false;
    }
    
    
    return mesh;
}

Mesh createMesh3D(const Mesh & mesh, const RVector & z, int topLayer, int bottomLayer){
    Mesh mesh3(3);
    
    if (z.size() < 2){
        std::cout << "Warning!: " << WHERE_AM_I << "extrusion vector size need z be greater than 1" << std::endl;
    }
        
    bool first = true;
    for (Index iz = 0; iz < z.size(); iz ++){
        for (Index ic = 0; ic < mesh.nodeCount(); ic ++){
            int marker = 0;
            if (first) marker = mesh.node(ic).marker();
            
            mesh3.createNode(mesh.node(ic).pos() + RVector3(0.0, 0.0, z[iz]), marker);
        }
        first = false;
    }
    
    std::vector < Node * > nodes;
    
    for (Index iz = 1; iz < z.size(); iz ++){
        first = true;
        for (Index ic = 0; ic < mesh.cellCount(); ic ++){
            Index nC = mesh.cell(ic).nodeCount();
            nodes.resize(nC * 2) ;
            
            for (Index k = 0; k < nC; k ++){
                nodes[k] = & mesh3.node((iz-1) * mesh.nodeCount() + mesh.cell(ic).node(k).id());
            }
            for (Index k = 0; k < nC; k ++){
                nodes[nC + k] = & mesh3.node(iz * mesh.nodeCount() + mesh.cell(ic).node(k).id());
            }
            mesh3.createCell(nodes, mesh.cell(ic).marker());
            
            if (iz == 1){
                // create top layer boundaries // in revers direction so the outer normal shows downward into the mesh
                std::vector < Node * > nBound(nC); for (Index k = 0; k < nC; k ++) nBound[nC - k - 1] = nodes[k];
                mesh3.createBoundary(nBound, topLayer);
            }
            if (iz == z.size()-1){
                // create bottom layer boundaries
                std::vector < Node * > nBound(nC); for (Index k = 0; k < nC; k ++) nBound[k] = nodes[nC + k];
                mesh3.createBoundary(nBound, bottomLayer);
            }
        }
        first = false;
    }
    
    nodes.resize(4);
    for (Index iz = 1; iz < z.size(); iz ++){
        for (Index ib = 0; ib < mesh.boundaryCount(); ib ++){
            if (mesh.boundary(ib).marker() != 0){
                nodes[0] = & mesh3.node((iz-1) * mesh.nodeCount() + mesh.boundary(ib).node(0).id());
                nodes[1] = & mesh3.node((iz-1) * mesh.nodeCount() + mesh.boundary(ib).node(1).id());
                nodes[3] = & mesh3.node(iz * mesh.nodeCount() + mesh.boundary(ib).node(0).id());
                nodes[2] = & mesh3.node(iz * mesh.nodeCount() + mesh.boundary(ib).node(1).id());
                mesh3.createBoundary(nodes, mesh.boundary(ib).marker());
            }
        }
    }
    
    return mesh3;
}

bool addTriangleBoundary(Mesh & mesh, double xBoundary, double yBoundary, int cellMarker, 
                          bool save){
    DEPRECATED //17.09.215
    int boundMarker = -5;
    
    mesh.createNeighbourInfos(true);
    Boundary *b = NULL;
    for (std::vector < Boundary * >::const_iterator it = mesh.boundaries().begin();
        it != mesh.boundaries().end(); it ++){
        b = (*it);
        if (! b->leftCell() || ! b->rightCell()) {
            if (b->marker() != MARKER_BOUND_HOMOGEN_NEUMANN){
                b->setMarker(boundMarker);
            }
        }
    }
    
    std::vector < Boundary * > markedBoundaries(mesh.findBoundaryByMarker(boundMarker));
    Boundary *start(markedBoundaries.front());
    
    std::list < Node * > boundNodes;
    
    boundNodes.push_back(& start->node(0));
    Boundary * thisBoundary = start;
    
    while (thisBoundary != NULL){
        Node * thisNode = boundNodes.back();
        Boundary * nextBoundary = NULL;
    
        for (std::set < Boundary * >::iterator it = thisNode->boundSet().begin();
                it != thisNode->boundSet().end(); it++){
        
            if ((*it)->marker() == boundMarker && (*it) != thisBoundary){
                nextBoundary = (*it);
                break;
            }
        }
        if (nextBoundary){
            Node * nextNode = NULL; 
            if (&nextBoundary->node(0) != thisNode){
                nextNode = &nextBoundary->node(0);
            } else {
                nextNode = &nextBoundary->node(1);
            }
            //** Check closed polygones here 
            if (find(boundNodes.begin(), boundNodes.end(), nextNode) == boundNodes.end()) {
                boundNodes.push_back(nextNode);
            } else {
                break;
            }
        }
        thisBoundary = nextBoundary;
    }
    boundNodes.push_front(& start->node(1));
    thisBoundary = start;
    
    while (thisBoundary != NULL){
        Node * thisNode = boundNodes.front();
        Boundary * nextBoundary = NULL;
    
        for (std::set < Boundary * >::iterator it = thisNode->boundSet().begin();
                it != thisNode->boundSet().end(); it++){
        
            if ((*it)->marker() == boundMarker && (*it) != thisBoundary){
                nextBoundary = (*it);
                break;
            }
        }
        if (nextBoundary){
            Node * nextNode = NULL; 
            if (& nextBoundary->node(0) != thisNode){
                nextNode = & nextBoundary->node(0);
            } else {
                nextNode = & nextBoundary->node(1);
            }

            //** Check closed polygones here 
            if (find(boundNodes.begin(), boundNodes.end(), nextNode) == boundNodes.end()) {
                boundNodes.push_front(nextNode);
            } else {
                break;
            }
        }
        thisBoundary = nextBoundary;
    }
    
    Mesh poly;
    std::vector < Node * > innerBound;
    for (std::list < Node * >::iterator it = boundNodes.begin(); 
        it != boundNodes.end(); it ++){
        innerBound.push_back(poly.createNode((*it)->pos()));
    }
        
    //** looking for polygon ends
    Node * upperLeftSurface = innerBound.front(); 
    Node * upperRightSurface = innerBound.back();
    if (upperLeftSurface->pos()[0] > upperRightSurface->pos()[0]){
        upperLeftSurface = innerBound.back(); 
        upperRightSurface = innerBound.front();
    }
    
    std::vector < Node * > outerBound;
    outerBound.push_back(upperLeftSurface);
    
    Node *n1 = poly.createNode(upperLeftSurface->pos() - RVector3(xBoundary, 0.0));
    Node *n2 = poly.createNode(upperLeftSurface->pos() - RVector3(xBoundary, - mesh.ymin() + yBoundary));
    Node *n3 = poly.createNode(upperRightSurface->pos() - RVector3(-xBoundary, - mesh.ymin() + yBoundary));
    Node *n4 = poly.createNode(upperRightSurface->pos() - RVector3(-xBoundary, 0.0));
    
    outerBound.push_back(upperLeftSurface);   
    RVector dx(increasingRange(4.0, xBoundary, 20));
    
    for (Index i = 1; i < dx.size()-1; i ++) {
        outerBound.push_back(poly.createNode(upperLeftSurface->pos() 
                            + RVector3(-dx[i], 0)));
    }
    outerBound.push_back(n1);   
    for (Index i = 0; i < 9; i ++) {
        outerBound.push_back(poly.createNode(n1->pos() 
                            + (n2->pos() - n1->pos()) / 10 * (i + 1)));
    }
    outerBound.push_back(n2);   
    for (Index i = 0; i < 9; i ++) {
        outerBound.push_back(poly.createNode(n2->pos() 
                            + (n3->pos() - n2->pos()) / 10 * (i + 1)));
    }
    outerBound.push_back(n3);   
    for (Index i = 0; i < 9; i ++) {
        outerBound.push_back(poly.createNode(n3->pos() 
                            + (n4->pos() - n3->pos()) / 10 * (i + 1)));
    }
    outerBound.push_back(n4);   
    for (Index i = dx.size()-2; i > 0; i --) {
        outerBound.push_back(poly.createNode(upperRightSurface->pos() 
                            + RVector3(dx[i], 0)));
    }
    outerBound.push_back(upperRightSurface);   
    
    for (Index i = 0; i < outerBound.size() -1; i ++){
        poly.createEdge(*outerBound[i], *outerBound[i + 1], -1);
    }
    
    for (Index i = 0; i < innerBound.size() -1; i ++){
        poly.createEdge(*innerBound[i], *innerBound[i + 1], boundMarker);
    }
    
    Mesh boundaryMesh;
    TriangleWrapper(poly, boundaryMesh, "-YpzeAfaq" + str(33));
    
    std::vector < Boundary * > markedBoundaries2(boundaryMesh.findBoundaryByMarker(boundMarker));
    if (markedBoundaries.size() == markedBoundaries2.size()){
        for (std::vector < Cell * >::const_iterator it = boundaryMesh.cells().begin(); 
                it != boundaryMesh.cells().end(); it ++){
            mesh.copyCell(*(*it))->setMarker(cellMarker);
        }
    } else {
        return false;
        mesh.save("inMesh");
        boundaryMesh.save("boundaryMesh");
        throwError(1, WHERE_AM_I + " Sorry!! Boundary mesh is not consistent to input mesh boundary. " 
                        + toStr(markedBoundaries.size()) + " != " + toStr(markedBoundaries2.size()));
    }
    
    if (save){
        mesh.save("shortFOPinMesh");
        boundaryMesh.save("shortFOPboundaryMesh");
    }
    
    //** Fix boundary marker
    mesh.createNeighbourInfos(true);
    b = NULL;
    for (std::vector < Boundary * >::const_iterator it = mesh.boundaries().begin();
        it != mesh.boundaries().end(); it ++){
        b = (*it);
        if (! b->leftCell() || ! b->rightCell()) {
            if (b->center().y() == mesh.ymax()){
                b->setMarker(MARKER_BOUND_HOMOGEN_NEUMANN);
            } else {
                b->setMarker(MARKER_BOUND_MIXED);
            }
        }
    }
    return true;
}

} // namespace GIMLI
