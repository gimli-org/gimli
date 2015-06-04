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

#include "mesh.h"

#include "kdtreeWrapper.h"
#include "memwatch.h"
#include "meshentities.h"
#include "node.h"
#include "shape.h"
#include "sparsematrix.h"
#include "stopwatch.h"

#include <boost/bind.hpp>

#include <map>

namespace GIMLI{

std::ostream & operator << (std::ostream & str, const Mesh & mesh){
    str << "\tNodes: " << mesh.nodeCount() << "\tCells: " << mesh.cellCount() << "\tBoundaries: " << mesh.boundaryCount();
    return str;
}

Mesh::Mesh(Index dim)
    : dimension_(dim),
    rangesKnown_(false), 
    neighboursKnown_(false), 
    tree_(NULL),
    staticGeometry_(true){
        
    oldTet10NumberingStyle_ = true;
    cellToBoundaryInterpolationCache_ = 0;
}

Mesh::Mesh(const std::string & filename)
    : rangesKnown_(false),
    neighboursKnown_(false),
    tree_(NULL),
    staticGeometry_(true){
        
    dimension_ = 3;
    oldTet10NumberingStyle_ = true;
    cellToBoundaryInterpolationCache_ = 0;
    load(filename);
}

Mesh::Mesh(const Mesh & mesh)
    : rangesKnown_(false),
    neighboursKnown_(false),
    tree_(NULL),
    staticGeometry_(true){
        
    oldTet10NumberingStyle_ = true;
    cellToBoundaryInterpolationCache_ = 0;
    copy_(mesh);
}

Mesh & Mesh::operator = (const Mesh & mesh){
    if (this != & mesh){
        copy_(mesh);
    } return *this;
}

void Mesh::copy_(const Mesh & mesh){
    clear();
    rangesKnown_ = false;
    
    setStaticGeometry(mesh.staticGeometry());
    dimension_ = mesh.dim();
    nodeVector_.reserve(mesh.nodeCount());
    
    for (Index i = 0; i < mesh.nodeCount(); i ++){
        this->createNode(mesh.node(i));
    }

    boundaryVector_.reserve(mesh.boundaryCount());
    for (Index i = 0; i < mesh.boundaryCount(); i ++){
        this->createBoundary(mesh.boundary(i));
    }

    cellVector_.reserve(mesh.cellCount());
    for (Index i = 0; i < mesh.cellCount(); i ++){
        this->createCell(mesh.cell(i));
    }

    for (Index i = 0; i < mesh.regionMarker().size(); i ++){
        this->addRegionMarker(mesh.regionMarker()[i]);
    }
    for (Index i = 0; i < mesh.holeMarker().size(); i ++){
        this->addHoleMarker(mesh.holeMarker()[i]);
    }
    
    setExportDataMap(mesh.exportDataMap());
    setCellAttributes(mesh.cellAttributes());
    
    if (mesh.neighboursKnown()){
        this->createNeighbourInfos(true);
    }
//     std::cout << "COPY mesh " << mesh.cell(0) << " " << cell(0) << std::endl;
}

Mesh::~Mesh(){
    clear();
}

void Mesh::setStaticGeometry(bool stat){
    staticGeometry_ = stat;
}

void Mesh::clear(){
    if (tree_) {
        deletePtr()(tree_);
        tree_ = NULL;
    }

    for_each(cellVector_.begin(), cellVector_.end(), deletePtr());
    cellVector_.clear();

    for_each(boundaryVector_.begin(), boundaryVector_.end(), deletePtr());
    boundaryVector_.clear();

    for_each(nodeVector_.begin(), nodeVector_.end(), deletePtr());
    nodeVector_.clear();

    if (cellToBoundaryInterpolationCache_){
        delete cellToBoundaryInterpolationCache_;
    }

    rangesKnown_ = false;
    neighboursKnown_ = false;
}

Node * Mesh::createNode_(const RVector3 & pos, int marker, int id){
    if (id == -1) id = nodeCount();
    nodeVector_.push_back(new Node(pos));
    nodeVector_.back()->setMarker(marker);
    nodeVector_.back()->setId(id);
    return nodeVector_.back();
}

Node * Mesh::createNode(const Node & node){
    return createNode_(node.pos(), node.marker(), -1);
}

Node * Mesh::createNode(double x, double y, double z, int marker){
    return createNode_(RVector3(x, y, z), marker, -1);
}

Node * Mesh::createNode(const RVector3 & pos, int marker){
    return createNode_(pos, marker, -1);
}

Node * Mesh::createNodeWithCheck(const RVector3 & pos, double tol){
    fillKDTree_();

    Node * refNode = tree_->nearest(pos);
    if (refNode){
        if (pos.distance(refNode->pos()) < tol) return refNode;
    }

//     for (Index i = 0, imax = nodeVector_.size(); i < imax; i++){
//         if (pos.distance(nodeVector_[i]->pos()) < 1e-6) return nodeVector_[i];
//     }

    Node * newNode = createNode(pos);
    tree_->insert(newNode);
    return newNode;
}

Boundary * Mesh::createBoundary(std::vector < Node * > & nodes, int marker){
    switch (nodes.size()){
      case 1: return createBoundaryChecked_< NodeBoundary >(nodes, marker); break;
      case 2: return createBoundaryChecked_< Edge >(nodes, marker); break;
      case 3: {
        if (dimension_ == 2)
            return createBoundaryChecked_< Edge3 >(nodes, marker);
        return createBoundaryChecked_< TriangleFace >(nodes, marker); } break;
      case 4: return createBoundaryChecked_< QuadrangleFace >(nodes, marker); break;
      case 6: return createBoundaryChecked_< Triangle6Face >(nodes, marker); break;
      case 8: return createBoundaryChecked_< Quadrangle8Face >(nodes, marker); break;
    }
    std::cout << WHERE_AM_I << "WHERE_AM_I << cannot determine boundary for nodes: " << nodes.size() << std::endl;
    return NULL;
}

Boundary * Mesh::createBoundary(const Boundary & bound){
    std::vector < Node * > nodes(bound.nodeCount());
    for (Index i = 0; i < bound.nodeCount(); i ++) nodes[i] = &node(bound.node(i).id());
    return createBoundary(nodes, bound.marker());
}

Boundary * Mesh::createBoundary(const Cell & cell){
    std::vector < Node * > nodes(cell.nodeCount());
    for (Index i = 0; i < cell.nodeCount(); i ++) nodes[i] = &node(cell.node(i).id());
    return createBoundary(nodes, cell.marker());
}

Boundary * Mesh::createNodeBoundary(Node & n1, int marker){
    std::vector < Node * > nodes(1); nodes[0] = & n1;
    return createBoundaryChecked_< NodeBoundary >(nodes, marker);
}

Boundary * Mesh::createEdge(Node & n1, Node & n2, int marker){
    std::vector < Node * > nodes(2);  nodes[0] = & n1; nodes[1] = & n2;
    return createBoundaryChecked_< Edge >(nodes, marker);
}

Boundary * Mesh::createEdge3(Node & n1, Node & n2, Node & n3, int marker){
    std::vector < Node * > nodes(3); nodes[0] = & n1; nodes[1] = & n2; nodes[2] = & n3;
    return createBoundaryChecked_< Edge3 >(nodes, marker);
}

Boundary * Mesh::createTriangleFace(Node & n1, Node & n2, Node & n3, int marker){
    std::vector < Node * > nodes(3);  nodes[0] = & n1; nodes[1] = & n2; nodes[2] = & n3;
    return createBoundaryChecked_< TriangleFace >(nodes, marker);
}

Boundary * Mesh::createQuadrangleFace(Node & n1, Node & n2, Node & n3, Node & n4, int marker){
    std::vector < Node * > nodes(4);  nodes[0] = & n1; nodes[1] = & n2; nodes[2] = & n3, nodes[3] = & n4;
    return createBoundaryChecked_< QuadrangleFace >(nodes, marker);
}

Cell * Mesh::createCell(int marker){
    std::vector < Node * > nodes(0);
    return createCell_< Cell >(nodes, marker, cellCount());
}

Cell * Mesh::createCell(std::vector < Node * > & nodes, int marker){
    switch (nodes.size()){
        case 0: return createCell_< Cell >(nodes, marker, cellCount()); break;
        case 2: return createCell_< EdgeCell >(nodes, marker, cellCount()); break;
        case 3:
            switch (dimension_){
                case 1: return createCell_< Edge3Cell >(nodes, marker, cellCount()); break;
                case 2: return createCell_< Triangle >(nodes, marker, cellCount()); break;
            }
            break;
        case 4:
            switch (dimension_){
                case 2: return createCell_< Quadrangle >(nodes, marker, cellCount()); break;
                case 3: return createCell_< Tetrahedron >(nodes, marker, cellCount()); break;
            }
            break;
        case 5: return createCell_< Pyramid >(nodes, marker, cellCount()); break;
        case 6:
            switch (dimension_){
                case 2: return createCell_< Triangle6 >(nodes, marker, cellCount()); break;
                case 3: return createCell_< TriPrism >(nodes, marker, cellCount()); break;
            }
            break;
        case 8:
            switch (dimension_){
                case 2: return createCell_< Quadrangle8 >(nodes, marker, cellCount()); break;
                case 3: return createCell_< Hexahedron >(nodes, marker, cellCount()); break;
            }
            break;
        case 10: return createCell_< Tetrahedron10 >(nodes, marker, cellCount()); break;
        case 13: return createCell_< Pyramid13 >(nodes, marker, cellCount()); break;
        case 15: return createCell_< TriPrism15 >(nodes, marker, cellCount()); break;
        case 20: return createCell_< Hexahedron20 >(nodes, marker, cellCount()); break;

    }
    std::cout << WHERE_AM_I << "WHERE_AM_I << cannot determine cell for nodes: " << nodes.size() << " for dim: " << dimension_ << std::endl;
    return NULL;
}

Cell * Mesh::createCell(const Cell & cell){
    std::vector < Node * > nodes(cell.nodeCount());
    for (Index i = 0; i < cell.nodeCount(); i ++) nodes[i] = &node(cell.node(i).id());
    return createCell(nodes, cell.marker());
}

Cell * Mesh::createTriangle(Node & n1, Node & n2, Node & n3, int marker){
    std::vector < Node * > nodes(3);  nodes[0] = & n1; nodes[1] = & n2; nodes[2] = & n3;
    return createCell_< Triangle >(nodes, marker, cellCount());
}

Cell * Mesh::createQuadrangle(Node & n1, Node & n2, Node & n3, Node & n4, int marker){
    std::vector < Node * > nodes(4);
    nodes[0] = & n1; nodes[1] = & n2; nodes[2] = & n3; nodes[3] = & n4;
    return createCell_< Quadrangle >(nodes, marker, cellCount());
}

Cell * Mesh::createTetrahedron(Node & n1, Node & n2, Node & n3, Node & n4, int marker){
    std::vector < Node * > nodes(4);
    nodes[0] = & n1; nodes[1] = & n2; nodes[2] = & n3; nodes[3] = & n4;
    return createCell_< Tetrahedron >(nodes, marker, cellCount());
}

Cell * Mesh::copyCell(const Cell & cell){
    std::vector < Node * > nodes(cell.nodeCount());
    for (Index i = 0; i < nodes.size(); i ++) {
        nodes[i] = createNodeWithCheck(cell.node(i).pos());
        nodes[i]->setMarker(cell.node(i).marker());
    }
    Cell * c = createCell(nodes);

    c->setMarker(cell.marker());
    c->setAttribute(cell.attribute());
    return c;
}

Boundary * Mesh::copyBoundary(const Boundary & bound){
    std::vector < Node * > nodes(bound.nodeCount());
    for (Index i = 0; i < nodes.size(); i ++) {
        nodes[i] = createNodeWithCheck(bound.node(i).pos());
        nodes[i]->setMarker(bound.node(i).marker());
    }
    Boundary * b = createBoundary(nodes);

    b->setMarker(bound.marker());
    return b;
}

Node & Mesh::node(Index i) {
    if (i > nodeCount() - 1){
        std::cerr << WHERE_AM_I << " requested node: " << i << " does not exist." << std::endl;
        exit(EXIT_MESH_NO_NODE);
    } return *nodeVector_[i];
}

Node & Mesh::node(Index i) const {
    if (i > nodeCount() - 1){
        std::cerr << WHERE_AM_I << " requested node: " << i << " does not exist." << std::endl;
        exit(EXIT_MESH_NO_NODE);
    } return *nodeVector_[i];
}

Cell & Mesh::cell(Index i) const {
    if (i > cellCount() - 1){
      std::cerr << WHERE_AM_I << " requested cell: " << i << " does not exist." << std::endl;
      exit(EXIT_MESH_NO_CELL);
    } return *cellVector_[i];
}

Cell & Mesh::cell(Index i) {
    if (i > cellCount() - 1){
      std::cerr << WHERE_AM_I << " requested cell: " << i << " does not exist." << std::endl;
      exit(EXIT_MESH_NO_CELL);
    } return *cellVector_[i];
}

Boundary & Mesh::boundary(Index i) const {
    if (i > boundaryCount() - 1){
      std::cerr << WHERE_AM_I << " requested boundary: " << i << " does not exist." << std::endl;
      exit(EXIT_MESH_NO_BOUNDARY);
    } return *boundaryVector_[i];
}

Boundary & Mesh::boundary(Index i) {
    if (i > boundaryCount() - 1){
      std::cerr << WHERE_AM_I << " requested boundary: " << i << " does not exist." << std::endl;
      exit(EXIT_MESH_NO_BOUNDARY);
    } return *boundaryVector_[i];
}

void Mesh::findRange_() const{
    if (!rangesKnown_ || !staticGeometry_){
        minRange_ = RVector3(MAX_DOUBLE, MAX_DOUBLE, MAX_DOUBLE);
        maxRange_ = -minRange_;
        for (Index i=0; i < nodeVector_.size(); i ++){
            for (Index j=0; j < 3; j ++){
                minRange_[j] = min(nodeVector_[i]->pos()[j], minRange_[j]);
                maxRange_[j] = max(nodeVector_[i]->pos()[j], maxRange_[j]);
            }
        }
        rangesKnown_ = true;
    }
}

void Mesh::createHull(const Mesh & mesh){
    if (this->dim() == 3 && mesh.dim() == 2){
        clear();
        rangesKnown_ = false;
        nodeVector_.reserve(mesh.nodeCount());
        for (Index i = 0; i < mesh.nodeCount(); i ++) createNode(mesh.node(i));

        boundaryVector_.reserve(mesh.cellCount());
        for (Index i = 0; i < mesh.cellCount(); i ++) createBoundary(mesh.cell(i));
    } else {
        std::cerr << WHERE_AM_I << " increasing dimension fails, you should set the dimension for this mesh to 3" << std::endl;
    }
}

Index Mesh::findNearestNode(const RVector3 & pos){
    fillKDTree_();
    return tree_->nearest(pos)->id();
}

IndexArray cellIDX__;

Cell * Mesh::findCellBySlopeSearch_(const RVector3 & pos, Cell * start, 
                                    size_t & count, bool useTagging) const {
    Cell * cell = start;
    
    Index cellCounter = 0; //** for avoiding infinite loop
    do {
        if (cell->tagged() && useTagging) {
            cell = NULL;
        } else {
            cell->tag();
            cellIDX__.push_back(cell->id());
            RVector sf;

//             std::cout << cellIDX__.size() << " testpos: " << pos << std::endl;
//             std::cout << "cell: " << *cell << " touch: " << cell->shape().isInside(pos, true) << std::endl;
//             for (Index i = 0; i < cell->nodeCount() ; i ++){
//                 std::cout << cell->node(i)<< std::endl;
//             }

            if (cell->shape().isInside(pos, sf, false)) {
                return cell;
            } else {
                
                if (!neighboursKnown_){
                    const_cast<Mesh*>(this)->createNeighbourInfosCell_(cell);
//                     for (Index j = 0; j < cell->neighbourCellCount(); j++){
//                          cell->findNeighbourCell(j);
//                     }
                }
                
//                 for (Index i = 0; i < cell->neighbourCellCount(); i ++ ){
//                     if (cell->neighbourCell(i)){
//                         std::cout << "\t " << i << " " << *cell->neighbourCell(i) << std::endl;
//                     } else {
//                         std::cout << "\t " << i << " " << 0 << std::endl;
//                     }
//                 }
                
                cell = cell->neighbourCell(sf);

//                 std::cout << "sf: " << sf << std::endl;
//                 std::cout << "neighCell " << cell << std::endl;
            }
            count++;
            if (count == 1000){
//                 std::cout << "testpos: " << pos << std::endl;
//                 std::cout << "cell: " << this->cell(cellIDX__.back()) << std::endl;
//                 for (Index i = 0;
//                      i < this->cell(cellIDX__.back()).nodeCount(); i ++){
//                     std::cout << this->cell(cellIDX__.back()).node(i)<< std::endl;
//                 }

                std::cout << WHERE_AM_I << " exit with submesh " << cellIDX__.size() << std::endl;
                std::cout << "probably cant find a cell for " << pos << std::endl;
                
                Mesh subMesh; subMesh.createMeshByCellIdx(*this, cellIDX__);

                subMesh.exportVTK("submesh");
                this->exportVTK("submeshParent");
                return NULL;
                exit(0);
            }
        }
    } while (++cellCounter < cellCount() && cell);

    return NULL;
}

Cell * Mesh::findCell(const RVector3 & pos, size_t & count,
                      bool extensive) const {
    bool bruteForce = false;
    Cell * cell = NULL;
            
    if (bruteForce){
        for (Index i = 0; i < this->cellCount(); i ++) {
            count++;
            if (cellVector_[i]->shape().isInside(pos)){
                cell = cellVector_[i];
                
//                 std::cout << "testpos: " << pos << std::endl;
//                 std::cout << "cell: " << *cell<< std::endl;
                break;
            }
        }
    } else {
        Stopwatch swatch(true);
        cellIDX__.clear();
        count = 0;
        fillKDTree_();
        Node * refNode = tree_->nearest(pos);
        
        if (!refNode){
            std::cout << "pos: " << pos << std::endl;
            throwError(1, WHERE_AM_I + 
                       " no nearest node to pos. This is a empty mesh");
        }
        if (refNode->cellSet().empty()){
            std::cout << "Node: " << *refNode << std::endl;
            throwError(1, WHERE_AM_I + 
                       " no cells for this node. This is a corrupt mesh");
        }
            //std::cout << "Node: " << *refNode << std::endl;
//         for (std::set< Cell * >::iterator it = refNode->cellSet().begin();
//              it != refNode->cellSet().end(); it ++){
//             std::cout << (*it)->id() << std::endl;
//         }
        cell = findCellBySlopeSearch_(pos, *refNode->cellSet().begin(),
                                      count, false);
        if (cell) return cell;

//         exportVTK("slopesearch");
//         exit(0);
        if (extensive || 0){
//             __M
//             std::cout << "More expensive test here" << std::endl;
            cellIDX__.clear();
            std::for_each(cellVector_.begin(), cellVector_.end(), std::mem_fun(&Cell::untag));
            //!** *sigh, no luck with simple kd-tree search, try more expensive full slope search
            count = 0;
            for (Index i = 0; i < this->cellCount(); i ++) {
                cell = cellVector_[i];
                cell = findCellBySlopeSearch_(pos, cell, count, true);
                if (cell) {
                    
                    break;
                }
            }
        } else {
            return NULL;
        }
  //      std::cout << " n: " << count;
    }
    return cell;
}

std::vector < Boundary * > Mesh::findBoundaryByMarker(int marker) const {
    return findBoundaryByMarker(marker, marker + 1);
}

std::vector < Boundary * > Mesh::findBoundaryByMarker(int from, int to) const {
    std::vector < Boundary * > vBounds;
    vBounds.reserve(boundaryCount());

    for(std::vector< Boundary * >::const_iterator it = boundaryVector_.begin();
                                               it != boundaryVector_.end(); it++){
        if ((*it)->marker() >= from && (*it)->marker() < to) vBounds.push_back((*it));
    }

    return vBounds;
}

void Mesh::setBoundaryMarker(const IndexArray & ids, int marker){
    for(IndexArray::const_iterator it = ids.begin(); it != ids.end(); it++){
        if (*it < boundaryCount()){
            boundaryVector_[*it]->setMarker(marker);
        }
    }
}
    

std::vector < Cell * > Mesh::findCellByMarker(int from, int to) const {
    if (to == -1) to = MAX_INT;
    else if (to == 0) to = from + 1;

    std::vector < Cell * > vCell;
    vCell.reserve(cellCount());
    for(std::vector< Cell * >::const_iterator it = cellVector_.begin();
                                               it != cellVector_.end(); it++){
        if ((*it)->marker() >= from && (*it)->marker() < to) vCell.push_back((*it));
    }
    return vCell;
}

std::vector < Cell * > Mesh::findCellByAttribute(double from, double to) const {
    std::vector < Cell * > vCell;
    vCell.reserve(cellCount());

    if (to < TOLERANCE){
        for (Index i = 0; i < cellCount(); i ++){
            if ((cell(i).attribute() - from) < TOLERANCE) vCell.push_back(cellVector_[i]);
        }
    } else {
        if (to == -1) to = MAX_DOUBLE;
        for (Index i = 0; i < cellCount(); i ++){
            if (cell(i).attribute() >= from && cell(i).attribute() < to)
                vCell.push_back(cellVector_[i]);
        }
    }
    return vCell;
}

std::vector< Cell * > Mesh::cells(const IndexArray & ids) const {
    std::vector < Cell * > v(ids.size());
    for (Index i = 0; i < ids.size(); i ++) v[i] = cellVector_[ids[i]];
    return v;
}

std::vector< Node * > Mesh::nodes(const IndexArray & ids) const{
    std::vector < Node * > v(ids.size());
    for (Index i = 0; i < ids.size(); i ++) v[i] = nodeVector_[ids[i]];
    return v;
}
    
std::vector< Boundary * > Mesh::boundaries(const IndexArray & ids) const{
    std::vector < Boundary * > v(ids.size());
    for (Index i = 0; i < ids.size(); i ++) v[i] = boundaryVector_[ids[i]];
    return v;
}
    
void Mesh::setCellMarker(const IndexArray & ids, int marker){
    for(IndexArray::const_iterator it = ids.begin(); it != ids.end(); it++){
        if (*it < cellCount()){
            cellVector_[*it]->setMarker(marker);
        }
    }
}

void Mesh::setCellMarker(const RVector & attribute){
    if (attribute.size() >= cellVector_.size()){
        for (Index i = 0; i < cellVector_.size(); i ++){
            cellVector_[i]->setMarker(int(attribute[i]));
        }
    } else {
        throwError(1,"Mesh::setCellMarker: attribute size to small: " + 
            str(attribute.size()) + " < " + str(cellCount()));
    }
}
    
IndexArray Mesh::findNodesIdxByMarker(int marker) const {
    return find(this->nodeMarker() == marker);
}

// std::list < Index > Mesh::findListNodesIdxByMarker(int marker) const {
//     std::list < Index > idx;
//     for (Index i = 0; i < nodeCount(); i ++) {
//         if (node(i).marker() == marker) idx.push_back(i);
//     }
//     return idx;
// }

std::vector < RVector3 > Mesh::positions() const {
    IndexArray idx(this->nodeCount());
    std::generate(idx.begin(), idx.end(), IncrementSequence< Index >(0));
    return this->positions(idx);
}

std::vector < RVector3 > Mesh::positions(const IndexArray & idx) const {
    std::vector < RVector3 > pos; pos.reserve(idx.size());
    for (Index i = 0; i < idx.size(); i ++) {
        pos.push_back(node(idx[i]).pos());
    }
    return pos;
}

R3Vector Mesh::nodeCenters() const {
    R3Vector pos(this->cellCount());
    std::transform(cellVector_.begin(), cellVector_.end(), pos.begin(),
                   std::mem_fun(&Cell::center));
    return pos;
}

std::vector < RVector3 > Mesh::cellCenters() const {
    std::vector < RVector3 > pos(this->cellCount());
    std::transform(cellVector_.begin(), cellVector_.end(), pos.begin(),
                   std::mem_fun(&Cell::center));
    return pos;
}

std::vector < RVector3 > Mesh::boundaryCenters() const {
    std::vector < RVector3 > pos(this->boundaryCount());
    std::transform(boundaryVector_.begin(), boundaryVector_.end(), pos.begin(),
                   std::mem_fun(&Boundary::center));
    return pos;
}

RVector & Mesh::cellSizes() const{
    
    if (cellSizesCache_.size() != cellCount()){
        cellSizesCache_.resize(cellCount());
        
        std::transform(cellVector_.begin(), cellVector_.end(),
                       cellSizesCache_.begin(),
                       std::mem_fun(&Cell::size));
    } else {
        if (!staticGeometry_){
            cellSizesCache_.resize(0);
            return this->cellSizes();
        }
    }
    
    return cellSizesCache_;
}

RVector & Mesh::boundarySizes() const{
    if (boundarySizesCache_.size() != boundaryCount()){
        boundarySizesCache_.resize(boundaryCount());
        
        std::transform(boundaryVector_.begin(), boundaryVector_.end(), 
                       boundarySizesCache_.begin(),
                       std::mem_fun(&MeshEntity::size));
    } else {
        if (!staticGeometry_){
            boundarySizesCache_.resize(0);
            return this->boundarySizes();
        }
    }
    
    return boundarySizesCache_;
}

std::vector< RVector3 > & Mesh::boundarySizedNormals() const {
    if (boundarySizedNormCache_.size() != boundaryCount()){
        boundarySizedNormCache_.resize(boundaryCount());
        
        for (Index i = 0; i < boundaryCount(); i ++){
            boundarySizedNormCache_[i] = boundaryVector_[i]->norm() * boundaryVector_[i]->size();
        }
        
    } else {
        if (!staticGeometry_){
            boundarySizedNormCache_.resize(0);
            return this->boundarySizedNormals();
        }
    }
    
    return boundarySizedNormCache_;
}

void Mesh::sortNodes(const IndexArray & perm){

    for (Index i = 0; i < nodeVector_.size(); i ++) nodeVector_[i]->setId(perm[i]);
  //    sort(nodeVector_.begin(), nodeVector_.end(), std::less< int >(mem_fun(&BaseEntity::id)));
    sort(nodeVector_.begin(), nodeVector_.end(), lesserId< Node >);
    recountNodes();
}

void Mesh::recountNodes(){
    for (Index i = 0; i < nodeVector_.size(); i ++) nodeVector_[i]->setId(i);
}

void Mesh::createClosedGeometry(const std::vector < RVector3 > & vPos, int nSegments, double dxInner){
    THROW_TO_IMPL
    //this function should not be part of mesh
    //EAMeshWrapper eamesh(vPos, nSegments, dxInner, *this);
}

void Mesh::createClosedGeometryParaMesh(const std::vector < RVector3 > & vPos, int nSegments, double dxInner){
    createClosedGeometry(vPos, nSegments, dxInner);
    createNeighbourInfos();
    for (Index i = 0; i < cellCount(); i ++) cell(i).setMarker(i);
}

void Mesh::createClosedGeometryParaMesh(const std::vector < RVector3 > & vPos, int nSegments,
                                         double dxInner, const std::vector < RVector3 > & addit){
    THROW_TO_IMPL
    //this function should not be part of meshEntities
//   EAMeshWrapper eamesh;
//   eamesh.createMesh(vPos, nSegments, dxInner, addit);
//   eamesh.mesh(*this);
//   createNeighbourInfos();
//   for (Index i = 0; i < cellCount(); i ++) cell(i).setMarker(i);
}

Mesh Mesh::createH2() const {
    Mesh ret(this->dimension());
    ret.createRefined_(*this, false, true);
    ret.setCellAttributes(ret.cellMarker());
    return ret;
}

Mesh Mesh::createP2() const {
    Mesh ret(this->dimension());
    ret.createRefined_(*this, true, false);
    return ret;
}

int markerT(Node * n0, Node * n1){
    if (n0->marker() == -99 && n1->marker() == -99) return -1;
    if (n0->marker() == -99) return n1->marker();
    if (n1->marker() == -99) return n0->marker();
    if (n0->marker() == n1->marker()) return n1->marker();
    else return 0;
}

Node * Mesh::createRefinementNode_(Node * n0, Node * n1, std::map< std::pair < Index, Index >, Node * > & nodeMatrix){
    Node * n = nodeMatrix[std::make_pair(n0->id(), n1->id())];

    if (!n){
        if (n0 == n1){
            n = this->createNode(n0->pos(), n0->marker()) ;
            nodeMatrix[std::make_pair(n0->id(), n0->id())] = n;
        } else {
            n = this->createNode((n0->pos() + n1->pos()) / 2.0, markerT(n0, n1));
            nodeMatrix[std::make_pair(n0->id(), n1->id())] = n;
            nodeMatrix[std::make_pair(n1->id(), n0->id())] = n;
        }
    }
    return n;
}

void Mesh::createRefined_(const Mesh & mesh, bool p2, bool h2){

    this->clear();

    std::map< std::pair < Index, Index >, Node * > nodeMatrix;

    for (Index i = 0, imax = mesh.nodeCount(); i < imax; i ++) {
        this->createRefinementNode_(&mesh.node(i), &mesh.node(i), nodeMatrix);
    }

    std::vector < Node * > n;
    for (Index i = 0, imax = mesh.cellCount(); i < imax; i ++){
        const Cell & c = mesh.cell(i);
        Index cID = i;

        switch (c.rtti()){
            case MESH_EDGE_CELL_RTTI:
                n.resize(3);
                n[0] = &node(mesh.cell(i).node(0).id());
                n[1] = &node(mesh.cell(i).node(1).id());
                n[2] = createRefinementNode_(n[0], n[1], nodeMatrix);

                if (h2){
                    std::vector < Node * > e1(2);
                    e1[0] = n[0]; e1[1] = n[2];
                    this->createCell(e1, cID);
                    e1[0] = n[2]; e1[1] = n[1];
                    this->createCell(e1, cID);
                }
                break;
            case MESH_TRIANGLE_RTTI:
                n.resize(6);
                n[0] = &node(mesh.cell(i).node(0).id());
                n[1] = &node(mesh.cell(i).node(1).id());
                n[2] = &node(mesh.cell(i).node(2).id());

                n[3] = createRefinementNode_(n[0], n[1], nodeMatrix);
                n[4] = createRefinementNode_(n[1], n[2], nodeMatrix);
                n[5] = createRefinementNode_(n[2], n[0], nodeMatrix);

                if (h2){
                    this->createTriangle(*n[0], *n[3], *n[5], cID);
                    this->createTriangle(*n[1], *n[4], *n[3], cID);
                    this->createTriangle(*n[2], *n[5], *n[4], cID);
                    this->createTriangle(*n[3], *n[4], *n[5], cID);
                }

                break;
            case MESH_QUADRANGLE_RTTI:
                n.resize(8);
                n[0] = &node(mesh.cell(i).node(0).id());
                n[1] = &node(mesh.cell(i).node(1).id());
                n[2] = &node(mesh.cell(i).node(2).id());
                n[3] = &node(mesh.cell(i).node(3).id());

                n[4] = createRefinementNode_(n[0], n[1], nodeMatrix);
                n[5] = createRefinementNode_(n[1], n[2], nodeMatrix);
                n[6] = createRefinementNode_(n[2], n[3], nodeMatrix);
                n[7] = createRefinementNode_(n[3], n[0], nodeMatrix);

                if (h2){
                    Node *n8 = this->createNode(c.shape().xyz(RVector3(0.5, 0.5)));
                    this->createQuadrangle(*n[0], *n[4], *n8, *n[7], cID);
                    this->createQuadrangle(*n[1], *n[5], *n8, *n[4], cID);
                    this->createQuadrangle(*n[2], *n[6], *n8, *n[5], cID);
                    this->createQuadrangle(*n[3], *n[7], *n8, *n[6], cID);
                }

                break;
            case MESH_TETRAHEDRON_RTTI:
                n.resize(10);
                if (oldTet10NumberingStyle_){

                    for (Index j = 0; j < n.size(); j ++) {

                        n[j] = createRefinementNode_(& this->node(c.node(Tet10NodeSplitZienk[j][0]).id()),
                                                        & this->node(c.node(Tet10NodeSplitZienk[j][1]).id()),
                                                        nodeMatrix);
                    }

                    if (h2){
                        this->createTetrahedron(*n[4], *n[6], *n[5], *n[0], cID);
                        this->createTetrahedron(*n[4], *n[5], *n[6], *n[9], cID);
                        this->createTetrahedron(*n[7], *n[9], *n[4], *n[1], cID);
                        this->createTetrahedron(*n[7], *n[4], *n[9], *n[5], cID);
                        this->createTetrahedron(*n[8], *n[7], *n[5], *n[2], cID);
                        this->createTetrahedron(*n[8], *n[5], *n[7], *n[9], cID);
                        this->createTetrahedron(*n[6], *n[9], *n[8], *n[3], cID);
                        this->createTetrahedron(*n[6], *n[8], *n[9], *n[5], cID);
                    }

                } else {
                    for (Index j = 0; j < n.size(); j ++) {
                        n[j] = createRefinementNode_(& this->node(c.node(Tet10NodeSplit[j][0]).id()),
                                                        & this->node(c.node(Tet10NodeSplit[j][1]).id()),
                                                        nodeMatrix);
                    }

                    if (h2){
                        THROW_TO_IMPL
                    }
                }

                break;
            case MESH_HEXAHEDRON_RTTI:
                n.resize(20);
                for (Index j = 0; j < n.size(); j ++) {
                    n[j] = createRefinementNode_(& this->node(c.node(Hex20NodeSplit[j][0]).id()),
                                                 & this->node(c.node(Hex20NodeSplit[j][1]).id()),
                                                 nodeMatrix);
                }
                if (h2){
/* 27 new nodes 3 x 9 = 8 nodes + 12 edges + 6 facets + 1 center                     
          7-----14------6  \n
         /|            /|  \n
        / |           / |  \n
      15 19  -21-    13 18 \n
      /   |     24  /   |  \n
     /    |        /    |  \n
    4-----12------5  23 |  \n
    | 25  3-----10|-----2  \n
    |    /        |    /   \n
   16   / -22-    17  /    \n
    | 11    -20-  |  9     \n
    | /           | /      \n
    |/            |/       \n
    0------8------1        \n

*/
                    Node *n20 = createRefinementNode_(n[8], n[10], nodeMatrix);
                    Node *n21 = createRefinementNode_(n[12], n[14], nodeMatrix);
                    Node *n22 = createRefinementNode_(n[8], n[12], nodeMatrix);
                    Node *n23 = createRefinementNode_(n[9], n[13], nodeMatrix);
                    Node *n24 = createRefinementNode_(n[10], n[14], nodeMatrix);
                    Node *n25 = createRefinementNode_(n[11], n[15], nodeMatrix);

                    Node *n26 = createRefinementNode_(n20, n21, nodeMatrix);

                    std::vector < Node* > ns(8);
                    Node *n1_[]={ n[0], n[8], n20, n[11], n[16], n22, n26, n25 }; std::copy(&n1_[0], &n1_[8], &ns[0]);
                    this->createCell(ns, cID);

                    Node *n2_[]= { n[8], n[1], n[9], n20, n22, n[17], n23, n26 };  std::copy(&n2_[0], &n2_[8], &ns[0]);
                    this->createCell(ns, cID);

                    Node *n3_[]= { n[11], n20, n[10], n[3], n25, n26, n24, n[19] };  std::copy(&n3_[0], &n3_[8], &ns[0]);
                    this->createCell(ns, cID);

                    Node *n4_[]= { n20, n[9], n[2], n[10], n26, n23, n[18], n24 };  std::copy(&n4_[0], &n4_[8], &ns[0]);
                    this->createCell(ns, cID);

                    Node *n5_[]= { n[16], n22, n26, n25, n[4], n[12], n21, n[15] };  std::copy(&n5_[0], &n5_[8], &ns[0]);
                    this->createCell(ns, cID);

                    Node *n6_[]= { n22, n[17], n23, n26, n[12], n[5], n[13], n21 };  std::copy(&n6_[0], &n6_[8], &ns[0]);
                    this->createCell(ns, cID);

                    Node *n7_[]= { n25, n26, n24, n[19], n[15], n21, n[14], n[7] };  std::copy(&n7_[0], &n7_[8], &ns[0]);
                    this->createCell(ns, cID);

                    Node *n8_[]= { n26, n23, n[18], n24, n21, n[13], n[6], n[14] };  std::copy(&n8_[0], &n8_[8], &ns[0]);
                    this->createCell(ns, cID);

                }

                break;
            case MESH_TRIPRISM_RTTI:
                n.resize(15);
                for (Index j = 0; j < n.size(); j ++) {
                    n[j] = createRefinementNode_(& this->node(c.node(Prism15NodeSplit[j][0]).id()),
                                                 & this->node(c.node(Prism15NodeSplit[j][1]).id()),
                                                 nodeMatrix);
                }
                if (h2){

                    Node *nf1 = createRefinementNode_(n[6], n[9], nodeMatrix);
                    Node *nf2 = createRefinementNode_(n[7], n[10], nodeMatrix);
                    Node *nf3 = createRefinementNode_(n[8], n[11], nodeMatrix);

                    std::vector < Node* > ns(6);
                    Node *n1_[]={ n[0], n[6], n[8], n[12], nf1, nf3 }; std::copy(&n1_[0], &n1_[6], &ns[0]);
                    this->createCell(ns, cID);
                    Node *n2_[]={ n[1], n[7], n[6], n[13], nf2, nf1 }; std::copy(&n2_[0], &n2_[6], &ns[0]);
                    this->createCell(ns, cID);
                    Node *n3_[]={ n[2], n[8], n[7], n[14], nf3, nf2 }; std::copy(&n3_[0], &n3_[6], &ns[0]);
                    this->createCell(ns, cID);
                    Node *n4_[]={ n[6], n[7], n[8], nf1, nf2, nf3 }; std::copy(&n4_[0], &n4_[6], &ns[0]);
                    this->createCell(ns, cID);

                    Node *n5_[]={ n[12], nf1, nf3, n[3], n[9], n[11] }; std::copy(&n5_[0], &n5_[6], &ns[0]);
                    this->createCell(ns, cID);
                    Node *n6_[]={ n[13], nf2, nf1, n[4], n[10], n[9] }; std::copy(&n6_[0], &n6_[6], &ns[0]);
                    this->createCell(ns, cID);
                    Node *n7_[]={ n[14], nf3, nf2, n[5], n[11], n[10] }; std::copy(&n7_[0], &n7_[6], &ns[0]);
                    this->createCell(ns, cID);
                    Node *n8_[]={ nf1, nf2, nf3, n[9], n[10], n[11] }; std::copy(&n8_[0], &n8_[6], &ns[0]);
                    this->createCell(ns, cID);

                }
                break;
            case MESH_PYRAMID_RTTI:
                n.resize(13);
                for (Index j = 0; j < n.size(); j ++) {
                    n[j] = createRefinementNode_(& this->node(c.node(Pyramid13NodeSplit[j][0]).id()),
                                                    & this->node(c.node(Pyramid13NodeSplit[j][1]).id()),
                                                    nodeMatrix);
                }
                if (h2){
                    THROW_TO_IMPL
                }
                break;
            default: std::cerr << c.rtti() <<" " << std::endl; THROW_TO_IMPL  break;
        }

        if (p2 && !h2){
            createCell(n, cID);
        }

    } // for_each cell

    for (Index i = 0; i < mesh.boundaryCount(); i++){

        const Boundary & b = mesh.boundary(i);

        switch (b.rtti()){
            case MESH_BOUNDARY_NODE_RTTI:
                n.resize(1);
                n[0] = &node(b.node(0).id());
                break;
            case MESH_EDGE_RTTI:
                n.resize(3);
                n[0] = &node(b.node(0).id());
                n[1] = &node(b.node(1).id());
                n[2] = createRefinementNode_(n[0], n[1], nodeMatrix);
                break;
            case MESH_TRIANGLEFACE_RTTI:
                n.resize(6);
                n[0] = &node(b.node(0).id());
                n[1] = &node(b.node(1).id());
                n[2] = &node(b.node(2).id());

                n[3] = createRefinementNode_(n[0], n[1], nodeMatrix);
                n[4] = createRefinementNode_(n[1], n[2], nodeMatrix);
                n[5] = createRefinementNode_(n[2], n[0], nodeMatrix);
                break;
            case MESH_QUADRANGLEFACE_RTTI:
                n.resize(8);
                n[0] = &node(b.node(0).id());
                n[1] = &node(b.node(1).id());
                n[2] = &node(b.node(2).id());
                n[3] = &node(b.node(3).id());

                n[4] = createRefinementNode_(n[0], n[1], nodeMatrix);
                n[5] = createRefinementNode_(n[1], n[2], nodeMatrix);
                n[6] = createRefinementNode_(n[2], n[3], nodeMatrix);
                n[7] = createRefinementNode_(n[3], n[0], nodeMatrix);
                break;
            default: std::cerr << b.rtti() <<" " << std::endl; THROW_TO_IMPL  break;
        }

        if (p2 && !h2){
            createBoundary(n, b.marker());
        } else {
            switch (b.rtti()){
                case MESH_BOUNDARY_NODE_RTTI:
                    this->createBoundary(n, b.marker());
                    break;
                case MESH_EDGE_RTTI:
                    this->createEdge(*n[0], *n[2], b.marker());
                    this->createEdge(*n[2], *n[1], b.marker());
                    break;
                case MESH_TRIANGLEFACE_RTTI:
                    this->createTriangleFace(*n[0], *n[3], *n[5], b.marker());
                    this->createTriangleFace(*n[1], *n[4], *n[3], b.marker());
                    this->createTriangleFace(*n[2], *n[5], *n[4], b.marker());
                    this->createTriangleFace(*n[3], *n[4], *n[5], b.marker());
                    break;
                case MESH_QUADRANGLEFACE_RTTI:
                    /* 
                     * 3---6---2
                     * |   |   |            
                     * 7---8---5
                     * |   |   |
                     * 0---4---1
                    */
                    Node *n8 = nodeMatrix[std::make_pair(n[4]->id(), n[6]->id())];
                    if (!n8) n8 = createRefinementNode_(n[5], n[7], nodeMatrix);

                    this->createQuadrangleFace(*n[0], *n[4], *n8, *n[7], b.marker());
                    this->createQuadrangleFace(*n[1], *n[5], *n8, *n[4], b.marker());
                    this->createQuadrangleFace(*n[2], *n[6], *n8, *n[5], b.marker());
                    this->createQuadrangleFace(*n[3], *n[7], *n8, *n[6], b.marker());

                    break;

            }
        } // if not p2
    } // for_each boundary
    
    //! Copy data if available
    for (std::map < std::string, RVector >::const_iterator 
         it = mesh.dataMap().begin(); it != mesh.dataMap().end(); it ++){
        
        if (it->second.size() == mesh.cellCount()){
            addData(it->first, it->second(this->cellMarker()));
        } else if (it->second.size() == mesh.nodeCount()){
            THROW_TO_IMPL
        } else if (it->second.size() == mesh.boundaryCount()){
            THROW_TO_IMPL
        }
    }
    //! copy original marker into the new mesh
    this->setCellMarker(RVector(mesh.cellMarker())(this->cellMarker()));
        
}
    
void Mesh::cleanNeighbourInfos(){
    //std::cout << "Mesh::cleanNeighbourInfos()"<< std::endl;
    for (Index i = 0; i < cellCount(); i ++){
        cell(i).cleanNeighbourInfos();
    }
    for (Index i = 0; i < boundaryCount(); i ++){
        boundary(i).setLeftCell(NULL);
        boundary(i).setRightCell(NULL);
    }
}

void Mesh::createNeighbourInfos(bool force){
//     double med = 0.;
//     __MS(neighboursKnown_ << " " <<force)
    if (!neighboursKnown_ || force){
        this->cleanNeighbourInfos();

//         Stopwatch sw(true);
        
        for (Index i = 0; i < cellCount(); i ++){
//            if (i%10000 ==0) __MS(i); 
            createNeighbourInfosCell_(&cell(i));
//             med+=sw.duration(true);
        }
        neighboursKnown_ = true;
    } else {
//         __M
    }
    
//     std::cout << med << " " << med/cellCount() << std::endl;
}

void Mesh::createNeighbourInfosCell_(Cell *c){
    
    for (Index j = 0; j < c->boundaryCount(); j++){
        if (c->neighbourCell(j)) continue;
        
        c->findNeighbourCell(j);
        std::vector < Node * > nodes(c->boundaryNodes(j));
//                 std::cout << findBoundary(nodes) << std::endl;

        Boundary * bound = createBoundary(nodes, 0);

//         Boundary * bound = findBoundary(*nodes[0], *nodes[1], *nodes[2]);
//         Boundary * bound = findBoundary(nodes);
//         if (!bound) {
//             bound = createBoundary_< TriangleFace >(nodes, 0, boundaryCount());
//         }
            
        bool cellIsLeft = true;
        if (bound->shape().nodeCount() == 2) {
            cellIsLeft = (c->boundaryNodes(j)[0]->id() == bound->node(0).id());
        } else if (bound->shape().nodeCount() > 2) {
            // normal vector of boundary shows outside for left cell ... every boundary needs a left cell
            if (bound->normShowsOutside(*c)){
                cellIsLeft = true;
            } else {
                cellIsLeft = false;
            }
        }

        if (bound->leftCell() == NULL && cellIsLeft) {
            if (bound->rightCell() == c){
                //* we were already here .. no need to do it again
                continue;
            }
            bound->setLeftCell(c);
            if (c->neighbourCell(j) && bound->rightCell() == NULL) bound->setRightCell(c->neighbourCell(j));
            
        } else if (bound->rightCell() == NULL){
            if (bound->leftCell() == c){
                //* we were already here .. no need to do it again
                continue;
            } else { 
                bound->setRightCell(c);
                if (c->neighbourCell(j) && bound->leftCell() == NULL ) bound->setLeftCell(c->neighbourCell(j));
            }
        }
                   
//         if (!bound->leftCell()){
//             std::cout << *bound << " " << bound->leftCell() << " " << *bound->rightCell() << std::endl;
//             throwError(1, WHERE + " Ooops, crosscheck -- every boundary need left cell.");
//         }
                   
//         std::cout << bound->id() << " " << bound->leftCell() << " " << bound->rightCell() << std::endl;
            
        //** cross check;
        if (((bound->leftCell() != c) && (bound->rightCell() != c)) || 
            (bound->leftCell() == bound->rightCell())){
            std::cerr << *c << std::endl;
            std::cerr << *bound << std::endl;
            std::cerr << bound->leftCell() << " " << bound->rightCell() << std::endl;
            if (bound->leftCell()){
                std::cerr << *bound->leftCell() << std::endl;
            }
            if (bound->rightCell()){
                std::cerr << *bound->rightCell() << std::endl;
            }
                
            
            throwError(1, WHERE + " Ooops, crosscheck --this should not happen.");
        } else {
//                     std::cout << nBounds << std::endl;
//                     std::cerr << bound->leftCell() << " " << bound->rightCell() << std::endl;
        }
    } // for_each boundary in cell
}

void Mesh::create1DGrid(const RVector & x){
    this->clear();
    this->setDimension(1);
    if (unique(sort(x)).size() != x.size()) {
        std::cerr << WHERE_AM_I << "Warning! there are non-unique values in pos" << std::endl;
    }

    if (x.size() > 1){
        this->createNode(x[0], 0.0, 0.0);
        for (Index i = 1; i < x.size(); i ++){
            this->createNode(x[i], 0.0, 0.0);
            std::vector < Node * > nodes(2);
            nodes[0] = & node(nodeCount() - 2);
            nodes[1] = & node(nodeCount() - 1);
            this->createCell(nodes);
        }
        this->createNeighbourInfos();

        for (Index i = 0; i < boundaryCount(); i ++){
            if (boundary(i).leftCell() == NULL || boundary(i).rightCell() == NULL){
                if (boundary(i).node(0).pos()[0] == x[0]) boundary(i).setMarker(1);
                else if (boundary(i).node(0).pos()[0] == x[x.size()-1]) boundary(i).setMarker(2);
            }
        }

    } else {
        std::cerr << WHERE_AM_I << "Warning! there are too few positions given: "
                << x.size() << std::endl;
    }
}

void Mesh::create2DGrid(const RVector & x, const RVector & y, int markerType){

    this->clear();
    this->setDimension(2);
    if (unique(sort(x)).size() != x.size()) {
        std::cerr << WHERE_AM_I << "Warning! there are non-unique values in pos" << std::endl;
    }

    if (unique(sort(y)).size() != y.size()) {
        std::cerr << WHERE_AM_I << "Warning! there are non-unique values in pos" << std::endl;
    }

    int marker = 0;
    if (x.size() > 1 && y.size() > 1){
        for (Index i = 0; i < y.size(); i ++){
            if (i > 0 && markerType == 2) marker++;

            for (Index j = 0; j < x.size(); j ++){
                this->createNode(x[j], y[i], 0.0);

                if (i > 0 && j > 0){
                    if (markerType == 1 || markerType == 12) marker++;
                    std::vector < Node * > nodes(4);
                    nodes[3] = & node(this->nodeCount() - 2);
                    nodes[2] = & node(this->nodeCount() - 1);
                    nodes[1] = & node(this->nodeCount() - 1 - x.size());
                    nodes[0] = & node(this->nodeCount() - 2 - x.size());

                    this->createCell(nodes, marker);
                  
//                     this->createTriangle(*nodes[1], *nodes[2], *nodes[3], marker);
//                     this->createTriangle(*nodes[0], *nodes[1], *nodes[3], marker);
                    
                }
            }
            if (markerType == 1) marker = 0;
        }
        this->createNeighbourInfos();

        for (Index i = 0; i < boundaryCount(); i ++){
            if (boundary(i).leftCell() == NULL || boundary(i).rightCell() == NULL){
                // Left
                if (std::abs(boundary(i).norm()[0] + 1.0) < TOLERANCE) boundary(i).setMarker(1);
                // Right
                else if (std::abs(boundary(i).norm()[0] - 1.0) < TOLERANCE) boundary(i).setMarker(2);
                // Top
                else if (std::abs(boundary(i).norm()[1] - 1.0) < TOLERANCE) boundary(i).setMarker(3);
                // Bottom
                else if (std::abs(boundary(i).norm()[1] + 1.0) < TOLERANCE) boundary(i).setMarker(4);
            }
        }

    } else {
        std::cerr << WHERE_AM_I << "Warning! there are too few positions given: "
            << x.size() << " " << y.size() << std::endl;
    }
}

void Mesh::create3DGrid(const RVector & x, const RVector & y, const RVector & z, int markerType){

    this->clear();
    this->setDimension(3);
    if (unique(sort(x)).size() != x.size()) {
        std::cerr << WHERE_AM_I << "Warning! there are non-unique values in pos" << std::endl;
    }

    if (unique(sort(y)).size() != y.size()) {
        std::cerr << WHERE_AM_I << "Warning! there are non-unique values in pos" << std::endl;
    }

    if (unique(sort(z)).size() != z.size()) {
        std::cerr << WHERE_AM_I << "Warning! there are non-unique values in pos" << std::endl;
    }

    int marker = 0;
    if (x.size() > 1 && y.size() > 1 && z.size() > 1){
        for (Index k = 0; k < z.size(); k ++){
            
            if (k > 0 && markerType == 3) marker++; //** count only increasing z
            
            for (Index j = 0; j < y.size(); j ++){
                
                if (j > 0 && markerType == 2) marker++;  //** count increasing y or yz
                if (j > 0 && k > 0 && markerType == 23) marker++;  //** count increasing y or yz
                
                for (Index i = 0; i < x.size(); i ++){ //**  count increasing x, yz, xz or xyz
                    
                    this->createNode(x[i], y[j], z[k]);
                    
                    if (i > 0 && j > 0 && k > 0){
                        if (markerType == 1 || markerType == 12 || markerType == 13 || markerType == 123) marker++; //** increasing y

                        std::vector < Node * > nodes(8);

                        nodes[7] = & node(this->nodeCount() - 2);
                        nodes[6] = & node(this->nodeCount() - 1);
                        nodes[5] = & node(this->nodeCount() - 1 - x.size());
                        nodes[4] = & node(this->nodeCount() - 2 - x.size());

                        nodes[3] = & node(this->nodeCount() - 2 - x.size() * y.size());
                        nodes[2] = & node(this->nodeCount() - 1 - x.size() * y.size());
                        nodes[1] = & node(this->nodeCount() - 1 - x.size() - x.size() * y.size());
                        nodes[0] = & node(this->nodeCount() - 2 - x.size() - x.size() * y.size());
                        this->createCell(nodes, marker);
                    } //** first row/column/layer
                } //** x loop (i)
                if (markerType == 1) marker = 0;
                if (j > 0 && k > 0 && markerType == 13) marker -= (x.size() - 1);
            } //** y loop (j)
//            if (k > 0 && markerType == 13) marker -= (x.size() - 1) * (y.size() - 1);
            if (k > 0 && markerType == 13) marker += (x.size() - 1);
            if (markerType == 2 || markerType == 12) marker = 0;
        } //** z loop (k)
        this->createNeighbourInfos();

         for (Index i = 0; i < boundaryCount(); i ++){
            if (boundary(i).leftCell() == NULL || boundary(i).rightCell() == NULL){
                // Left
                if (std::abs(boundary(i).norm()[0] + 1.0) < TOLERANCE) boundary(i).setMarker(1);
                // Right
                else if (std::abs(boundary(i).norm()[0] - 1.0) < TOLERANCE) boundary(i).setMarker(2);
                // Top
                else if (std::abs(boundary(i).norm()[1] - 1.0) < TOLERANCE) boundary(i).setMarker(3);
                // Bottom
                else if (std::abs(boundary(i).norm()[1] + 1.0) < TOLERANCE) boundary(i).setMarker(4);
                // Front
                else if (std::abs(boundary(i).norm()[1] - 1.0) < TOLERANCE) boundary(i).setMarker(5);
                // Back
                else if (std::abs(boundary(i).norm()[2] + 1.0) < TOLERANCE) boundary(i).setMarker(6);
            }
        }
        

    } else {
        std::cerr << WHERE_AM_I << "Warning! there are too few positions given: "
            << x.size() << " " << y.size() << " " << z.size() << std::endl;
    }
}

void Mesh::createMeshByBoundaries(const Mesh & mesh, const std::vector < Boundary * > & bounds){
    this->clear();
    this->setDimension(mesh.dim());

    std::map < int, Node* > nodeMap;

    //** Create new nodes
    for (size_t i = 0; i < bounds.size(); i ++){
        MeshEntity * ent = bounds[i];
        for (Index j = 0; j < ent->nodeCount(); j ++){
             if (nodeMap.count(ent->node(j).id()) == 0){
                 nodeMap[ent->node(j).id()] = this->createNode(ent->node(j));
             }
        }
    }

    //! Create new boundaries
    for (size_t i = 0; i < bounds.size(); i ++){
        MeshEntity * ent = bounds[i];
        std::vector < Node * > nodes(ent->nodeCount());
        for (Index j = 0; j < nodes.size(); j ++){
            nodes[j] = nodeMap[ent->node(j).id()];
        }

        createBoundary(nodes, bounds[i]->marker());
    }

}

void Mesh::createMeshByCellIdx(const Mesh & mesh, IndexArray & idxListIn){
    this->clear();
    this->setDimension(mesh.dim());

    std::map < int, Node* > nodeMap;
    
    IndexArray idxList = unique(sort(idxListIn));

    if (idxList.size() != idxListIn.size()){
        std::cerr << "This should not happen: double values in idxListIn: " 
                  << str(idxListIn.size()) << " " 
                  << str(idxList.size()) << std::endl;
    }
    
    //** Create new nodes
    for (Index i = 0; i < idxList.size(); i ++){
        
        Cell * cell = & mesh.cell(idxList[i]);
        
        for (Index j = 0; j < cell->nodeCount(); j ++){
            if (nodeMap.count(cell->node(j).id()) == 0){

                nodeMap[cell->node(j).id()] =
                        this->createNode(cell->node(j).pos(),
                                         cell->node(j).marker());
            }
        }
    }

    //! Create new cells
    for (Index i = 0; i < idxList.size(); i ++){
        Cell * cell = &mesh.cell(idxList[i]);
        std::vector < Node * > nodes(cell->nodeCount());

        for (Index j = 0; j < nodes.size(); j ++){
            nodes[j] = nodeMap[cell->node(j).id()];
        }
        
        createCell(nodes, cell->marker());
    }

    //! copy all boundary with marker != 0
    for (Index i = 0, imax = mesh.boundaryCount(); i < imax; i ++){
        Boundary * bound = &mesh.boundary(i);

        if (bound->marker() != 0){
            bool inside = true;
            std::vector < Node * > nodes(bound->nodeCount());

            for (Index j = 0, jmax = bound->nodeCount(); j < jmax; j ++){
                if (nodeMap.find(bound->node(j).id()) != nodeMap.end()) {
                    nodes[j] = nodeMap[bound->node(j).id()];
                } else {
                    inside = false;
                    break;
                }
            }
            if (inside){
                //! check that all nodes have a common cell
                if (findCommonCell(nodes, false)){
                    createBoundary(nodes, bound->marker());
                }
            }
        }
    }

    //! Copy data if available
    for (std::map < std::string, RVector >::const_iterator 
         it = mesh.dataMap().begin(); it != mesh.dataMap().end(); it ++){
        
        if (it->second.size() == mesh.cellCount()){
            addData(it->first, it->second(idxList));
        } else if (it->second.size() == mesh.nodeCount()){
            THROW_TO_IMPL
        } else if (it->second.size() == mesh.boundaryCount()){
            THROW_TO_IMPL
        }
    }
        
    //! Create all remaining boundaries
    createNeighbourInfos();
}

void Mesh::createMeshByMarker(const Mesh & mesh, int from, int to){
    if (to == -1) to = MAX_INT;
    else if (to == 0) to = from + 1;

    IndexArray cellIdx;

    for (Index i = 0; i < mesh.cellCount(); i ++){
        if (mesh.cell(i).marker() >= from && mesh.cell(i).marker() < to){
            cellIdx.push_back(i);
        }
    }
    createMeshByCellIdx(mesh, cellIdx);
}

void Mesh::addExportData(const std::string & name, const RVector & data) {
  //  std::cout << "add export Data: " << name << " " << min(data) << " "  << max(data) << std::endl;
    if (exportDataMap_.count(name)){
        exportDataMap_[name] = data;
    } else {
        exportDataMap_.insert(std::make_pair(name, data));
    }
}

RVector Mesh::exportData(const std::string & name) const {
    if (exportDataMap_.count(name)){
        return exportDataMap_.find(name)->second;
    } else {
        throwError(1, " Warning!! requested export 'data' vector " + name +
        " does not exist.");
    }
    return RVector(0);
}

void Mesh::clearExportData(){
    exportDataMap_.clear();
}

void Mesh::dataInfo() const{
    if (exportDataMap_.empty()){
        std::cout << "No data." << std::endl;
    } else {
        for (std::map < std::string, RVector >::const_iterator
            it = exportDataMap_.begin(); it != exportDataMap_.end(); it ++){
            std::cout << it->first << ": " << str(it->second.size()) << std::endl;
        }
    }
}

IVector Mesh::nodeMarker() const {
    IVector tmp(nodeCount());
    std::transform(nodeVector_.begin(), nodeVector_.end(), tmp.begin(), std::mem_fun(& Node::marker));
    return tmp;
}

IVector Mesh::boundaryMarker() const {
    IVector tmp(boundaryCount());
    std::transform(boundaryVector_.begin(), boundaryVector_.end(), tmp.begin(),
                   std::mem_fun(&Boundary::marker));
    return tmp;
}

IVector Mesh::cellMarker() const{
    IVector tmp(cellCount());
    std::transform(cellVector_.begin(), cellVector_.end(), tmp.begin(),
                    std::mem_fun(&Cell::marker));
    return tmp;
}

RVector Mesh::cellAttributes() const{
    #ifdef _MSC_VER
	std::vector < double > tmp(cellCount());
    std::transform(cellVector_.begin(), cellVector_.end(), tmp.begin(),
                    std::mem_fun(&Cell::attribute));
	return tmp;
	#else
	RVector tmp(cellCount());
    std::transform(cellVector_.begin(), cellVector_.end(), tmp.begin(),
                    std::mem_fun(&Cell::attribute));
    return tmp;
	#endif
}

void Mesh::setCellAttributes(const RVector & attr){
    if (attr.size() != (uint)cellCount()){
        throwError(1, WHERE_AM_I + " std::vector attr.size() != cellCount()" + toStr(attr.size()) + " " + toStr(cellCount()));
    }
    for (Index i = 0; i < cellCount(); i ++) cell(i).setAttribute(attr[i]);
}

void Mesh::setCellAttributes(double attr){
    for (Index i = 0; i < cellCount(); i ++) cell(i).setAttribute(attr);
}

void Mesh::mapCellAttributes(const std::map < float, float > & aMap){
    std::map< float, float >::const_iterator itm;

    if (aMap.size() != 0){
        for (Index i = 0, imax = cellCount(); i < imax; i++){
            itm = aMap.find(float(cell(i).marker()));
            if (itm != aMap.end()) cell(i).setAttribute((*itm).second);
        }
    }
}

void Mesh::mapAttributeToParameter(const IndexArray & cellMapIndex,
                                    const RVector & attributeMap, double defaultVal){
    DEPRECATED
//     RVector mapModel(attributeMap.size() + 2);
//     mapModel[1] = defaultVal;
//   //mapModel[std::slice(2, (int)attributeMap.size(), 1)] = attributeMap;
//     for (Index i = 0; i < attributeMap.size(); i ++) mapModel[i + 2] = attributeMap[i];
//
//     std::vector< Cell * > emptyList;
//
//     for (Index i = 0, imax = cellCount(); i < imax; i ++){
//         if (cellMapIndex[i] > (int)mapModel.size() -1){
//             std::cerr << WHERE_AM_I << " cellMapIndex[i] > attributeMap " << cellMapIndex[i]
//                     << " " << mapModel.size() << std::endl;
//         }
//
//         cell(i).setAttribute(mapModel[cellMapIndex[i]]);
//
//         if (mapModel[cellMapIndex[i]] < TOLERANCE){
//             emptyList.push_back(&cell(i));
//         }
//     }

    //fillEmptyCells(emptyList);
}

void Mesh::mapBoundaryMarker(const std::map < int, int > & aMap){
    std::map< int, int >::const_iterator itm;
    if (aMap.size() != 0){
        for (Index i = 0, imax = boundaryCount(); i < imax; i++){
            itm = aMap.find(boundary(i).marker());
            if (itm != aMap.end()){
	       boundary(i).setMarker((*itm).second);
            }
        }
    }
}

void Mesh::fillEmptyCells(const std::vector< Cell * > & emptyList, double background){
    if (emptyList.size() == 0) return;

    if (background != -1.0){
        for (size_t i = 0; i < emptyList.size(); i ++) emptyList[i]->setAttribute(background);
        return;
    }

    bool smooth = false;
    bool horizontalWeight = true;
    
    createNeighbourInfos();
    if (emptyList.size() > 0){
        if (debug())std::cout << "Prolongate " << emptyList.size() << " empty cells. (" << this->cellCount() << ")" << std::endl;
        std::vector< Cell * > nextVector;
        Cell * cell;

        std::map< Cell*, double > prolongationMap;

        RVector3 XY(1., 1., 0.);
        if (this->dim() == 2) XY[1] = 0.0;
        
        for (size_t i = 0; i < emptyList.size(); i ++){
            cell = emptyList[i];

            double weight = 0.0;
            double val = 0.0;
            for (Index j = 0; j < cell->neighbourCellCount(); j ++){
                Cell * nCell = cell->neighbourCell(j);
                if (nCell){
                    if (nCell->attribute() > TOLERANCE){
                        if (horizontalWeight){
                            Boundary * b=findCommonBoundary(*nCell, *cell);
                            if (b){
                                double zWeight = (b->norm()*XY).abs() + 1e-6;
                                val += nCell->attribute() * zWeight;
                                weight += zWeight;
                            }
                        } else {
                            val += nCell->attribute();
                            weight += 1.0;
                        }
                    }
                }
            }
            if (weight < 1e-8) {
                nextVector.push_back(cell);
            } else {
                if (smooth){
                    cell->setAttribute(val / weight);
                } else {
                    prolongationMap[cell] = val / weight;
                }
            }
        }

        if (!smooth){//**apply std::map< uint, val > prolongationMap;
            for (std::map< Cell *, double >::iterator it= prolongationMap.begin();
                it != prolongationMap.end(); it ++){
                    it->first->setAttribute(it->second);
            }
        }
        if (emptyList.size() == nextVector.size()){
            this->exportVTK("fillEmptyCellsFail");
            std::cerr << WHERE_AM_I << " WARNING!! cannot fill emptyList: see fillEmptyCellsFail.vtk"<< std::endl;
            std::cerr << "trying to fix"<< std::endl;

            for (size_t i = 0; i < emptyList.size(); i ++) emptyList[i]->setAttribute(mean(this->cellAttributes()));
            nextVector.clear();
        }
        fillEmptyCells(nextVector, background);
    } //** if emptyList.size() > 0
}

Mesh & Mesh::scale(const RVector3 & s){
    std::for_each(nodeVector_.begin(), nodeVector_.end(),
                  boost::bind(& Node::scale, _1, boost::ref(s)));
    std::for_each(holeMarker_.begin(), holeMarker_.end(),
                  boost::bind(& RVector3::scale, _1, boost::ref(s)));
    std::for_each(regionMarker_.begin(), regionMarker_.end(),
                  boost::bind(& RVector3::scale, _1, boost::ref(s)));
    
    rangesKnown_ = false;
    return *this;
}

Mesh & Mesh::translate(const RVector3 & t){
    std::for_each(nodeVector_.begin(), nodeVector_.end(),
                   boost::bind(& Node::translate, _1, boost::ref(t)));
    std::for_each(holeMarker_.begin(), holeMarker_.end(),
                  boost::bind(& RVector3::translate, _1, boost::ref(t)));
    std::for_each(regionMarker_.begin(), regionMarker_.end(),
                  boost::bind(& RVector3::translate, _1, boost::ref(t)));
    
    rangesKnown_ = false;
    return *this;
}

Mesh & Mesh::rotate(const RVector3 & r){
    std::for_each(nodeVector_.begin(), nodeVector_.end(),
                   boost::bind(& Node::rotate, _1, boost::ref(r)));
    std::for_each(holeMarker_.begin(), holeMarker_.end(),
                  boost::bind(& RVector3::rotate, _1, boost::ref(r)));
    std::for_each(regionMarker_.begin(), regionMarker_.end(),
                  boost::bind(& RVector3::rotate, _1, boost::ref(r)));
    
    rangesKnown_ = false;
    return *this;
}

void Mesh::relax(){
   THROW_TO_IMPL
   //  int E = 0;

//     for (int T = 6; T >= 3; T--){
//       for (int s = 0; s < Ns; s++){
// 	if (side[s].mark == 0){
// 	  if ((node[side[s].a].mark == 0) &&
// 	       (node[side[s].b].mark == 0) &&
// 	       (node[side[s].c].mark == 0) &&
// 	       (node[side[s].d].mark == 0)) {
// 	    E = node[side[s].c].Nne + node[side[s].d].Nne - node[side[s].a].Nne - node[side[s].b].Nne;

// 	    if (E == T) {
// 	      node[side[s].a].Nne++; node[side[s].b].Nne++;
// 	      node[side[s].c].Nne--; node[side[s].d].Nne--;
// 	      swap_side(s);
// 	    }
// 	  }
// 	}
//       }
//     }
  }

void Mesh::smooth(bool nodeMoving, bool edgeSliding, uint smoothFunction, uint smoothIteration){
    createNeighbourInfos();

    for (Index j = 0; j < smoothIteration; j++){
//         if (edgeSwapping) {
//             for (Index i = 0; i < boundaryCount(); i++) dynamic_cast< Edge & >(boundary(i)).swap(1);
//         }
        if (nodeMoving) {
            for (Index i = 0; i < nodeCount(); i++) {
                bool forbidMove = (node(i).marker() != 0);

                std::pair < Boundary *, Boundary * > slide(0, 0);
                bool noSlide = false;
                for (std::set< Boundary * >::iterator it = node(i).boundSet().begin();
                     it != node(i).boundSet().end(); it ++){
                   
                    if ((*it)->marker() != 0){
                        if (slide.first == 0){
                            slide.first = (*it);
                        } else {
                            
                            if (slide.second == 0){
                                if (slide.first->norm() == (*it)->norm()){
                                    slide.second = (*it);
                                } else {
                                    // two marker bounds with different norm -> corner
                                    noSlide = true;
                                }
                            } else {
                                // more than two marker bounds -> corner
                                noSlide = true;
                            }
                        }
                    }
                
                    if (edgeSliding) {
                        forbidMove = forbidMove || noSlide;
                    }else {
                        forbidMove = forbidMove || (*it)->marker() != 0;
                    }
                
                    if (!edgeSliding){
                        forbidMove = forbidMove || ((*it)->leftCell() == NULL || (*it)->rightCell() == NULL);
                    }
                    
                    if (forbidMove) break;
                }
                    
                if (!forbidMove) {
                    if (slide.first && slide.second){
                    // move weighted with itself as double weight .. results in slight slide
                        node(i).setPos((
                                slide.first->node(0).pos() + 
                                slide.first->node(1).pos() +
                                slide.second->node(0).pos() +
                                slide.second->node(1).pos()) / 4.0);
                    } else {
                        node(i).smooth(smoothFunction);
                    }
                }
            }
        }
    }
}

void Mesh::fillKDTree_() const {
MEMINFO
    if (!tree_) tree_ = new KDTreeWrapper();

    if (tree_->size() != nodeCount()){
        if (tree_->size() == 0){
//            for (Index i = 0; i < nodeCount(); i ++) tree_->insert(nodeVector_[i]);
            for_each(nodeVector_.begin(), nodeVector_.end(), boost::bind(&KDTreeWrapper::insert, tree_, _1));

            tree_->tree()->optimize();
        } else {
            throwError(1, WHERE_AM_I + toStr(this) + " kd-tree is only partially filled: this should no happen: nodeCount = " + toStr(nodeCount())
                                      + " tree-size() " + toStr(tree_->size()));
        }
    }
MEMINFO
}

void Mesh::addRegionMarker(const RegionMarker & reg){
    regionMarker_.push_back(reg);
}
    
void Mesh::addRegionMarker(const RVector3 & pos, int marker, double area){
    if (area < 0) {
        addHoleMarker(pos);
    } else {
        regionMarker_.push_back(RegionMarker(pos, marker, area));
    }
}
 
void Mesh::addHoleMarker(const RVector3 & pos){
    holeMarker_.push_back(pos);
}
    
RSparseMapMatrix & Mesh::cellToBoundaryInterpolation() const {
    if (!cellToBoundaryInterpolationCache_){
        if (!neighboursKnown_){
            throwError(1, "Please call once createNeighbourInfos() for the given mesh.");
        }
        
        cellToBoundaryInterpolationCache_ = new RSparseMapMatrix(this->boundaryCount(), 
                                                                 this->cellCount());    
        for (Index i = 0; i < boundaryCount(); i ++){
            Boundary * b = boundaryVector_[i];
            Cell * lC = b->leftCell();
            Cell * rC = b->rightCell();
    
            double df1 = 0.0;
            double df2 = 0.0;
            
            bool harmonic = false;
            if (lC) df1 = b->center().distance(lC->center());
            if (rC) df2 = b->center().distance(rC->center());
            double d12 = (df1 + df2);
            
            if (lC && rC){
                if (harmonic){
                } else {
                    cellToBoundaryInterpolationCache_->addVal(b->id(), lC->id(), df2/d12);
                    cellToBoundaryInterpolationCache_->addVal(b->id(), rC->id(), -df2/d12 + 1.0);
                }
            } else if (lC){
                cellToBoundaryInterpolationCache_->addVal(b->id(), lC->id(), 1.0);
            } else {
                throwError(1, WHERE_AM_I + " this should not happen");
            }
        }
    } else {
        if (!staticGeometry_){
            delete cellToBoundaryInterpolationCache_;
            cellToBoundaryInterpolationCache_ = 0;
            return this->cellToBoundaryInterpolation();
        }
    }
        
    return *cellToBoundaryInterpolationCache_;
}

std::vector< RVector3 > 
Mesh::cellDataToBoundaryGradient(const RVector & cellData) const {
    return cellDataToBoundaryGradient(cellData, 
      boundaryDataToCellGradient(this->cellToBoundaryInterpolation()*cellData));
}

std::vector< RVector3 > 
Mesh::cellDataToBoundaryGradient(const RVector & cellData,
        const std::vector< RVector3 > & cellGrad) const{
    if (!neighboursKnown_){
        throwError(1, "Please call once createNeighbourInfos() for the given mesh.");
    }
    std::vector< RVector3 > ret(boundaryCount());
    
    for (Index i = 0; i < boundaryCount(); i ++){
        Boundary * b = boundaryVector_[i];
        Cell * lC = b->leftCell();
        Cell * rC = b->rightCell();
    
        RVector3 grad(0.0, 0.0, 0.0);
        RVector3 tangent((b->node(1).pos() - b->node(0).pos()).norm());
        
        if (lC && rC){
            double df1 = b->center().distance(lC->center());
            double df2 = b->center().distance(rC->center());
                    
            ret[b->id()] = b->norm()*(cellData[rC->id()] - cellData[lC->id()]) / (df1+df2);
                        
            ret[b->id()] += tangent * (tangent.dot(cellGrad[lC->id()]) + 
                                       tangent.dot(cellGrad[rC->id()])) * 0.5;
        } else if(lC){
            ret[b->id()] = tangent * tangent.dot(cellGrad[lC->id()]);
        }
    }                 
    return ret;
}

std::vector< RVector3 > 
Mesh::boundaryDataToCellGradient(const RVector & v) const{
    if (!neighboursKnown_){
        throwError(1, "Please call once createNeighbourInfos() for the given mesh.");
    }
    std::vector < RVector3 > ret(this->cellCount());

    const std::vector < RVector3 > & flow = this->boundarySizedNormals();
    RVector3 vec(0.0, 0.0, 0.0);
    for (Index i = 0; i < this->boundaryCount(); i ++ ){
        Boundary * b = this->boundaryVector_[i];
        vec = flow[b->id()] * v[b->id()];
        
        if (b->leftCell()){
            ret[b->leftCell()->id()] += vec;
        } 
        if (b->rightCell()){
            ret[b->rightCell()->id()] -= vec;
        } 
    }
    for (Index i = 0; i < ret.size(); i ++ ){
        ret[i] /= cellVector_[i]->size();
    }
    return ret;
}
        
RVector Mesh::divergence(const std::vector < RVector3 > & V) const{
    RVector ret(this->cellCount());
    if (!neighboursKnown_){
        throwError(1, "Please call once createNeighbourInfos() for the given mesh.");
    }
    
    const std::vector < RVector3 > & flow = this->boundarySizedNormals();
    double vec = 0;
    for (Index i = 0; i < this->boundaryCount(); i ++ ){
        Boundary * b = this->boundaryVector_[i];
        vec = flow[b->id()].dot(V[b->id()]);
        
        if (b->leftCell()){
            ret[b->leftCell()->id()] += vec;
        } 
        if (b->rightCell()){
            ret[b->rightCell()->id()] -= vec;
        } 
    }
    return ret / this->cellSizes();
}

        
    
} // namespace GIMLI
