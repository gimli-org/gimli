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

#ifndef _GIMLI_MESH__H
#define _GIMLI_MESH__H

#include "gimli.h"
#include "meshentities.h"
#include "node.h"
#include "pos.h"

#include <vector>
#include <list>
#include <set>
#include <map>
#include <fstream>

namespace GIMLI{

class KDTreeWrapper;

template < class T > class DLLEXPORT BoundingBox;
typedef BoundingBox< double > RBoundingBox;
typedef BoundingBox< double > IBoundingBox;

//! A BoundingBox
/*! A BoundingBox which contains a min and max Vector3< ValueType >. IBoundingBox and RBoundingBox typedefs an Integer and Real Boundingbox, respectivly. */
template < class ValueType > class DLLEXPORT BoundingBox{
public:
    /*! Default constructor, with BoundingBox[0.0, 0.0, 0.0; 1.0, 1.0, 1.0] */
    BoundingBox(const Pos < ValueType > & min = Pos < double >(0, 0, 0), const Pos < ValueType > & max = Pos < double >(1.0, 1.0, 1.0))
        : min_(min), max_(max){
    }

    /*! Construct BBox from position vector */
    BoundingBox(const std::vector < Pos < ValueType > > & vPos){
        min_ = Pos < ValueType >((ValueType)MAX_DOUBLE, (ValueType)MAX_DOUBLE, (ValueType)MAX_DOUBLE);
        max_ = min_ * -1.0;
        for (uint i = 0; i < vPos.size(); i ++){
            min_[ 0 ] = std::min(vPos[ i ][ 0 ], min_[ 0 ]);
            min_[ 1 ] = std::min(vPos[ i ][ 1 ], min_[ 1 ]);
            min_[ 2 ] = std::min(vPos[ i ][ 2 ], min_[ 2 ]);
            max_[ 0 ] = std::max(vPos[ i ][ 0 ], max_[ 0 ]);
            max_[ 1 ] = std::max(vPos[ i ][ 1 ], max_[ 1 ]);
            max_[ 2 ] = std::max(vPos[ i ][ 2 ], max_[ 2 ]);
        }
    }

    /*! Copy constructor. */
    BoundingBox(const BoundingBox < ValueType > & bbox){ this->copy_(bbox); }

    /*! Assignment operator. */
    BoundingBox < ValueType > & operator = (const BoundingBox < ValueType > bbox){
    	if (this != & bbox){
            this->copy_(bbox);
        } return * this;
    }

    /*! Set minimum Position. */
    void setMin(const Pos < ValueType > & min) { min_ = min; }
    /*! Returns a copy of the minimum position. */
    const Pos < ValueType > min() const { return min_; }

    /*! Set maximum Position. */
    void setMax(const Pos < ValueType > & max) { max_ = max; }
    /*! Returns a copy of the maximum position. */
    const Pos < ValueType > max() const { return max_; }

protected:

    void copy_(const BoundingBox < ValueType > & bbox){
    	min_ = bbox.min();
        max_ = bbox.max();
    }

    Pos < ValueType > min_;
    Pos < ValueType > max_;
};

template < class ValueType > std::ostream & operator << (std::ostream & str, const BoundingBox< ValueType > & bb){
    str << "BoundingBox [ " << bb.min() << ", " << bb.max() << " ]" ;
    return str;
}

DLLEXPORT std::ostream & operator << (std::ostream & str, const Mesh & mesh);

class DLLEXPORT Mesh {
public:

    /*! Default constructor, create empty mesh with dimension dim */
    Mesh(uint dim=2);

    /*! Constructor, read mesh from filename */
    Mesh(const std::string & filename);

    /*! Copy constructor. */
    Mesh(const Mesh & mesh);

    /*! Copy assignment operator. */
    Mesh & operator = (const Mesh & mesh);

    /*! Default destructor. */
    ~Mesh();

    void clear();

    //** start creation stuff
    Node * createNode(double x, double y, double z, int marker = 0);

    Node * createNode(const Node & node);

    Node * createNode(const RVector3 & pos, int marker = 0);

    Node * createNodeWithCheck(const RVector3 & pos, double tol = 1e-6);

    Boundary * createBoundary(std::vector < Node * > & nodes, int marker = 0);
    Boundary * createBoundary(const Boundary & bound);
    Boundary * createBoundary(const Cell & cell);
    Boundary * createNodeBoundary(Node & n1, int marker = 0);
    Boundary * createEdge(Node & n1, Node & n2, int marker = 0);
    Boundary * createEdge3(Node & n1, Node & n2, Node & n3, int marker = 0);
    Boundary * createTriangleFace(Node & n1, Node & n2, Node & n3, int marker = 0);
    Boundary * createQuadrangleFace(Node & n1, Node & n2, Node & n3, Node & n4, int marker = 0);

    /*! Create empty cell without a node or a shape. */
    Cell * createCell(int marker = 0);
    Cell * createCell(std::vector < Node * > & nodes, int marker = 0);
    Cell * createCell(const Cell & cell);
    Cell * createTriangle(Node & n1, Node & n2, Node & n3, int marker = 0);
    Cell * createQuadrangle(Node & n1, Node & n2, Node & n3, Node & n4, int marker = 0);
    Cell * createTetrahedron(Node & n1, Node & n2, Node & n3, Node & n4, int marker = 0);

    /*! Create a cell as a copy of a cell from an alternative mesh. Each Node of cell will be created with check. */
    Cell * copyCell(const Cell & cell);
    /*! Create a boundary as a copy of a boundary from an alternative mesh. Each Node of cell will be created with check. */
    Boundary * copyBoundary(const Boundary & bound);

    void create1DGrid(const RVector & x);

    void create2DGrid(const RVector & x, const RVector & y, int markerType = 0);

    void create3DGrid(const RVector & x, const RVector & y, const RVector & z, int markerType = 0);

    /*! Create one dimensional grid. Boundary on the domain border will get marker = 1 .*/
    void createGrid(const RVector & x) { create1DGrid(x); }

    /*! Create two dimensional grid. Boundary on the domain border will get marker = 1 .*/
    void createGrid(const RVector & x, const RVector & y, int markerType=0) {
        create2DGrid(x, y, markerType);
    }
    /*! Create three dimensional grid. Boundary on the domain border will get marker = 1 .*/
    void createGrid(const RVector & x, const RVector & y, const RVector & z, int markerType=0){
        create3DGrid(x, y, z, markerType);
    }

    /*! Fill this 3D mesh with 3D boundary elements from the 2D mesh cells. Increase mesh dimension. Mesh should contain 2D cells. */
    void createHull(const Mesh & mesh);

    void createClosedGeometry(const std::vector < RVector3 > & vPos, int nSegments, double dxInner);

    void createClosedGeometryParaMesh(const std::vector < RVector3 > & vPos, int nSegments, double dxInner);

    void createClosedGeometryParaMesh(const std::vector < RVector3 > & vPos, int nSegments, double dxInner,
                                        const std::vector < RVector3 > & addit);

    /*! Create and copy global H2 mesh of this mesh.*/
    Mesh createH2() const;

    /*! Create and copy global P2 mesh of this mesh.*/
    Mesh createP2() const;

    /*! Create a partly mesh from mesh, based on cell-ids */
    void createMeshByCellIdx(const Mesh & mesh, std::vector < int > & idxList);

    /*! Create a partly mesh without cells from mesh, based on a vector of ptrs to boundaries */
    void createMeshByBoundaries(const Mesh & mesh, const std::vector < Boundary * > & bounds);

    /*! Create a partly mesh from mesh, based on meshs attributes. For a single attribute set to to 0, for unlimited set to to -1 */
    void createMeshByMarker(const Mesh & mesh, int from, int to = -1);
    //** end creation stuff

    //! Show some infos
    void showInfos(){ std::cout << *this << std::endl; }

    const std::vector< Node * > & nodes() const { return nodeVector_; }

    const std::vector< Cell * > & cells() const { return cellVector_; }

    const std::vector< Boundary * > & boundaries() const { return boundaryVector_; }

    inline uint nodeCount() const { return nodeVector_.size(); }
    Node & node(uint i) const;
    Node & node(uint i);

    uint cellCount() const { return cellVector_.size(); }
    Cell & cell(uint i) const;
    Cell & cell(uint i);

    uint boundaryCount() const { return boundaryVector_.size(); }
    Boundary & boundary(uint i) const;
    Boundary & boundary(uint i);

    /*! Returns a vector of all node positions */
    std::vector < RVector3 > positions() const;

    /*! Returns a vector of node positions for an index vector */
    std::vector < RVector3 > positions(const IndexArray & idx) const;

    /*! Returns a vector of all cell center positions*/
    std::vector < RVector3 > cellCenters() const;
    std::vector < RVector3 > cellCenter() const { return cellCenters(); }

    /*! Returns a RVector of all cell sizes */
    RVector cellSizes() const;

    /*! Returns a vector of all boundary marker */
    std::vector < int > boundaryMarker() const;

    /*! Returns a vector of all cell marker */
    std::vector < int > cellMarker() const;

    /*! Returns a vector of all node marker */
    std::vector < int > nodeMarker() const;

    /*! Returns an index vector of all nodes that match the marker */
    IndexArray findNodesIdxByMarker(int marker) const;

//     /*! Returns an index list of all nodes that match the marker */
//     std::list < uint > findListNodesIdxByMarker(int marker) const;

    /*! Returns a vector of boundary ptrs with the boundary marker equal marker.*/
    std::vector < Boundary * > findBoundaryByMarker(int marker) const;

    /*! Returns a vector of boundary ptrs with the boundary marker between [from and to). \n
        for to equal open end set to = MAX_INT */
    std::vector < Boundary * > findBoundaryByMarker(int from, int to) const;

    /*! Return ptr to the cell that match position pos, counter holds amount of touch tests.
        Searching is done first by nearest-neighbour-kd-tree search,
        followed by slope-search if extensive is set. Return NULL if no cell can be found. */
    Cell * findCell(const RVector3 & pos, size_t & counter, bool extensive) const ;

    /*! Shortcut for \ref findCell(const RVector3 & pos, size_t & counter, bool extensive) */
    Cell * findCell(const RVector3 & pos, bool extensive = true) const {
        size_t counter; return findCell(pos, counter, extensive); }

    /*! Return the index to the node of this mesh with the smallest distance to pos. */
    uint findNearestNode(const RVector3 & pos);

    /*! Returns vector of cell ptrs with marker match the range [from .. to). \n
        For single marker match to is set to 0, for open end set to = -1 */
    std::vector < Cell * > findCellByMarker(int from, int to = 0) const;

    /*! Returns vector of cell ptrs with attribute match the range [from .. to). \n
        For single attribute match to is set to 0.0, for open end set to = -1.0 */
    std::vector < Cell * > findCellByAttribute(double from, double to = 0.0) const;

//     void setAttributes(const RVector & atts) { attributes_ = atts; }
//     const RVector & attributes() const { return attributes_; }
//     RVector & attributes() { return attributes_; }
    //** end get infos stuff

    //** start mesh modification stuff
    /*! Prolongate the attribute of each cell in emptyList by the attribute from neighbouring cells.
    The attributes have to be lower than \ref TOLERANCE. This function is called recursivly until all zero-attribute-cells in emptyList are filles with an attribute greater than Zero. */
    void fillEmptyCells(const std::vector< Cell * > & emptyList, double background = -1.0);

    void recountNodes();

    void sortNodes(const std::vector < int > & perm);

    /*! Return true if createNeighbourInfos is called once */
    inline bool neighboursKnown() const { return neighboursKnown_; }
    
    /*! Remove from each boundary the ptr to the corresponding left and right cell*/
    void cleanNeighbourInfos();

    /*! Search and set to each boundary the corresponding left and right cell.*/
    void createNeighbourInfos(bool force = false);

    /*! Create and store boundaries and neighboring information for this cell.*/
    void createNeighbourInfosCell_(Cell *c);
        
    void relax();

    void smooth(bool nodeMoving, bool edgeSwapping, uint smoothFunction, uint smoothIteration);

    /*! Scale mesh with \ref RVector3 s*/
    Mesh & scale(const RVector3 & s);

    /*! Translate mesh with \ref RVector3 t*/
    Mesh & translate(const RVector3 & t);

    /*! Rotate mesh with \ref RVector3 r, r in radian, If you want to rotate in degree, use \ref degToRad(const RVector3 & deg). */
    Mesh & rotate(const RVector3 & r);
    //** end mesh modification stuff

    /*! apply a 4x4 transformation matrix to the whole mesh*/
    template < class Matrix > Mesh & transform(const Matrix & mat){
//         std::for_each(nodeVector_.begin(), nodeVector_.end(),
//                        bind2nd(std::mem_fun(&Node::pos().transform), mat));
        for (uint i = 0; i < nodeVector_.size(); i ++) nodeVector_[ i ]->pos().transform(mat);
        rangesKnown_ = false;
        return *this;
    }

    //** start I/O stuff
    int save(const std::string & fileName, IOFormat format = Binary) const;
    int saveAscii(const std::string & fileName) const;

    /*! Be carefull with interchanging binary meshs between 32-64bit architecture. Atm we save fixed int for counter and idx.
  We have to replace and test it with uint32 or uint16 */
    int saveBinary(const std::string & fileName) const;

    /*! Load Mesh from file and try to import fileformat regarding file suffix.*/
    void load(const std::string & fileName, IOFormat format = Binary);

    void loadAscii(const std::string & fileName);

    void importMod(const std::string & fileName);

    void importVTK(const std::string & fbody);

    void importVTU(const std::string & fbody);

    /*! Import Ascii STL, and save triangles as \ref Boundary(ies). */
    void importSTL(const std::string & fileName, bool isBinary = false);

    /*! Be carefull with interchanging binary meshs between 32-64bit architecture. Atm we save fixed int for counter and idx.
    We have to replace and test it with uint32 or uint16 */
    void loadBinary(const std::string & fileName);

    /*! Save mesh in binary format v.2.0. should be possible to interchange on all little endian platforms. Contains all export data, marker an attribute.
        If something goes wrong while writing, an exception is thrown.
    sizes uint8 (1 byte), uint16 (2 byte), uint32 (4 byte), uint64 (8 byte)
    Format: */
    void saveBinaryV2(const std::string & fbody) const;

    /*! Load mesh in binary format v.2.0. should be possible to interchange on all little endian platforms. Format see \ref saveBinaryV2.
        If something goes wrong while reading, an exception is thrown. */
    void loadBinaryV2(const std::string & fbody);

    int exportSimple(const std::string & fbody, const RVector & data) const ;

    /*! Very simple export filter. Write to file fileName:
        x_0 y_0 [z_0] data1_0 [data2]_0 \n
        x_n y_n [z_n] data1_n [data2]_n \n
        with n = number of cells, and x y [z] == cell center
    */
    int exportMidCellValue(const std::string & fileName, const RVector & data1,
                            const RVector & data2 = RVector(0)) const ;

    void exportVTK(const std::string & fbody,
                   const std::map< std::string, RVector > & data,
                   const std::vector < RVector3 > & vec,
                   bool writeCells=true) const;

    void exportVTK(const std::string & fbody,
                   const std::map< std::string, RVector > & data,
                   bool writeCells=true) const;

    /*! Export mesh and whole exportData map */
    void exportVTK(const std::string & fbody, bool writeCells=true) const;

    /*! Export mesh and whole exportData map and vector data in vec*/
    void exportVTK(const std::string & fbody,
                   const std::vector < RVector3 > & vec,
                   bool writeCells=true) const;

    void readVTKPoints_(std::fstream & file, const std::vector < std::string > & row);
    void readVTKCells_(std::fstream & file, const std::vector < std::string > & row);
    void readVTKScalars_(std::fstream & file, const std::vector < std::string > & row);
    void readVTKPolygons_(std::fstream & file, const std::vector < std::string > & row);

    /*! Export the mesh in filename using vtu format:
    Visualization Toolkit Unstructured Points Data (http://www.vtk.org)
    Set binary to true writes the data content in binary format.
    The file suffix .vtu will be added or substituted if .vtu or .vtk is found.
    \ref exportData, cell.marker and cell.attribute will be exported as data. */
    void exportVTU(const std::string & filename, bool binary = false) const ;

    /*! Export the boundary of this mesh in vtu format: Visualization Toolkit Unstructured Points Data (http://www.vtk.org) Set Binary to true writes the datacontent in binary format. The file suffix .vtu will be added or substituted if .vtu or .vtk is found. */
    void exportBoundaryVTU(const std::string & fbody, bool binary = false) const ;

    /*! Internal function for exporting VTU */
    void addVTUPiece_(std::fstream & file, const Mesh & mesh,
                        const std::map < std::string, RVector > & data) const;

    void exportAsTetgenPolyFile(const std::string & filename);
    //** end I/O stuff


    void addExportData(const std::string & name, const RVector & data);

    RVector exportData(const std::string & name) const;

    const std::map< std::string, RVector > & exportDataMap() const { return exportDataMap_; }

    void setExportDataMap(const std::map< std::string, RVector > & eMap) { exportDataMap_ = eMap; }

    void clearExportData();

    /*!Set the comment for VTK Ascii export headline.*/
    void setCommentString(const std::string & commentString) {commentString_ = commentString;}

    /*!Return comment for VTK Ascii export headline.*/
    const std::string & commentString() const {return commentString_;}

    //** probably deprecated
    void mapCellAttributes(const std::map < float, float > & aMap);

    void mapAttributeToParameter(const std::vector< int > & cellMapIndex,
                                    const RVector & attributeMap, double defaultVal);

    //void mapParameterToAttribute(const std::vector< int > & cellMapIndex);
    /*! Change all boundary marker that match bMap.first to bMap.second. */
    void mapBoundaryMarker(const std::map < int, int > & aMap);

    /*! Set all cell attributes to the valaues in vector attribute.*/
    void setCellAttributes(const RVector & attribute);

    /*! Set all cell attributes to the scalar value: attribute */
    void setCellAttributes(double attribute);

    /*! Return a RVector of all cell attributes */
    RVector cellAttributes() const;

    //** probably deprecated
    double xmin() const { findRange_(); return minRange_[0]; }
    double ymin() const { findRange_(); return minRange_[1]; }
    double zmin() const { findRange_(); return minRange_[2]; }
    double xmax() const { findRange_(); return maxRange_[0]; }
    double ymax() const { findRange_(); return maxRange_[1]; }
    double zmax() const { findRange_(); return maxRange_[2]; }

    const RBoundingBox boundingBox() const { findRange_(); return RBoundingBox(minRange_, maxRange_);}

    void setDimension(uint dim){ dimension_ = dim;}
    uint dimension() const { return dimension_; }
    uint dim() const { return dimension_; }

protected:
    void copy_(const Mesh & mesh);

    void findRange_() const ;

    Node * createNode_(const RVector3 & pos, int marker, int id);

    template < class B > Boundary * createBoundary_(std::vector < Node * > & nodes, int marker, int id){
        if (id == -1) id = boundaryCount();
        boundaryVector_.push_back(new B(nodes));
        boundaryVector_.back()->setMarker(marker);
        boundaryVector_.back()->setId(id);
        return boundaryVector_.back();
    }

    template < class B > Boundary * createBoundaryChecked_(std::vector < Node * > & nodes, int marker){
        Boundary * b = findBoundary(nodes);
        if (!b) {
            b = createBoundary_< B >(nodes, marker, boundaryCount());
        } else {
            if (marker != 0) b->setMarker(marker);
        }
        return b;
    }

    template < class C > Cell * createCell_(std::vector < Node * > & nodes, int marker, int id){
        if (id == -1) id = cellCount();
        cellVector_.push_back(new C(nodes));
        cellVector_.back()->setMarker(marker);
        cellVector_.back()->setId(id);
        return cellVector_.back();
    }

//    Node * createRefinementNode_(Node * n0, Node * n1, SparseMapMatrix < Node *, Index > & nodeMatrix);
    Node * createRefinementNode_(Node * n0, Node * n1, std::map< std::pair < Index, Index >, Node * > & nodeMatrix);

    void createRefined_(const Mesh & mesh, bool p2, bool r2);

    Cell * findCellBySlopeSearch_(const RVector3 & pos, Cell * start, size_t & count, bool tagging) const;

    void fillKDTree_() const;

    std::vector< Node * >     nodeVector_;
    std::vector< Boundary * > boundaryVector_;
    std::vector< Cell * >     cellVector_;

    uint dimension_;

   // RVector attributes_;

    mutable RVector3 minRange_;
    mutable RVector3 maxRange_;
    mutable bool rangesKnown_;
    bool neighboursKnown_;

    std::map< std::string, RVector > exportDataMap_;

    mutable KDTreeWrapper * tree_;

    bool oldTet10NumberingStyle_;

    std::string commentString_;

}; // class Mesh

} // namespace GIMLI;

#endif // _GIMLI_MESH__H

