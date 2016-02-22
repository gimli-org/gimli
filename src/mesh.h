/***************************************************************************
 *   Copyright (C) 2006-2016 by the resistivity.net development team       *
 *   Carsten Rücker carsten@gimli.org                                      *
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
    BoundingBox(const Pos < ValueType > & min=Pos < double >(0, 0, 0),
                const Pos < ValueType > & max=Pos < double >(1.0, 1.0, 1.0))
        : min_(min), max_(max){
    }

    /*! Construct BBox from position vector */
    BoundingBox(const std::vector < Pos < ValueType > > & vPos){
        min_ = Pos < ValueType >((ValueType)MAX_DOUBLE, (ValueType)MAX_DOUBLE, (ValueType)MAX_DOUBLE);
        max_ = min_ * -1.0;
        for (uint i = 0; i < vPos.size(); i ++){
            min_[0] = std::min(vPos[i][0], min_[0]);
            min_[1] = std::min(vPos[i][1], min_[1]);
            min_[2] = std::min(vPos[i][2], min_[2]);
            max_[0] = std::max(vPos[i][0], max_[0]);
            max_[1] = std::max(vPos[i][1], max_[1]);
            max_[2] = std::max(vPos[i][2], max_[2]);
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

    /*! Check if a point lie inside (with boundary). */
    bool isInside(const Pos < ValueType > & p){
        return ((p[0] <= max_[0] && p[0] >= min_[0]) && 
                (p[1] <= max_[1] && p[1] >= min_[1]) && 
                (p[2] <= max_[2] && p[2] >= min_[2]));
    }
    
    /*! Set minimum Position. */
    void setMin(const Pos < ValueType > & min) { min_ = min; }
    /*! Return a copy of the minimum position. */
    const Pos < ValueType > & min() const { return min_; }

    /*! Set maximum Position. */
    void setMax(const Pos < ValueType > & max) { max_ = max; }
    /*! Return a copy of the maximum position. */
    const Pos < ValueType > & max() const { return max_; }

protected:

    void copy_(const BoundingBox < ValueType > & bbox){
    	min_ = bbox.min();
        max_ = bbox.max();
    }

    Pos < ValueType > min_;
    Pos < ValueType > max_;
};

template < class ValueType > std::ostream & operator << (std::ostream & str, const BoundingBox< ValueType > & bb){
    str << "BoundingBox [" << bb.min() << ", " << bb.max() << "]" ;
    return str;
}

DLLEXPORT std::ostream & operator << (std::ostream & str, const Mesh & mesh);


class DLLEXPORT RegionMarkerPLC : public RVector3{
public:
    RegionMarkerPLC(const RVector3 & pos, int marker, double area=0.0)
    : RVector3(pos), marker_(marker), area_(area){}
    
    ~RegionMarkerPLC(){}
    
    inline int marker() const {return marker_;}
    inline double area() const {return area_;}
    
protected:
    int marker_;
    double area_;
};
    

class DLLEXPORT Mesh {

public:
    typedef RegionMarkerPLC RegionMarker;
    typedef std::vector< RegionMarker > RegionMarkerList;
    typedef RVector3 HoleMarker;
    typedef std::vector< RVector3 > HoleMarkerList;
    
    /*! Default constructor, create empty mesh with dimension dim */
    Mesh(Index dim=2);

    /*! Constructor, read mesh from filename */
    Mesh(const std::string & filename, bool createNeighbourInfos=true);

    /*! Copy constructor. */
    Mesh(const Mesh & mesh);

    /*! Copy assignment operator. */
    Mesh & operator = (const Mesh & mesh);

    /*! Default destructor. */
    ~Mesh();

    /*! Clear all data, inclusive all caches.*/
    void clear();

    /*!If the mesh is static in geometry and shape some useful informations are cached. 
     * (cell sizes, boundary sizes, ...)
     For dynamic meshes, i.e., node positions can be moved, you have to set staticGeometry to false to avoid any caching.*/
    void setStaticGeometry(bool stat);
    
    /*! Return true if this mesh have static geometry. [Default=True]*/
    inline bool staticGeometry() const { return staticGeometry_; }
        
    /*! Set the dimension of the mesh. [Default = 2] */
    void setDimension(uint dim){ dimension_ = dim;}
    
    /*! Return the dimension of this mesh.*/
    uint dimension() const { return dimension_; }
    
    /*! Shortcut for \ref dimension.*/
    uint dim() const { return dimension_; }

    //** start creation stuff
    Node * createNode(double x, double y, double z, int marker=0);

    Node * createNode(const Node & node);

    Node * createNode(const RVector3 & pos, int marker=0);

    Node * createNodeWithCheck(const RVector3 & pos, double tol=1e-6);

    Boundary * createBoundary(std::vector < Node * > & nodes, int marker=0);
    /*! Create a boundary from the given node indieces */
    Boundary * createBoundary(const IndexArray & nodes, int marker=0);
    Boundary * createBoundary(const Boundary & bound);
    Boundary * createBoundary(const Cell & cell);
    Boundary * createNodeBoundary(Node & n1, int marker = 0);
    Boundary * createEdge(Node & n1, Node & n2, int marker = 0);
    Boundary * createEdge3(Node & n1, Node & n2, Node & n3, int marker = 0);
    Boundary * createTriangleFace(Node & n1, Node & n2, Node & n3, int marker = 0);
    Boundary * createQuadrangleFace(Node & n1, Node & n2, Node & n3, Node & n4, int marker = 0);

    /*! Create empty cell without a node or a shape. */
    Cell * createCell(int marker = 0);
    Cell * createCell(std::vector < Node * > & nodes, int marker=0);
    /*! Create a cell from the given node indieces */
    Cell * createCell(const IndexArray & nodes, int marker=0);
    Cell * createCell(const Cell & cell);
    Cell * createTriangle(Node & n1, Node & n2, Node & n3, int marker=0);
    Cell * createQuadrangle(Node & n1, Node & n2, Node & n3, Node & n4, int marker=0);
    Cell * createTetrahedron(Node & n1, Node & n2, Node & n3, Node & n4, int marker=0);

    /*! Create a cell as a copy of a cell from an alternative mesh. Each Node of cell will be created with check. */
    Cell * copyCell(const Cell & cell);
    /*! Create a boundary as a copy of a boundary from an alternative mesh. Each Node of cell will be created with check. */
    Boundary * copyBoundary(const Boundary & bound);

    /*! Delete all given cells from the given mesh. Warning will be really deleted.*/
    void deleteCells(const std::vector < Cell * > & cells);
    
    void create1DGrid(const RVector & x);

    void create2DGrid(const RVector & x, const RVector & y, int markerType=0);

    void create3DGrid(const RVector & x, const RVector & y, const RVector & z,
                      int markerType=0);

    /*! Create one dimensional grid. 
     * Boundary on the domain border will get 
     * marker: 1,2 for: left, right.*/
    void createGrid(const RVector & x) { create1DGrid(x); }

    /*! Create two dimensional grid. Boundary on the domain border will get the
     * marker: 1,2,3,4 for: left, right, top, bottom*/
    void createGrid(const RVector & x, const RVector & y, int markerType=0) {
        create2DGrid(x, y, markerType);
    }
    /*! Create three dimensional grid. Boundary on the domain border will get the
     * marker: 1,2,3,4,5,6 for: left, right, top, bottom, front, back*/
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

    /*! Create a new mesh that is a part from this mesh, based on cell-ids */
    Mesh createMeshByCellIdx(const IndexArray & idxList);
    
    /*! Create a partly mesh from mesh, based on cell-ids */
    void createMeshByCellIdx(const Mesh & mesh, const IndexArray & idxList);

    /*! Create a partly mesh without cells from mesh, based on a vector of ptrs to boundaries */
    void createMeshByBoundaries(const Mesh & mesh, const std::vector < Boundary * > & bounds);

    /*! Create a partly mesh from mesh, based on meshs attributes. For a single attribute set to to 0, for unlimited set to to -1 */
    void createMeshByMarker(const Mesh & mesh, int from, int to=-1);
    //** end creation stuff

    //! Show some infos
    void showInfos(){ std::cout << *this << std::endl; }

    const std::vector< Node * > & nodes() const { return nodeVector_; }

    const std::vector< Cell * > & cells() const { return cellVector_; }
    
    /*! Return const reference to all boundaries*/
    const std::vector< Boundary * > & boundaries() const { return boundaryVector_; }

    /*! Return vector of node from index list */
    std::vector< Node * > nodes(const IndexArray & ids) const;
    
    /*! Return vector of cells from index list */
    std::vector< Cell * > cells(const IndexArray & ids) const;
    
    /*! Return vector of boundaries from index list */
    std::vector< Boundary * > boundaries(const IndexArray & ids) const;
    
    
    
    inline Index nodeCount() const { return nodeVector_.size(); }
    Node & node(Index i) const;
    Node & node(Index i);

    Index cellCount() const { return cellVector_.size(); }
    Cell & cell(Index i) const;
    Cell & cell(Index i);

    Index boundaryCount() const { return boundaryVector_.size(); }
    Boundary & boundary(Index i) const;
    Boundary & boundary(Index i);

    /*! Return a vector of all node positions */
    R3Vector positions() const;
    
    /*! Return a vector of node positions for an index vector */
    R3Vector positions(const IndexArray & idx) const;
    
    /*! Return all node positions. */
    R3Vector nodeCenters() const;

    /*! Return a vector of all cell center positions*/
    R3Vector cellCenters() const;
    R3Vector cellCenter() const { return cellCenters(); }

    /*! Return a vector of all center positions for all boundaries */
    R3Vector boundaryCenters() const;
        
    /*! Return the reference to a RVector of all cell sizes. Cached for static geometry.*/
    RVector & cellSizes() const;

    /*! Return the reference to a RVector of all boundary sizes. Cached for static geometry. */
    RVector & boundarySizes() const;

    /*! Return the reference to the vector of scaled normal directions for each boundary.
     * Cached for static geometry and will be build on first call. Not thread safe, perhaps not python GC safe. 
     Return \f$ \{A_i \vec{n}_i\} \forall i = [0..N_B]\f$.
     Where \f$ A_i\f$ is the size and \f$ \vec{n}_i\f$ the normal direction for the i-th boundary. 
     If you want to use this, i.e. for the calculation of inside or outside flow through the boundary, you need to recognize the orientation of this boundary to the cell the flow goes into or comes from.
     For the left cell neighbor the normal direction should be always the outer normal.*/
    R3Vector & boundarySizedNormals() const;
    
    
    /*! Return a vector of all boundary marker */
    IVector boundaryMarker() const;

    /*! Return a vector of all node marker */
    IVector nodeMarker() const;

    /*! Return an index vector of all nodes that match the marker */
    IndexArray findNodesIdxByMarker(int marker) const;

//     /*! Return an index list of all nodes that match the marker */
//     std::list < Index > findListNodesIdxByMarker(int marker) const;

    /*! Return a vector of boundary ptrs with the boundary marker equal marker.*/
    std::vector < Boundary * > findBoundaryByMarker(int marker) const;

    /*! Return a vector of boundary ptrs with the boundary marker between [from and to). \n
        for to equal open end set to = MAX_INT */
    std::vector < Boundary * > findBoundaryByMarker(int from, int to) const;
    
    /*! Set the marker to all boundaries in index array. */
    void setBoundaryMarker(const IndexArray & ids, int marker);
    
    /*! Return ptr to the cell that match position pos, counter holds amount of touch tests.
        Searching is done first by nearest-neighbour-kd-tree search,
        followed by slope-search if extensive is set. Return NULL if no cell can be found. */
    Cell * findCell(const RVector3 & pos, size_t & counter, bool extensive) const ;

    /*! Shortcut for \ref findCell(const RVector3 & pos, size_t & counter, bool extensive) */
    Cell * findCell(const RVector3 & pos, bool extensive=true) const {
        size_t counter; return findCell(pos, counter, extensive); }

    /*! Return the index to the node of this mesh with the smallest distance to pos. */
    Index findNearestNode(const RVector3 & pos);
    
    /*! Return vector of cell ptrs with marker match the range [from .. to). \n
        For single marker match to is set to 0, for open end set to = -1 */
    std::vector < Cell * > findCellByMarker(int from, int to=0) const;

    /*! Return vector of cell ptrs with attribute match the range [from .. to). \n
        For single attribute match to is set to 0.0, for open end set to = -1.0 */
    std::vector < Cell * > findCellByAttribute(double from, double to=0.0) const;

    //** end get infos stuff

    //** start mesh modification stuff
    /*! DEPRECATED Prolongate the attribute of each cell in emptyList by the attribute 
     * from neighbouring cells.
     * The attributes have to be lower than \ref TOLERANCE.
     * This function is called recursively until all zero-attribute-cells in
     * emptyList are filled with an attribute greater than Zero. */
    void fillEmptyCells(const std::vector< Cell * > & emptyList, 
                        double background=-1.0);

    /*! Prolongate the empty (lower than \ref TOLERANCE.) cell values in vals 
     * from its neighbouring cells.
     * This function is called recursively until all zero-attribute-values in
     * vals are filled with an attribute greater than Zero. 
     * RVector vals need to be of size \ref cellCount().
     * If Background is unequal -1.0 all empty values will be set to background.
     */
    void prolongateEmptyCellsValues(RVector & vals, double background=-1.0) const;
                            
    void recountNodes();

    void sortNodes(const IndexArray & perm);

    /*! Return true if createNeighbourInfos is called once */
    inline bool neighboursKnown() const { return neighboursKnown_; }
    
    /*! Remove from each boundary the ptr to the corresponding left and right cell*/
    void cleanNeighbourInfos();

    /*! Search and set to each boundary the corresponding left and right cell.*/
    void createNeighbourInfos(bool force=false);

    /*! Create and store boundaries and neighboring information for this cell.*/
    void createNeighbourInfosCell_(Cell *c);
        
    void relax();

    /*! Smooth the mesh via moving all free nodes into the average of all neighboring nodes. Repeat this smoothIteration times. There is currently only this smoothFunction. EdgeSwapping is deactivated.*/
    void smooth(bool nodeMoving, bool edgeSwapping, uint smoothFunction, uint smoothIteration);

    /*! Scale the mesh with \ref RVector3 s. And return a reference to the mesh (no copy)*/
    Mesh & scale(const RVector3 & s);

    /*! Translate the mesh with \ref RVector3 t. And return a reference to the mesh (no copy)*/
    Mesh & translate(const RVector3 & t);

    /*! Rotate mesh the with \ref RVector3 r, r in radian, If you want to rotate in degree, use \ref degToRad(const RVector3 & deg). 
     And return a reference to the mesh (no copy) */
    Mesh & rotate(const RVector3 & r);
    //** end mesh modification stuff

    /*! Swap coordinate i with j for i and j lower then dimension of the mesh*/
    void swapCoordinates(Index i, Index j);
    //** end mesh modification stuff
    
    /*! apply a 4x4 transformation matrix to the whole mesh*/
    template < class Matrix > Mesh & transform(const Matrix & mat){
//         std::for_each(nodeVector_.begin(), nodeVector_.end(),
//                        bind2nd(std::mem_fun(&Node::pos().transform), mat));
        for (uint i = 0; i < nodeVector_.size(); i ++) nodeVector_[i]->pos().transform(mat);
        rangesKnown_ = false;
        return *this;
    }

    //** start I/O stuff
    int save(const std::string & fileName, IOFormat format = Binary) const;
    int saveAscii(const std::string & fileName) const;

    /*! Be carefull with interchanging binary meshs between 32-64bit architecture. 
     * Atm we save fixed int for counter and idx.
     * We have to replace and test it with uint32 or uint16 */
    int saveBinary(const std::string & fileName) const;

    /*! Load Mesh from file and try to import fileformat regarding file suffix.
     * If createNeighbourInfos if set, the mesh is checked for consistency and 
     * missing boundaries will be created. */
    void load(const std::string & fileName, 
              bool createNeighbours=true, IOFormat format=Binary);

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
    int exportMidCellValue(const std::string & fileName, const RVector & data1, const RVector & data2) const ;
    // no default arg here .. pygimli@win64 linker bug
    int exportMidCellValue(const std::string & fileName, const RVector & data1) const{
        RVector data2(0);
        return exportMidCellValue(fileName, data1, data2);
    }
    
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

    /*!DEPRECATED will be removed. 
        Add data to the mesh that will be saved with by using the binary mesh 
     * format v.2. or exported with the appropriate name in VTK format,
     * if the size of data equals the amount of nodes, cells or boundaries. 
     */
    void addExportData(const std::string & name, const RVector & data);

    /*! Return the data with a given name. 
     * If there is no such data an exception is thrown.*/
    RVector exportData(const std::string & name) const;

    /*! Return the full data map read only. */
    const std::map< std::string, RVector > & exportDataMap() const {
        return exportDataMap_; }

    /*! Set the full data map.*/
    void setExportDataMap(const std::map< std::string, RVector > & eMap) { 
        exportDataMap_ = eMap; }

    /*! Empty the data map.*/
    void clearExportData();

    /*!*/
    void addData(const std::string & name, const RVector & data){ addExportData(name, data); }
    /*!*/
    RVector data(const std::string & name) const { return exportData(name); }
    /*!*/
    bool haveData(const std::string & name) const { return exportDataMap_.count(name) > 0; }
    /*!*/
    const std::map< std::string, RVector > & dataMap() const { return exportDataMap_; }
  
    /*! Print data map info.*/
    void dataInfo() const;
  
    /*! Set the comment for VTK Ascii export headline.*/
    void setCommentString(const std::string & commentString) {commentString_ = commentString;}

    /*! Return comment for VTK Ascii export headline.*/
    const std::string & commentString() const {return commentString_;}

    void mapCellAttributes(const std::map < float, float > & aMap);
    
    void mapAttributeToParameter(const IndexArray & cellMapIndex,
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

    /*! Set the cell marker of all indices in ids to marker. */
    void setCellMarkers(const IndexArray & ids, int marker);

    /*! Set all cell marker the values in attribute. */
    void setCellMarkers(const IVector & marker);
    
    /*! Set all cell marker the values in attribute (casting to int)*/
    void setCellMarkers(const RVector & attribute);
    
    /*! Return a vector of all cell marker */
    IVector cellMarkers() const;
    

    //** probably deprecated
    double xmin() const { findRange_(); return minRange_[0]; }
    double ymin() const { findRange_(); return minRange_[1]; }
    double zmin() const { findRange_(); return minRange_[2]; }
    double xmax() const { findRange_(); return maxRange_[0]; }
    double ymax() const { findRange_(); return maxRange_[1]; }
    double zmax() const { findRange_(); return maxRange_[2]; }

    const RBoundingBox boundingBox() const { findRange_(); return RBoundingBox(minRange_, maxRange_);}


    /*! Return the reference to the matrix for cell value to boundary value interpolation matrix. */
    RSparseMapMatrix & cellToBoundaryInterpolation() const;

    /*!Return the divergence for each cell of a given vector field for each 
     * boundary.
     * The divergence is calculated by simple 1 point boundary integration 
     * over each cell. 
     * \f$ d(cell) = \sum_boundaries V(boundary center) \cdot n(boundary)\f$
     * Higher order integration needs to be implemented. 
     * Contact the author if you need this.*/
    RVector divergence(const R3Vector & V) const;

    /*! Interpolate boundary based values to cell based gradients. */
    R3Vector boundaryDataToCellGradient(const RVector & boundaryData) const;
    
    /*! Interpolate cell based values to boundary based gradients. */
    R3Vector cellDataToBoundaryGradient(const RVector & cellData) const;
        
    /*! Interpolate cell based values to boundary based gradients with a given cell Gradient.*/
    R3Vector cellDataToBoundaryGradient(const RVector & cellData,
                                        const R3Vector & cellGradient) const;

    /*! Add a region marker for tetgen or triangle creation if the mesh 
     *is a PLC, if area is < 0 a hole is added. */    
    void addRegionMarker(const RVector3 & pos, int marker, double area=0);
    void addRegionMarker(const RegionMarker & reg);

    const RegionMarkerList & regionMarker() const { return regionMarker_; }
    
    /*! Add a hole marker for tetgen or triangle creation if the mesh
     * is a PLC */    
    void addHoleMarker(const RVector3 & pos);
    
    /*!Return read only reference for all defined hole regions. */
    const HoleMarkerList & holeMarker() const { return holeMarker_; }
    
    
protected:
    void copy_(const Mesh & mesh);

    void findRange_() const ;

    Node * createNode_(const RVector3 & pos, int marker, int id);

    template < class B > Boundary * createBoundary_(
        std::vector < Node * > & nodes, int marker, int id){
        
        if (id == -1) id = boundaryCount();
        boundaryVector_.push_back(new B(nodes));
        boundaryVector_.back()->setMarker(marker);
        boundaryVector_.back()->setId(id);
        return boundaryVector_.back();
    }

    template < class B > Boundary * createBoundaryChecked_(
        std::vector < Node * > & nodes, int marker){
        
        Boundary * b = findBoundary(nodes);
        if (!b) {
            b = createBoundary_< B >(nodes, marker, boundaryCount());
        } else {
            if (marker != 0) b->setMarker(marker);
        }
        return b;
    }

    template < class C > Cell * createCell_(
        std::vector < Node * > & nodes, int marker, int id){
        
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

    mutable RVector3 minRange_;
    mutable RVector3 maxRange_;
    mutable bool rangesKnown_;

    bool neighboursKnown_;
    
    mutable KDTreeWrapper * tree_;

    /*! A static geometry mesh caches geometry informations. */
    bool staticGeometry_;
    mutable RVector cellSizesCache_;
    mutable RVector boundarySizesCache_;
    mutable R3Vector boundarySizedNormCache_;
    
    mutable RSparseMapMatrix * cellToBoundaryInterpolationCache_;
    
    bool oldTet10NumberingStyle_;

    std::map< std::string, RVector > exportDataMap_;
    std::string commentString_;

    // for PLC creation
    RegionMarkerList regionMarker_;
    HoleMarkerList holeMarker_;
    
}; // class Mesh

} // namespace GIMLI;

#endif // _GIMLI_MESH__H

