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

#include "mesh.h"
#include "node.h"
#include "matrix.h"
#include "vectortemplates.h"

#include <map>
#include <fstream>

namespace GIMLI{

void Mesh::load( const std::string & fbody, IOFormat format ){

    if ( fbody.find( ".mod" ) != std::string::npos ){
        importMod( fbody );
    } else if ( fbody.find( ".vtk" ) != std::string::npos ){
        importVTK( fbody );
    } else if ( fbody.find( ".vtu" ) != std::string::npos ){
        importVTU( fbody );
    } else if ( format == Binary || fbody.find( MESHBINSUFFIX ) != std::string::npos ) {
        try {
             return loadBinary( fbody );
        } catch( std::exception & e ){
            std::cout << "Failed to loadBinary " << e.what() << std::endl;
            std::cout << "try load bms.v2" << std::endl;
            loadBinaryV2( fbody );
        }
    } else {
        loadAscii( fbody );
    }
}

void Mesh::loadAscii( const std::string & fbody ){
    THROW_TO_IMPL
}

int Mesh::save( const std::string & fbody, IOFormat format ) const {
  if ( format == Binary || fbody.find( MESHBINSUFFIX ) != std::string::npos ) {
    return saveBinary( fbody );
  } else return saveAscii( fbody );
  return 1;
}

int Mesh::saveAscii( const std::string & fbody ) const {
  //** save nodes
  std::fstream file; if ( !openOutFile( fbody + ".n", & file ) ) {  return 0; }

  file.precision( 14 );
  for ( uint i = 0; i < nodeCount(); i ++ ){
    for ( uint j = 0; j < dim(); j ++ ) file << node( i ).pos()[ j ] << "\t";
    file << node( i ).marker() << std::endl;
  }
  file.close();

  //** save cells;
  if ( !openOutFile( fbody + ".e", & file ) ) { return 0; }

  for ( int i = 0, imax = cellCount(); i < imax; i++){
    for ( int j = 0, jmax = cell( i ).nodeCount(); j < jmax; j++){
      file << cell( i ).node( j ).id() << "\t" ;
    }
    file << cell( i ).marker() << std::endl;
  }
  file.close();

  //** save boundarys;
  if ( !openOutFile( fbody + ".s", & file ) ) { return 0; }

  for ( int i = 0, imax = boundaryCount(); i < imax; i++){
    for ( int j = 0, jmax = boundary( i ).nodeCount(); j < jmax; j++){
      file << boundary( i ).node( j ).id() << "\t" ;
    }
    file << "-33\t-33\t" << boundary( i ).marker() << std::endl;
  }

  file.close();
  return 1;
}

int Mesh::saveBinary( const std::string & fbody ) const {
    std::string fileName( fbody.substr( 0, fbody.rfind( MESHBINSUFFIX ) ) + MESHBINSUFFIX );

    FILE *file;
    file = fopen( fileName.c_str(), "w+b" );

    int dimension = dim();
    if ( !fwrite( &dimension, sizeof( int ), 1, file ) ){
    }

    //! write vertex dummy-infos;
    int dummy[ 127 ]; memset( dummy, 0, 127 * sizeof(int) );
    if ( !fwrite( dummy, sizeof(int), 127, file ) ) {
    }

    int nVerts = nodeCount();
    if ( !fwrite( &nVerts, sizeof(int), 1, file ) ){
    }

    double * koord = new double[ dimension * nVerts ];
    for ( int i = 0; i < nVerts; i ++ ){
        for ( int j = 0; j < dimension; j ++ ){
            koord[ i * dimension + j ] = node( i ).pos()[ j ];
        }
    }

    if ( !fwrite( koord, sizeof(double), dimension * nVerts, file ) ){

    }
    delete [] koord;

    int * marker = new int[ nVerts ];
    for ( int i = 0; i < nVerts; i ++ ) marker[ i ] = node( i ).marker();

    if ( !fwrite( marker, sizeof( int ), nVerts, file ) ){

    }
    delete [] marker;

    //! write cell dummy-infos
    if ( !fwrite( dummy, sizeof( int ), 127, file ) ){

    }

    int nCells = cellCount();
    if ( !fwrite( &nCells, sizeof( int ), 1, file ) ){

    }

    int * cellVerts = new int[ nCells ];
    int nCellIdx = 0;

    for ( int i = 0; i < nCells; i ++ ){
        cellVerts[ i ] = cell( i ).nodeCount();
        nCellIdx += cellVerts[ i ];
    }
    if ( !fwrite( cellVerts, sizeof(int), nCells, file ) ){

    }

    int * idx = new int[ nCellIdx ];

    int count = 0;
    for ( int i = 0; i < nCells; i ++ ){
        for ( int j = 0; j < cellVerts[ i ]; j ++ ){
            idx[ count ] = cell( i ).node( j ).id();
            count ++;
        }
    }

    if ( !fwrite( idx, sizeof(int), nCellIdx, file ) ){

    }
    delete [] idx;
    delete [] cellVerts;

    double * attribute = new double[ nCells ];
    for ( int i = 0; i < nCells; i ++ ) attribute[ i ] = cell( i ).marker();

    if ( !fwrite( attribute, sizeof( double ), nCells, file ) ){
    }

    delete [] attribute;

  //** write boundary dummy-infos
    if ( !fwrite( dummy, sizeof(int), 127, file ) ){

    }

    int nBounds = boundaryCount();
    if ( !fwrite( &nBounds, sizeof(int), 1, file ) ){

    }

    int * boundVerts = new int[ nBounds ];
    int nBoundIdx = 0;
    for ( int i = 0; i < nBounds; i ++ ){
        boundVerts[ i ] = boundary( i ).nodeCount();
        nBoundIdx += boundVerts[ i ];
    }
    if ( !fwrite( boundVerts, sizeof(int), nBounds, file ) ){

    }

    idx = new int[ nBoundIdx ];
    count = 0;
    for ( int i = 0; i < nBounds; i ++ ){
        for ( int j = 0; j < boundVerts[ i ]; j ++ ){
            idx[ count ] = boundary( i ).node( j ).id();
            count ++;
        }
    }

    if ( !fwrite( idx, sizeof(int), nBoundIdx, file ) ){

    }

  delete [] idx;
  delete [] boundVerts;

  marker = new int[ nBounds ];
  for ( int i = 0; i < nBounds; i ++ ) marker[ i ] = boundary( i ).marker();
 if (!fwrite( marker, sizeof( int ), nBounds, file )){}
  delete [] marker;

  idx = new int[ nBounds ];
  for ( int i = 0; i < nBounds; i ++ ){
    if ( boundary( i ).leftCell() != NULL ){
      idx[ i ] = boundary( i ).leftCell()->id();
    } else {
      idx[ i ] = -1;
    }
  }
  if ( !fwrite( idx, sizeof( int ), nBounds, file ) ){}

  for ( int i = 0; i < nBounds; i ++ ){
    if ( boundary( i ).rightCell() != NULL ){
      idx[ i ] = boundary( i ).rightCell()->id();
    } else {
      idx[ i ] = -1;
    }
  }
  if ( !fwrite( idx, sizeof( int ), nBounds, file )) {}
  delete [] idx;

  fclose( file );
  return 1;
}

void Mesh::loadBinary( const std::string & fbody ){
//   sizeof( int ) = 4 byte
//   int[ 1 ] dimension
//   int[ 127 ] dummy vertices information
//   int[ 1 ] nVerts, number of vertices
//   double[ dimension * nVerts ]; coordinates,  dimension == 2 ( x, y ), dimension == 3 (x, y, z )
//   int[ nVerts ] vertex markers
//   int[ 127 ] dummy cell information
//   int[ 1 ] nCells, number of cell
//   int[ nCells ] cellVerts; number of nodes for each cell
//   int[ sum( cellVerts ) ] cellsidx
//   double[ nCells ] attribute, cell attributes
//   int[ 127 ] dummy boundary information
//   int[ 1 ] nBounds, number of boundarys
//   int[ nBounds ] boundVerts; numbers of nodes for each boundary
//   int[ sum( boundVerts ) ] boundIdx
//   int[ nBounds ] boundary markers
//   int[ nBounds ] leftNeighbor idx (-1) if no neighbor present or info unavailable
//   int[ nBounds ] rightNeighbor idx (-1) if no neighbor present or info unavailable
    //std::cout << "load binary " << fbody << std::endl;
    clear();
    std::string fileName( fbody.substr( 0, fbody.rfind( MESHBINSUFFIX ) ) + MESHBINSUFFIX );

    FILE *file; file = fopen( fileName.c_str(), "r+b" );
    if ( !file ) {
        throwError( EXIT_OPEN_FILE, WHERE_AM_I + " " + fileName + ": " + strerror( errno ) );
    }

    int dim = 0;
    uint ret = fread( &dim, sizeof( int ), 1, file );

    if ( ( dim !=2 && dim !=3 ) || ( ret == 0 ) ){
        throwError( 1, WHERE_AM_I + " cannot determine dimension " + toStr( dim ) );
    }
    this->setDimension( dim );
    //** read vertex dummy-infos
    int dummy[ 127 ]; ret = fread( dummy, sizeof( int ), 127, file );
    int nVerts; ret = fread( &nVerts, sizeof( int ), 1, file );
    double * koords = new double[ dimension_ * nVerts ];
    ret = fread( koords, sizeof( double ), dimension_ * nVerts, file );

    int * nodeMarker = new int[ nVerts ]; ret = fread( nodeMarker, sizeof( int ), nVerts, file );
    //** read cell dummy-infos
    ret = fread( dummy, sizeof( int ), 127, file );

  int nCells; ret = fread( &nCells, sizeof( int ), 1, file );
  int * cellVerts = new int[ nCells ]; ret = fread( cellVerts, sizeof( int ), nCells, file );
  int nCellIdx = 0; for ( int i = 0; i < nCells; i ++ ) nCellIdx += cellVerts[ i ];
  int * cellIdx = new int[ nCellIdx ]; ret = fread( cellIdx, sizeof( int ), nCellIdx, file );
  double * attribute = new double[ nCells ]; ret = fread( attribute, sizeof( double ), nCells, file );

  //** read boundary dummy-infos
  ret = fread( dummy, sizeof( int ), 127, file );

  int nBounds; ret = fread( &nBounds, sizeof(int), 1, file );
  int * boundVerts = new int[ nBounds ]; ret = fread( boundVerts, sizeof( int ), nBounds, file );
  int nBoundIdx = 0; for ( int i = 0; i < nBounds; i ++ ) nBoundIdx += boundVerts[ i ];
  int * boundIdx = new int[ nBoundIdx ]; ret = fread( boundIdx, sizeof( int ), nBoundIdx, file );
  int * boundMarker = new int[ nBounds ]; ret = fread( boundMarker, sizeof( int ), nBounds, file );
  int * left = new int[ nBounds ]; ret = fread( left, sizeof( int ), nBounds, file );
  int * right = new int[ nBounds ]; ret = fread( right, sizeof( int ), nBounds, file );

    //** create Nodes;
    nodeVector_.reserve( nVerts );
    for ( int i = 0; i < nVerts; i ++ ){
        createNode( RVector3( 0.0, 0.0, 0.0 ) );
        for ( uint j = 0; j < dimension_; j ++ ){
            node( i ).pos()[ j ] = koords[ i * dimension_ + j ];
        }
        node( i ).setMarker( nodeMarker[ i ] );
    }

    //** create Cells;
    int count = 0;
    std::vector < Node * > pNodeVector;

    cellVector_.reserve( nCells );
    for ( int i = 0; i < nCells; i ++ ){
        std::vector < Node * > nodes( cellVerts[ i ] );
        for ( uint j = 0; j < nodes.size(); j ++ ) {
            nodes[ j ] = & node( cellIdx[ count + j ] );
        }
        createCell( nodes );
        count += cellVerts[ i ];
    }

    for ( int i = 0; i < nCells; i ++ ) {
        cell( i ).setMarker( (int)rint(attribute[ i ]) );
        cell( i ).setAttribute( attribute[ i ] );
    }

    //** create Boundaries;
    count = 0;
    boundaryVector_.reserve( nBounds );
    for ( int i = 0; i < nBounds; i ++ ){
        std::vector < Node * > nodes( boundVerts[ i ] );
        for ( uint j = 0; j < nodes.size(); j ++ ) nodes[ j ] = & node( boundIdx[ count + j ] );
        createBoundary( nodes );
        count += boundVerts[ i ];
    }

    for ( int i = 0; i < nBounds; i ++ ){
        boundary( i ).setMarker( boundMarker[ i ] );
        if ( left[ i ] != -1 ) {
            boundary( i ).setLeftCell( &cell( left[ i ] ) );
        }
        if ( right[ i ] != -1 ) boundary( i ).setRightCell( &cell( right[ i ] ) );
    }

    delete [] koords;
    delete [] nodeMarker;

    delete [] cellIdx;
    delete [] cellVerts;
    delete [] attribute;

    delete [] boundIdx;
    delete [] boundVerts;
    delete [] boundMarker;
    delete [] left;
    delete [] right;

    fclose( file );

}

template < class ValueType > void writeToFile(FILE * file, const ValueType & v, int count=1){
    if (!fwrite(&v, sizeof(ValueType), count, file)){
        throwError(EXIT_OPEN_FILE, WHERE_AM_I + strerror(errno)  + str(errno));
    }
}

template < class ValueType > void readFromFile(FILE * file, ValueType & v, int count=1){
    uint ret = fread(&v, sizeof( ValueType ), count, file);

    if (ret && ferror(file)){
        throwError(EXIT_OPEN_FILE, WHERE_AM_I + strerror(errno) + " " + str(errno));
    }
}

void Mesh::saveBinaryV2( const std::string & fbody ) const {
//     std::cout << sizeof( uint8 ) << " "  << sizeof( uint16 ) << " " << sizeof( uint32 ) << " " << sizeof( uint64 ) << std::endl;

//   uint8[ 1 ] dimension
//   uint8[ 1 ] file format version
//   uint32[ 1 ] nVerts, number of vertices, max 2 ^ 32 (4e9)
//   double[ 3 * nVerts ]; coordinates, dimension == 3 (x, y, z )
//   int32[ nVerts ] vertex markers [-2e9, .. , 2e9 ]
//   uint32[ 1 ] nCells, number of cell
//   uint8[ nCells ] cell nodeCount;
//   uint32[ sum( cell nodeCount) ] cellsidx

//   int32[ nCells ] cellMarkers [-2e9, .. , 2e9 ]
//   uint32[ 1 ] nBounds, number of boundarys
//   uint8[ nBounds ] bound nodeCount
//   uint32[ sum( bound nodeCount ) ] boundsidx
//   int32[ nBounds ] boundaryMarkers [-2e9, .. , 2e9 ]
//   int32[ nBounds ] leftNeighbour idx (-1) if no neighbor present or info unavailable
//   int32[ nBounds ] rightNeighbour idx (-1) if no neighbor present or info unavailable

    std::string fileName( fbody.substr( 0, fbody.rfind( MESHBINSUFFIX ) ) + MESHBINSUFFIX );

    FILE *file;
    file = fopen( fileName.c_str(), "w+b" );
    if ( !file ) {
        throwError( EXIT_OPEN_FILE, WHERE_AM_I + " " + fileName + ": " + strerror( errno ) );
    }

    //** write preample
    writeToFile( file, uint8( this->dimension() ) );
    writeToFile( file, uint8( 2 ) );

    //** write nodes
    double * coord = new double[ 3 * this->nodeCount() ];
    for ( uint i = 0; i < this->nodeCount(); i ++ ){
        for ( uint j = 0; j < 3; j ++ ){
            coord[ i * 3 + j ] = node( i ).pos()[ j ];
        }
    }

    int32 * marker = new int32[ this->nodeCount() ];
    for ( uint i = 0; i < this->nodeCount(); i ++ ) marker[ i ] = node( i ).marker();

    writeToFile( file, uint32( this->nodeCount() ) );
    writeToFile( file, coord[ 0 ], 3 * this->nodeCount() );
    writeToFile( file, marker[ 0 ], this->nodeCount() );

    //** write Cells
    uint32 nCells = this->cellCount();
    uint8 * cellVerts = new uint8[ nCells ];
    uint32 nCellIdx = 0;
    for ( uint i = 0; i < nCells; i ++ ){
        cellVerts[ i ] = (uint8)cell( i ).nodeCount();
        nCellIdx += cellVerts[ i ];
    }
    uint32 * cellIdx = new uint32[ nCellIdx ];
    uint count = 0;
    for ( uint i = 0; i < nCells; i ++ ){
        for ( uint j = 0; j < cellVerts[ i ]; j ++ ){
            cellIdx[ count ] = cell( i ).node( j ).id();
            count ++;
        }
    }
    int32 * cellMarker = new int32[ nCells ];
    for ( uint i = 0; i < nCells; i ++ ) cellMarker[ i ] = cell( i ).marker();

    writeToFile( file, nCells );
    writeToFile( file, cellVerts[ 0 ], nCells );
    writeToFile( file, cellIdx[ 0 ], nCellIdx );
    writeToFile( file, cellMarker[ 0 ], nCells );

    //** write boundarys
    uint32 nBound = this->boundaryCount();
    uint8 * boundVerts = new uint8[ nBound ];
    std::vector < uint32 > boundIdx;
    int32 * boundMarker = new int32[ nBound ];
    int32 * leftCells = new int32[ nBound ];
    int32 * rightCells = new int32[ nBound ];

    count = 0;
    for ( uint i = 0; i < nBound;  i ++ ){
        Boundary * bound = &this->boundary( i );

        boundVerts[ i ] = (uint8)bound->nodeCount();

        for ( uint j = 0; j < boundVerts[ i ]; j ++ ){
            boundIdx.push_back( (uint32)bound->node( j ).id() );
        }

        boundMarker[ i ] = bound->marker();

        if ( bound->leftCell() ) {
            leftCells[ i ] = bound->leftCell()->id();
        } else leftCells[ i ] = -1;

        if ( bound->rightCell() ) {
            rightCells[ i ] = bound->rightCell()->id();
        } else rightCells[ i ] = -1;
    }


    writeToFile( file, nBound );
    writeToFile( file, boundVerts[ 0 ], nBound );
    writeToFile( file, boundIdx[ 0 ], boundIdx.size() );
    writeToFile( file, boundMarker[ 0 ], nBound );
    writeToFile( file, leftCells[ 0 ], nBound );
    writeToFile( file, rightCells[ 0 ], nBound );

    writeToFile( file, exportDataMap_.size() );

    for ( std::map < std::string, RVector >::const_iterator it = exportDataMap_.begin(); it != exportDataMap_.end(); it ++ ){
        writeToFile( file, it->first.length() );
        writeToFile( file, it->first[ 0 ], it->first.length() );
        writeToFile( file, it->second.size() );
        writeToFile( file, it->second[ 0 ], it->second.size() );
    }

    fclose( file );
    delete [] coord;
    delete [] marker;
    delete [] cellVerts;
    delete [] cellIdx;
    delete [] cellMarker;
    delete [] boundVerts;
    delete [] boundMarker;
    delete [] leftCells;
    delete [] rightCells;
}

void Mesh::loadBinaryV2( const std::string & fbody ) {

    this->clear();
    std::string fileName( fbody.substr( 0, fbody.rfind( MESHBINSUFFIX ) ) + MESHBINSUFFIX );

    FILE *file;
    file = fopen( fileName.c_str(), "r+b" );

    if ( !file ) {
        throwError( EXIT_OPEN_FILE, WHERE_AM_I + " " + fileName + ": " + strerror( errno ) );
    }

    uint8 dim; readFromFile( file, dim );
    if ( dim !=2 && dim !=3 ){
        throwError( 1, WHERE_AM_I + " cannot determine dimension " + toStr( dim ) );
    }
    this->setDimension( dim );
    uint8 version; readFromFile( file, version );

    //** read nodes
    uint32 nVerts; readFromFile( file, nVerts );

    if ( nVerts > 1e9 ){
        throwError( 1, WHERE_AM_I + " probably something wrong: nVerts > 1e9 " + toStr( nVerts ) );
    }

    double * coord = new double[ 3 * nVerts ]; readFromFile( file, coord[ 0 ], 3 * nVerts );
    int32 * marker = new int32[ nVerts ]; readFromFile( file, marker[ 0 ], nVerts );

    nodeVector_.reserve( nVerts );
    for ( uint i = 0; i < nVerts; i ++ ) {
        this->createNode( coord[ i * 3 ], coord[ i * 3 + 1 ], coord[ i * 3 + 2 ], marker[ i ] );
    }

    //** read cells
    uint32 nCells; readFromFile( file, nCells );
    uint8 * cellVerts = new uint8[ nCells ]; readFromFile( file, cellVerts[ 0 ], nCells );
    uint nCellIdx = 0; for ( uint i = 0; i < nCells; i ++ ) nCellIdx += cellVerts[ i ];
    uint32 * cellIdx = new uint32[ nCellIdx ];   readFromFile( file, cellIdx[ 0 ], nCellIdx );
    int32 * cellMarker = new int32[ nCells ]; readFromFile( file, cellMarker[ 0 ], nCells );

    //** create cells
    uint count = 0;
    cellVector_.reserve( nCells );
    for ( uint i = 0; i < nCells; i ++ ){
        std::vector < Node * > nodes( cellVerts[ i ] );
        for ( uint j = 0; j < nodes.size(); j ++ ) nodes[ j ] = & node( cellIdx[ count + j ] );
        this->createCell( nodes, cellMarker[ i ] );
        count += cellVerts[ i ];
    }

    //** read bounds
    uint32 nBound;                                  readFromFile( file, nBound );
    uint8 * boundVerts = new uint8[ nBound ];       readFromFile( file, boundVerts[ 0 ], nBound );
    uint nBoundIdx = 0; for ( uint i = 0; i < nBound; i ++ ) nBoundIdx += boundVerts[ i ];
    uint32 * boundIdx = new uint32[ nBoundIdx ];    readFromFile( file, boundIdx[ 0 ], nBoundIdx );
    int32 * boundMarker = new int32[ nBound ];      readFromFile( file, boundMarker[ 0 ], nBound );
    int32 * leftCells = new int32[ nBound ];        readFromFile( file, leftCells[ 0 ], nBound );
    int32 * rightCells = new int32[ nBound ];       readFromFile( file, rightCells[ 0 ], nBound );

    //** create cells
    count = 0;
    boundaryVector_.reserve( nBound );
    for ( uint i = 0; i < nBound; i ++ ){
        std::vector < Node * > nodes( boundVerts[ i ] );
        for ( uint j = 0; j < nodes.size(); j ++ ) nodes[ j ] = & node( boundIdx[ count + j ] );

        Boundary * bound = this->createBoundary( nodes, boundMarker[ i ] );
        count += boundVerts[ i ];

        if ( leftCells[ i ] > -1 ) bound->setLeftCell( &this->cell( leftCells[ i ] ) );
        if ( rightCells[ i ] > -1 ) bound->setRightCell( &this->cell( rightCells[ i ] ) );
    }

    size_t nData; readFromFile( file, nData );

    for ( uint i = 0; i < nData; i ++ ){
        size_t strLen; readFromFile( file, strLen );
        std::string str; str.resize( strLen ); readFromFile( file, str[0], strLen );
        size_t datLen; readFromFile( file, datLen );

        RVector dat( datLen ); readFromFile( file, dat[0], datLen );
        this->addExportData( str, dat );
        //delete [] str;
    }

    fclose( file );

    delete [] coord;
    delete [] marker;
    delete [] cellVerts;
    delete [] cellIdx;
    delete [] cellMarker;
    delete [] boundVerts;
    delete [] boundIdx;
    delete [] boundMarker;
    delete [] leftCells;
    delete [] rightCells;

}

int Mesh::exportSimple( const std::string & fbody, const RVector & data ) const {
  //output x y x y x y rhoa file
  std::fstream file; if ( !openOutFile( fbody , & file ) ){ exit( EXIT_MESH_EXPORT_FAILS ); }

  for ( uint i = 0; i < cellCount(); i ++ ){
    for ( uint j = 0; j < 3; j ++ ){
      file << cell( i ).node( j ).x() << "\t" << cell( i ).node( j ).y() << "\t";
    }
    file << data[ i ] << std::endl;
  }
  file.close();
  return 1;
}

void Mesh::exportVTK( const std::string & fbody  ) const {
    return exportVTK( fbody, exportDataMap_, std::vector < RVector3 >( 0 ) );
}

void Mesh::exportVTK( const std::string & fbody, const std::vector < RVector3 > & vec  ) const {
    return exportVTK( fbody, exportDataMap_, vec );
}

void Mesh::exportVTK( const std::string & fbody, const std::map< std::string, RVector > & dataMap ) const {
    return exportVTK( fbody, dataMap, std::vector < RVector3 >( 0 ) );
}
                        
void Mesh::exportVTK( const std::string & fbody, const std::map< std::string, RVector > & dataMap, 
                      const std::vector < RVector3 > & vec ) const {
    bool verbose = false;
    
    if ( verbose ){
        std::cout << "Write vtk " << fbody + ".vtk" << std::endl;
    }

    std::fstream file; if ( ! openOutFile( fbody.substr( 0, fbody.rfind( ".vtk" ) ) + ".vtk", & file ) ) { return; }
    std::map< std::string, RVector > data( dataMap );

    if ( cellCount() > 0 ){
        RVector tmp( cellCount() );
        std::transform( cellVector_.begin(), cellVector_.end(), &tmp[0], std::mem_fun( &Cell::marker ) );
        if ( !data.count( "_Marker" ) ) data.insert( std::make_pair( "_Marker",  tmp ) );
        if ( !data.count( "_Attribute" ) ) data.insert( std::make_pair( "_Attribute",  cellAttributes() ) );
    }

    bool binary = false;
    file.precision( 14 );
    file << "# vtk DataFile Version 3.0" << std::endl;
    if (commentString_.size() > 0) {
        file << commentString_ << std::endl;
    } else {
        file << "created by " << WHERE_AM_I << std::endl;
    }
    if ( binary ){
        file << "BINARY" << std::endl;
    } else {
        file << "ASCII" << std::endl;
    }
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    //** write nodes
    file << "POINTS " << nodeCount() << " double" << std::endl;

    for ( uint i = 0; i < nodeCount(); i ++ ){
        if ( binary ){
            file.write( (char*)&node( i ).pos()[0], sizeof( double ) );
            file.write( (char*)&node( i ).pos()[1], sizeof( double ) );
            file.write( (char*)&node( i ).pos()[2], sizeof( double ) );
        } else {
            file << node( i ).pos().x() << "\t"
                 << node( i ).pos().y() << "\t"
                 << node( i ).pos().z() << std::endl;
        }
    }
    if ( binary ){ file << std::endl; }

    //** write cells
    bool cells = true;
    if ( cells && cellCount() > 0 ){
        int idxCount = 0;
        for ( uint i = 0; i < cellCount(); i ++ ) idxCount += cell( i ).nodeCount()+1;

        file << "CELLS " << cellCount() << " " << idxCount << std::endl;

        long iDummy;
        for ( uint i = 0, imax = cellCount(); i < imax; i ++){
            if ( binary ){
                iDummy = cell( i ).nodeCount();
                file.write( (char*)&iDummy, sizeof( iDummy ) );
            } else {
                file << cell( i ).nodeCount() << "\t";
            }
            if ( cell( i ).rtti() == MESH_TETRAHEDRON10_RTTI && oldTet10NumberingStyle_ ){
                file << cell( i ).node( 0 ).id() << "\t" << cell( i ).node( 1 ).id() << "\t"
                    << cell( i ).node( 2 ).id() << "\t" << cell( i ).node( 3 ).id() << "\t"
                    << cell( i ).node( 4 ).id() << "\t" << cell( i ).node( 7 ).id() << "\t"
                    << cell( i ).node( 5 ).id() << "\t" << cell( i ).node( 6 ).id() << "\t"
                    << cell( i ).node( 9 ).id() << "\t" << cell( i ).node( 8 ).id();
            } else {

                for ( uint j = 0, jmax = cell( i ).nodeCount(); j < jmax; j ++){
                    if ( binary ){
                        iDummy = cell( i ).node(j).id();
                        file.write( (char*)&iDummy, sizeof( iDummy ) );
                    } else {
                        file << cell( i ).node(j).id() << "\t";
                    }
                }
            }
            if ( !binary ) file << std::endl;
        }
        if ( binary ) file << std::endl;

        file << "CELL_TYPES " << cellCount() << std::endl;
        iDummy = 10;
        for ( uint i = 0, imax = cellCount(); i < imax; i ++) {
            if ( binary ){
                file.write( (char*)&iDummy, sizeof( iDummy ) );
            } else {
                switch ( cell(i).rtti() ){
                    case MESH_EDGE_CELL_RTTI:
                    case MESH_EDGE_RTTI: file          << "3 "; break;
                    case MESH_EDGE3_RTTI:
                    case MESH_EDGE3_CELL_RTTI: file    << "21 "; break;
                    case MESH_TRIANGLE_RTTI: file      << "5 "; break;
                    case MESH_TRIANGLE6_RTTI: file     << "22 "; break;
                    case MESH_QUADRANGLE_RTTI: file    << "9 "; break;
                    case MESH_QUADRANGLE8_RTTI: file   << "23 "; break;
                    case MESH_TETRAHEDRON_RTTI: file   << "10 "; break;
                    case MESH_TETRAHEDRON10_RTTI: file << "24 "; break;
                    case MESH_TRIPRISM_RTTI: file      << "13 "; break;
                    case MESH_TRIPRISM15_RTTI: file    << "13 "; break;
                    case MESH_PYRAMID_RTTI: file       << "14 "; break;
                    case MESH_PYRAMID13_RTTI: file     << "14 "; break;
                    case MESH_HEXAHEDRON_RTTI: file    << "12 "; break;
                    case MESH_HEXAHEDRON20_RTTI: file  << "25 "; break;
                    default: std::cerr << WHERE_AM_I << " nothing known about." << cell(i).rtti()
                            << std::endl;
                }
            }
        }
        file << std::endl;
        
        //** write cell data
        file << "CELL_DATA " << cellCount() << std::endl;
        for ( std::map < std::string, RVector >::iterator it = data.begin(); it != data.end(); it ++ ){
            if ( verbose ){
                std::cout << it->first << " " << it->second.size() << std::endl;
            }
            if ( it->second.size() == (uint)cellCount() ){
                file << "SCALARS " << strReplaceBlankWithUnderscore( it->first )
                        << " double 1" << std::endl;
                file << "LOOKUP_TABLE default" << std::endl;

                for ( uint i = 0, imax = it->second.size(); i < imax; i ++) {
                    if ( binary ){
                        //file.write( (char*)&scaledValues[ i ], sizeof( double ) );
                    } else {
                        file << it->second[ i ] << " ";
                    }
                }
                file << std::endl;
            }
        }
    } else {  //   if !( cells && cellCount() > 0 ){
        //** write boundaries
        if ( boundaryCount() > 0 ){
            int idxCount = 0;
            for ( uint i = 0; i < boundaryCount(); i ++ ) idxCount += boundary( i ).nodeCount()+1;
            file << "CELLS " << boundaryCount() << " " << idxCount << std::endl;

            long iDummy;
            for ( int i = 0, imax = boundaryCount(); i < imax; i ++ ){
                if ( binary ){
            //        iDummy = cell( i ).nodeCount();
            // 	      file.write( (char*)&iDummy, sizeof( iDummy ) );
                } else {
                    file << boundary( i ).nodeCount() << "\t";
                }

                for ( uint j = 0, jmax = boundary( i ).nodeCount(); j < jmax; j ++){
                    if ( binary ){
// 	          iDummy = cell( i ).node(j).id();
// 	          file.write( (char*)&iDummy, sizeof( iDummy ) );
                    } else {
                        file << boundary( i ).node( j ).id() << "\t";
                    }
                }
                if ( !binary ) file << std::endl;
            }
            if ( binary ) file << std::endl;

            file << "CELL_TYPES " << boundaryCount() << std::endl;
            iDummy = 10;
            for ( uint i = 0, imax = boundaryCount(); i < imax; i ++) {
                if ( binary ){
                    file.write( (char*)&iDummy, sizeof( iDummy ) );
                } else {
                    switch ( boundary( i ).rtti() ){
                        case MESH_BOUNDARY_NODE_RTTI: file     << "1 "; break;
                        case MESH_EDGE_RTTI: file              << "3 "; break;
                        case MESH_EDGE3_RTTI: file             << "21 "; break;
                        case MESH_TRIANGLEFACE_RTTI: file      << "5 "; break;
                        case MESH_TRIANGLEFACE6_RTTI: file     << "22 "; break;
                        case MESH_QUADRANGLEFACE_RTTI: file    << "9 "; break;
                        case MESH_QUADRANGLEFACE8_RTTI: file   << "23 "; break;
                    default: std::cerr << WHERE_AM_I << " nothing know about." << boundary( i ).rtti()
                            << std::endl;
                    }
                }
            }
            file << std::endl;
        
        
            RVector tmp( boundaryCount() );
            std::transform( boundaryVector_.begin(), boundaryVector_.end(),
                            &tmp[0], std::mem_fun( &Boundary::marker ) );
        
            if ( !data.count( "_Marker" ) ) data.insert( std::make_pair( "_Marker",  tmp ) );
        
            //** write boundary data
            file << "CELL_DATA " << boundaryCount() << std::endl;
            for ( std::map < std::string, RVector >::iterator it = data.begin(); it != data.end(); it ++ ){
                if ( verbose ){
                    std::cout << it->first << " " << it->second.size() << std::endl;
                }
                if ( it->second.size() == (uint)boundaryCount() ){
                    file << "SCALARS " << strReplaceBlankWithUnderscore( it->first )
                            << " double 1" << std::endl;
                    file << "LOOKUP_TABLE default" << std::endl;

                    for ( uint i = 0, imax = it->second.size(); i < imax; i ++) {
                        if ( binary ){
                            //file.write( (char*)&scaledValues[ i ], sizeof( double ) );
                        } else {
                            file << it->second[ i ] << " ";
                        }
                    }
                    file << std::endl;
                }
            }
        } // if ( boundaryCount() > 0 )
        
    } // else write boundaries
    
    //** write point data
    file << "POINT_DATA " << nodeCount() << std::endl;
    for ( std::map < std::string, RVector >::iterator it = data.begin(); it != data.end(); it ++ ){
        if ( it->second.size() == (uint)nodeCount() ){
            file << "SCALARS " << strReplaceBlankWithUnderscore( it->first )
                 << " double 1" << std::endl;
            file << "LOOKUP_TABLE default" << std::endl;

            for ( uint i = 0, imax = nodeCount(); i < imax; i ++) {
                file << it->second[ i ] << " ";
            }
            file << std::endl;
        }
    }
        
    //** write point vector data
    if ( vec.size() == nodeCount() ){
        if ( verbose ){
            std::cout << "write vector field data" << std::endl;
        }
        file << "VECTORS vec double" << std::endl;
            
        for ( uint i = 0; i < vec.size(); i ++ ){
            file << vec[i][0] << " " << vec[i][1] << " "<< vec[i][2] << " " << std::endl;
        }
    } else {
        if ( vec.size() > 0 ){
            std::cerr << "Vector data size does not match node size: " << vec.size() << " " << nodeCount() << std::endl;
        }
    }
    
    file.close();
}

void Mesh::importVTK( const std::string & fbody ) {
    this->clear();
    std::fstream file; openInFile( fbody.substr( 0, fbody.rfind( ".vtk") ) + ".vtk", &file );

    std::vector < std::string > row ;
    getline( file, commentString_ ); //** vtk version line
    getline( file, commentString_ ); //** comment line
    while ( !file.eof() ){
        row = getRowSubstrings( file );
      //  std::cout << row.size() << std::endl;
        if ( row.size() ){
            if ( row[ 0 ] == "ASCII" ){
                break;
            } else if ( row[ 0 ] == "BINARY" ){
                THROW_TO_IMPL
                break;
            }
        }
    }
    row = getRowSubstrings( file );
    //std::cout << row.size() << std::endl;
    if ( row.back() == "UNSTRUCTURED_GRID" ){
        
        //** End reading header

        while ( !file.eof() ){
            row = getRowSubstrings( file );
            if ( row.size() ){
                if ( row[ 0 ] == "POINTS" ) readVTKPoints_( file, row );
                else if ( row[ 0 ] == "CELLS" ) readVTKCells_( file, row );
                else if ( row[ 0 ] == "SCALARS" ) readVTKScalars_( file, row );
            }
        } 
        
    } else if ( row.back() == "POLYDATA" ){
        while ( !file.eof() ){
            row = getRowSubstrings( file );
            if ( row.size() ){
                if ( row[ 0 ] == "POINTS" ) readVTKPoints_( file, row );
                else if ( row[ 0 ] == "POLYGONS" ) readVTKPolygons_( file, row );
            }
        } 
    } else {
            THROW_TO_IMPL
    }
    
    this->showInfos();
    file.close();
}

void Mesh::readVTKPoints_( std::fstream & file, const std::vector < std::string > & row ){
    uint nVerts = toInt( row[ 1 ] );
    //std::cout << "nVerts: " << nVerts << std::endl;
    double x = 0.0, y = 0.0, z = 0.0;
    for ( uint i = 0; i < nVerts; i ++ ) {
        file >> x >> y >> z;
        this->createNode( x, y, z );
    }
}

void Mesh::readVTKCells_( std::fstream & file, const std::vector < std::string > & row ){
    uint nCells = toInt( row[ 1 ] );
    //std::cout << "nCells: " << nCells << std::endl;
    uint nNodes = 0;
    uint id;
    std::vector < Node * > nodes;
    for ( uint i = 0; i < nCells; i ++ ) {
        file >> nNodes;
        nodes.resize( nNodes );
        for ( uint j = 0; j < nNodes; j ++ ) {
            file >> id;
            nodes[ j ] = &this->node( id );
        }
        this->createCell( nodes );
    }
}

void Mesh::readVTKPolygons_( std::fstream & file, const std::vector < std::string > & row ){
    uint nPoly = toInt( row[ 1 ] );
    uint nNodes = 0;
    uint id;
    std::vector < Node * > nodes;
    for ( uint i = 0; i < nPoly; i ++ ) {
        file >> nNodes;
        nodes.resize( nNodes );
        for ( uint j = 0; j < nNodes; j ++ ) {
            file >> id;
            nodes[ j ] = &this->node( id );
        }
        this->createBoundary( nodes );
    }
}

// template <class Arg, class Result> class bind : std::unary_function< Arg, Result >{
// public:
//     Result operator() ( Arg & o ) { return toDouble( o ); }
// };
//   struct unary_function{
//     typedef Arg argument_type;
//     typedef Result result_type;
//   };
// template < class ret > class bind{
// public:
//     bind( ) {
//         funct_ = &toDouble;
//     }
//     template < class type > ret operator() ( type & o ) { return toDouble( o ); }
// protected:
//     //void (Region::*)
//     size_t &funct_;
//     //(toDouble*) funct_;
// };

void Mesh::readVTKScalars_( std::fstream & file, const std::vector < std::string > & row ){
    std::string name( row[ 1 ] );
    std::vector < std::string > r( getRowSubstrings( file ) );
    if ( r.size() ){
        if ( r[ 0 ] == "LOOKUP_TABLE" ){
            r = getRowSubstrings( file );
        }
    }
    RVector data( r.size() );
//     std::cout.precision( 14 );
//     std::cout << *r.begin() << std::endl;
//     std::cout << toDouble( *r.begin() ) << std::endl;
//     std::cout << std::bind1st( &toDouble )( *r.begin() ) << std::endl;
    //std::cout << bind< double, toDouble >( &toDouble )( *r.begin() ) << std::endl;

    //std::copy( r.begin(), r.end(), data.begin(), bind< double >( toDouble );
    for ( uint i = 0; i < data.size(); i ++ ) data[ i ] = toDouble( r[ i ] );
    addExportData( name, data );
}


void Mesh::importVTU( const std::string & fbody ) {
    this->clear();
    THROW_TO_IMPL
}

void Mesh::exportVTU( const std::string & fbody, bool binary ) const {
    std::string filename( fbody );
    if ( filename.rfind( ".vtu" ) == std::string::npos ){
        filename = fbody.substr( 0, filename.rfind( ".vtk" ) ) + ".vtu";
    }
    std::fstream file; if ( ! openOutFile( filename, & file ) ) { return ; }
    file.precision( 14 );
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    file << "<UnstructuredGrid>" << std::endl;

    std::map< std::string, RVector > data( exportDataMap_ );
    if ( cellCount() > 0 ){
        RVector tmp( cellCount() );
        std::transform( cellVector_.begin(), cellVector_.end(), &tmp[0], std::mem_fun( &Cell::marker ) );
        if ( !data.count( "_Marker" ) ) data.insert( std::make_pair( "_Marker",  tmp ) );
        if ( !data.count( "_Attribute" ) ) data.insert( std::make_pair( "_Attribute",  cellAttributes() ) );
    }
    addVTUPiece_( file, *this, data );

    file << "</UnstructuredGrid>" << std::endl;
    file << "</VTKFile>" << std::endl;
    file.close();
}

void Mesh::exportBoundaryVTU( const std::string & fbody, bool binary ) const {
    std::string filename( fbody );
    if ( filename.rfind( ".vtu" ) == std::string::npos ){
        filename = fbody.substr( 0, filename.rfind( ".vtk" ) ) + ".vtu";
    }
    std::fstream file; if ( ! openOutFile( filename, & file ) ) { return ; }


    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    file << "<UnstructuredGrid>" << std::endl;

    std::vector < Boundary * > bs;
    for ( uint i = 0; i < boundaryCount(); i ++ ) {
        if ( boundary( i ).marker() != 0.0 ) bs.push_back( &boundary( i ) );
//         std::cout << boundary(i)<< " " << boundary(i).center() << " " << boundary(i).marker() << std::endl;
    }

    Mesh boundMesh;
    boundMesh.createMeshByBoundaries( *this, bs );
//    boundMesh.showInfos();
    //for ( uint i =0; i < boundMesh.nodeCount(); i ++ ) std::cout << boundMesh.node( i )<< std::endl;
    std::map< std::string, RVector > boundData;

    RVector tmp( boundMesh.boundaryCount() );
    std::transform( boundMesh.boundaries().begin(), boundMesh.boundaries().end(), &tmp[ 0 ],
                    std::mem_fun( &Boundary::marker ) );
    if ( !boundData.count( "_BoundaryMarker" ) ) boundData.insert( std::make_pair( "_BoundaryMarker",  tmp ) );
    boundMesh.exportVTK( fbody, boundData );
    addVTUPiece_( file, boundMesh, boundData );

    file << "</UnstructuredGrid>" << std::endl;
    file << "</VTKFile>" << std::endl;
    file.close();
}

void Mesh::addVTUPiece_( std::fstream & file, const Mesh & mesh,
                  const std::map < std::string, RVector > & data ) const{

    bool binary = false;

    std::vector < MeshEntity * > cells;

    if ( mesh.cellCount() == 0 && mesh.boundaryCount () > 0 ){
        cells.reserve( mesh.boundaryCount() );
        for ( uint i = 0; i < mesh.boundaryCount(); i ++ ) cells.push_back( & mesh.boundary( i ) );
    } else {
        cells.reserve( mesh.cellCount() );
        for ( uint i = 0; i < mesh.cellCount(); i ++ ) cells.push_back( & mesh.cell( i ) );
    }

    uint nNodes = mesh.nodeCount();
    uint nCells = cells.size();

    file << "<Piece NumberOfPoints=\"" << nNodes << "\" NumberOfCells=\"" << nCells << "\">" << std::endl;

    file << "<Points>" << std::endl;
    file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    for ( uint i = 0; i < mesh.nodeCount(); i ++ ) {
        file << mesh.node( i ).x() << " " << mesh.node( i ).y() << " " << mesh.node( i ).z() << " ";
    }
    file << std::endl << "</DataArray>" << std::endl;
    file << "</Points>" << std::endl;

    file << "<Cells>" << std::endl;
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    for ( uint i = 0; i < nCells; i ++ ) {
        MeshEntity * cell = cells[ i ];
        if ( cell->rtti() == MESH_TETRAHEDRON10_RTTI ){

            file << cell->node( 0 ).id() << " " << cell->node( 1 ).id() << " "
                << cell->node( 2 ).id() << " " << cell->node( 3 ).id() << " "
                << cell->node( 4 ).id() << " " << cell->node( 7 ).id() << " "
                << cell->node( 5 ).id() << " " << cell->node( 6 ).id() << " "
                << cell->node( 9 ).id() << " " << cell->node( 8 ).id() << " ";
        } else for ( uint j = 0; j < cell->nodeCount(); j ++ ) file << cell->node( j ).id() << " ";
    }
    file << std::endl << "</DataArray>" << std::endl;
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    uint count = 0;
    for ( uint i = 0; i < nCells; i ++ ) {
        count += cells[ i ]->nodeCount();
        file << count << " ";
    }
    file << std::endl << "</DataArray>" << std::endl;
    file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
    for ( uint i = 0; i < nCells; i ++ ) {
        switch ( cells[ i ]->rtti() ){
            case MESH_BOUNDARY_NODE_RTTI: file     << "1 "; break;
            case MESH_EDGE_CELL_RTTI:
            case MESH_EDGE_RTTI: file          << "3 "; break;
            case MESH_EDGE3_CELL_RTTI:
            case MESH_EDGE3_RTTI: file         << "21 "; break;
            case MESH_TRIANGLEFACE_RTTI:
            case MESH_TRIANGLE_RTTI: file      << "5 "; break;
            case MESH_TRIANGLEFACE6_RTTI:
            case MESH_TRIANGLE6_RTTI: file     << "22 "; break;
            case MESH_QUADRANGLEFACE_RTTI:
            case MESH_QUADRANGLE_RTTI: file    << "9 "; break;
            case MESH_QUADRANGLEFACE8_RTTI:
            case MESH_QUADRANGLE8_RTTI: file   << "23 "; break;
            case MESH_TETRAHEDRON_RTTI: file   << "10 "; break;
            case MESH_TETRAHEDRON10_RTTI: file << "24 "; break;
            case MESH_HEXAHEDRON_RTTI: file    << "12 "; break;
            case MESH_HEXAHEDRON20_RTTI: file  << "25 "; break;
            default: std::cerr << WHERE_AM_I << " nothing know about." << cells[ i ]->rtti() << std::endl;
        }
    }
    file << std::endl << "</DataArray>" << std::endl;
    file << "</Cells>" << std::endl;

    uint nodeData = 0, cellData = 0;
    for ( std::map < std::string, RVector >::const_iterator it = data.begin(); it != data.end(); it ++ ){
        if ( it->second.size() == nNodes ) nodeData++;
        else if ( it->second.size() == nCells ) cellData++;
        else {
            std::cerr << WHERE_AM_I << " dont know how to handle data array: " << it->first
                        << " with size " << it->second.size() << " nodesize = " << nNodes
                        << " cellsize = " << nCells << std::endl;
        }
    }

    if ( nodeData > 0 ){
        file << "<PointData>" << std::endl;
        for ( std::map < std::string, RVector >::const_iterator it = data.begin();
              it != data.end(); it ++ ){
            if ( it->second.size() == nNodes ) {
                 file << "<DataArray type=\"Float32\" Name=\""
                      << it->first << "\" format=\"ascii\">" << std::endl;
                 for ( uint i = 0; i < it->second.size(); i ++ ) file << it->second[ i ] << " ";
                 file << std::endl << "</DataArray>" << std::endl;
            }
        }
        file << "</PointData>" << std::endl;
    }

     if ( cellData > 0 ){
        file << "<CellData>" << std::endl;
        for ( std::map < std::string, RVector >::const_iterator it = data.begin();
              it != data.end(); it ++ ){
            if ( it->second.size() == nCells ) {
                 file << "<DataArray type=\"Float64\" Name=\"" << it->first;
                 if ( !binary ){
                    file << "\" format=\"ascii\">" << std::endl;
                    for ( uint i = 0; i < it->second.size(); i ++ ) file << it->second[ i ] << " ";
                 } else {
                    file << "\" format=\"binary\">";
                    file.write( (char*)&it->second[ 0 ], it->second.size() * sizeof( double ) );
                 }


                 file << std::endl << "</DataArray>" << std::endl;
            }
        }
        file << "</CellData>" << std::endl;
    }

    file << "</Piece>" << std::endl;
}

void Mesh::importMod( const std::string & filename ){
    RMatrix mat;
    std::vector < std::string > comments;
    loadMatrixCol( mat, filename, comments );
//     std::cout << comments.size() << std::endl;
//     std::cout << mat.size() << std::endl;
//     std::cout << mat[ 0 ].size() << std::endl;

    RVector x( unique( sort( cat( mat[ 0 ], mat[ 1 ] ) ) ) );
    RVector y( unique( sort( cat( mat[ 2 ], mat[ 3 ] ) ) ) * -1.0 );
    create2DGrid( x, y );
    if ( comments.size() > 4 && mat.cols() > 4 ) addExportData( comments[ 4 ], mat[ 4 ] );
    if ( comments.size() > 5 && mat.cols() > 5 ) addExportData( comments[ 5 ], mat[ 5 ] );
}

void Mesh::importSTL( const std::string & fileName, bool isBinary  ){
    double tolerance = 1e-3;
    
    std::vector < RVector3 > allVerts;
    
    if ( !isBinary ){ // try import ascii
        std::fstream file; openInFile( fileName, & file );

        std::vector < std::string > row;

        row = getNonEmptyRow( file );
        if ( row[ 0 ] != "solid" ) {
            file.close();
            importSTL( fileName, true );
        }
        
        bool finish = false;
        while ( !finish ){
            row = getNonEmptyRow( file );
            if ( row[ 0 ] == "facet" && row[ 1 ] == "normal" ){ //** facet normal  0.0  0.0  0.0
                row = getNonEmptyRow( file );  //** outer loop
                row = getNonEmptyRow( file ); //** vertex x y z;
                allVerts.push_back( RVector3( toDouble( row[ 1 ] ), toDouble( row[ 2 ] ), toDouble( row[ 3 ] ) ) );
                row = getNonEmptyRow( file ); //** vertex x y z
                allVerts.push_back( RVector3( toDouble( row[ 1 ] ), toDouble( row[ 2 ] ), toDouble( row[ 3 ] ) ) );
                row = getNonEmptyRow( file ); //** vertex x y z
                allVerts.push_back( RVector3( toDouble( row[ 1 ] ), toDouble( row[ 2 ] ), toDouble( row[ 3 ] ) ) );
                row = getNonEmptyRow( file );  //** endloop
            row = getNonEmptyRow( file );  //** endfacet;
            } else finish = true;
        }
        file.close();
    } else { // import Binary Format
        UNTESTED
        
        FILE * file; file = fopen( fileName.c_str(), "r+b" );

        char header[ 80 ];
        Index ret = 0; 
        ret = fread( &header, 1, 80, file );
        if ( ret == 0 ) throwError( 1, WHERE_AM_I + " Oops" );
        
        int nFaces = 0;
        ret = fread( &nFaces, 4, 1, file );
        if ( ret == 0 ) throwError( 1, WHERE_AM_I + " Oops" );

        std::vector < RVector3 > allVerts;

        float rd[ 48 ];
        char padding[ 2 ];
        for ( int i = 0; i < nFaces; i ++ ){
            ret = fread( &rd, 4, 12, file );
            if ( ret == 0 ) throwError( 1, WHERE_AM_I + " Oops" );
            
            allVerts.push_back( RVector3( rd[ 3 ], rd[ 4 ],  rd[ 5 ] ).round( tolerance * 0.01 ) );
            allVerts.push_back( RVector3( rd[ 6 ], rd[ 7 ],  rd[ 8 ] ).round( tolerance * 0.01 ) );
            allVerts.push_back( RVector3( rd[ 9 ], rd[ 10 ], rd[ 11 ] ).round( tolerance * 0.01 ) );

            ret = fread( &padding, 1, 2, file );
            if ( ret == 0 ) throwError( 1, WHERE_AM_I + " Oops" );
        }
      
        fclose( file );
    } // end import binary STL format
    
    Node *n1, *n2, *n3;
    if ( allVerts.size() % 3 == 0 && allVerts.size() > 0 ){
        for ( uint i = 0; i < allVerts.size() / 3; i ++ ){
            n1 = createNodeWithCheck( allVerts[ i * 3 ], tolerance );
            n2 = createNodeWithCheck( allVerts[ i * 3 + 1 ], tolerance );
            n3 = createNodeWithCheck( allVerts[ i * 3 + 2 ], tolerance );
            this->createTriangleFace( *n1, *n2, *n3, 0);
        }
    } else {
        throwError(1,  WHERE_AM_I + " there is something wrong in ascii-stl-format "
                + toStr( allVerts.size() ) + " " + toStr( allVerts.size() % 3 ) );
    }
}

int Mesh::exportMidCellValue( const std::string & fileName,
                              const RVector & data1, const RVector & data2 ) const {

    RMatrix mat( dimension_, cellCount() );
    for ( uint j = 0; j < cellCount(); j ++ ) {
        for ( uint i = 0; i < dimension_; i ++ ) {
            mat[ i ][ j ] = cell( j ).center()[ i ];
        }
    }

    if ( data1.size() == cellCount() ) mat.push_back( data1 );
    if ( data2.size() == cellCount() ) mat.push_back( data2 );

    return saveMatrixCol( mat, fileName );
}

void Mesh::exportAsTetgenPolyFile( const std::string & filename ){
    std::fstream file; openOutFile( filename.substr( 0, filename.rfind( ".poly" ) ) + ".poly", & file );

    uint nverts = nodeCount();
    uint nfacets = boundaryCount();

    //  file << "# nverties dimension nattrib boolbndrymarker" << endl;
    file << nverts << "\t3\t0\t1" << std::endl;

    //file << "#    pointNr. x y z marker" << endl;
    file.setf( std::ios::scientific, std::ios::floatfield );
    file.precision( 12 );
    for ( uint i = 0; i < nverts; i++ ){
        file << i << "\t" << node( i ).x( )
                << "\t" << node( i ).y( )
                << "\t" << node( i ).z( )
                << "\t" << node( i ).marker() << std::endl;
    }

    //file << "# nfacets boolbndrymarker" << endl;
    size_t nPolyFacet = 0;
    size_t nNodesPoly = 0;
    file << nfacets << "\t1" << std::endl;
    for ( size_t i = 0, imax = nfacets;i < imax; i ++ ){

        nPolyFacet = 1;
        file << nPolyFacet << "\t0\t" << boundary( i ).marker() << std::endl;
        //file << nPolyFacet << "\t" << facet( i ).holeCount() << "\t" << facet( i ).boundaryMarker() << endl;

        nNodesPoly = boundary( i ).nodeCount();
        file << nNodesPoly << " \t";

        for ( size_t k = 0, kmax = nNodesPoly; k < kmax; k++ ){
            file << boundary( i ).node( k ).id() << "\t";
        }
        file << std::endl;
    }

    //file << "# nholes" << endl;
    file << 0 << std::endl;
    //file << "#    hole x y z" << endl;
//     for ( int i = 0, imax = vecpHoles_.size(); i < imax ; i ++ ){
//         file << i << "\t"
//         << vecpHoles_[ i ]->x() << "\t"
//         << vecpHoles_[ i ]->y() << "\t"
//         << vecpHoles_[ i ]->z() << endl;
//     }

  //file << "# nregions" << endl;
  file << 0 << std::endl;
  //file << "#    region x y z attribute maxarea" << endl;
//   for ( int i = 0, imax = regions_.size(); i < imax ; i ++ ){
//     file << i << "\t"
//      << regions_[ i ]->x() << "\t"
//      << regions_[ i ]->y() << "\t"
//          << regions_[ i ]->z() << "\t"
//      << regions_[ i ]->attribute() << "\t"
//      << regions_[ i ]->dx() << endl;
//   }
  file.close();

}

} // namespace GIMLI{
