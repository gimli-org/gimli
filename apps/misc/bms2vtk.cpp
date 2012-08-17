/***************************************************************************
 *   Copyright (C) 2006-2007 by the resistivity.net development team       *
 *   Carsten Rücker carsten@resistivity.net                                *
 *   Thomas Günther thomas@resistivity.net                                 *
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

#include <gimli.h>
#include <vectortemplates.h>
#include <meshgenerators.h>
#include <mesh.h>
#include <optionmap.h>

#include <string>
//using namespace std;
using namespace GIMLI;

#define vcout if ( verbose ) std::cout
#define dcout if ( debug ) std::cout
#define DEBUG if ( debug )

//** should be moved to some other place into the library
void setAllNeumannBoundaryConditions( Mesh & mesh ){
    RBoundingBox bbox( mesh.boundingBox() );
    mesh.createNeighbourInfos();
    for ( uint i = 0; i < mesh.boundaryCount(); i ++ ){
        RVector3 cen( mesh.boundary( i ).center() );

        for ( uint dim = 0; dim < mesh.dimension(); dim++){
            if ( cen[ dim ] == bbox.min()[ dim ] || cen[ dim ] == bbox.max()[ dim ] ){
                 mesh.boundary( i ).setMarker( MARKER_BOUND_HOMOGEN_NEUMANN );
            }
        }

    }
}

void setDefaultWorldBoundaryConditions( Mesh & mesh ){
    RBoundingBox bbox( mesh.boundingBox() );
    mesh.createNeighbourInfos();
    for ( uint i = 0; i < mesh.boundaryCount(); i ++ ){
        RVector3 cen( mesh.boundary( i ).center() );

        for ( uint dim = 0; dim < mesh.dimension(); dim++){
            if ( cen[ dim ] == bbox.min()[ dim ] || cen[ dim ] == bbox.max()[ dim ] ){
                 mesh.boundary( i ).setMarker( MARKER_BOUND_MIXED );
            }
        }

        if ( cen[ mesh.dimension() -1 ] == bbox.max()[ mesh.dimension()-1 ] ){
            mesh.boundary( i ).setMarker( MARKER_BOUND_HOMOGEN_NEUMANN );
        }
    }
}

int main( int argc, char *argv [] ){
    bool readVecs = false, allNeumannBoundary = false, saveBin = false, exportVTK = false;
    bool refineP = false, refineH = false, exportBoundary = false, is2dMesh = false, makeLayers = false;
    int equiBoundary = 0, prolongBoundary = 0, prolongFactor = 5;
    std::string meshFile, outFile = NOT_DEFINED, modelFileName = NOT_DEFINED;

    OptionMap oMap;
    oMap.setDescription("Description. BMS2VTK - Convert BMS to VTK or create new ones from vectors\n");
    oMap.addLastArg( meshFile, "Mesh file or vector basename" );
    oMap.add( is2dMesh,          "2"  , "is2dMesh"          , "mesh is 2d instead of 3d (with -R)" );
    oMap.add( readVecs,          "R"  , "readVectors"       , "Read vectors .x/y/z from files" );
    oMap.add( makeLayers,        "L"  , "makeLayers"        , "Give each layer another boundary" );
    oMap.add( allNeumannBoundary,"N"  , "allNeumannBoundary", "All boundaries are Neumann" );
    oMap.add( refineH,           "H"  , "refineH"           , "Refine mesh spatially" );
    oMap.add( exportBoundary,    "B"  , "exportBoundary"    , "Export boundary" );
    oMap.add( exportVTK,         "V"  , "exportVTK"         , "Export vtk file" );
    oMap.add( refineP,           "P"  , "refineP"           , "Refine mesh polynomially" );
    oMap.add( equiBoundary,      "e:" , "equiBoundary"      , "equidistant boundary around hex mesh" );
    oMap.add( prolongBoundary,   "p:" , "prolongBoundary"   , "prolongated boundary around equiBoundary" );
    oMap.add( prolongFactor,     "f:" , "prolongFactor"     , "prolongation factor for prolongBoundary" );
    oMap.add( outFile,           "o:" , "outFile"           , "filename for output" );
    oMap.add( modelFileName,     "a:" , "modelFile"         , "model file to include" );
    oMap.parse( argc, argv );

    if ( outFile == NOT_DEFINED ) outFile = meshFile;

    Mesh mesh;
    if ( readVecs ) { //! create hexahedral mesh from given x/y/z vectors
        RVector x, y, z;
        load( x, meshFile + ".x" );
        load( y, meshFile + ".y" );
        load( z, meshFile + ".z" );
        z = sort( z );
        double xmin = x[ 0 ], xmax = x[ x.size() -1 ];
        double ymin = y[ 0 ], ymax = y[ y.size() -1 ];
        double zmin = z[ 0 ];
        if ( equiBoundary > 0 || prolongBoundary > 0 ) {
          //! create prolongation vector
          RVector addVec( equiBoundary + prolongBoundary , 1.0 );
          for ( size_t i = 1 ; i < addVec.size() ; i++ ) addVec[ i ] += addVec[ i - 1 ];
          if ( prolongBoundary > 0 ) {
              double dx = 1;
              for ( int i = 0 ; i < prolongBoundary ; i ++ ) {
                  dx *= prolongFactor;
                  addVec[ equiBoundary + i ] = addVec[ equiBoundary + i -1 ] + dx;
              }
          }
          std::cout << addVec << std::endl;
          RVector revVec = addVec;
          for ( size_t i = 0 ; i < revVec.size() ; i++ ) revVec[ i ] = addVec[ revVec.size() - i - 1 ];
          std::cout << revVec << std::endl;
          x = cat( cat( RVector( revVec * ( x[ 0 ] - x[ 1 ]) + x[ 0 ] ), x ),
              RVector( addVec * ( x[ x.size() - 1 ] - x[ x.size() - 2 ] ) + x[ x.size()-1 ] ) );
          std::cout << x << std::endl;
          y = cat( cat( RVector( revVec * ( y[ 0 ] - y[ 1 ]) + y[ 0 ] ), y ),
              RVector( addVec * ( y[ y.size() - 1 ] - y[ y.size() - 2 ] ) + y[ y.size()-1 ] ) );
          std::cout << y << std::endl;
          z = cat( RVector( revVec * ( z[ 0 ] - z[ 1 ] ) + z[ 0 ] ), z );
          std::cout << z << std::endl;
        }
        if ( is2dMesh ) {
            mesh = createMesh2D( x, z );
        } else {
            mesh = createMesh3D( x, y, z );
        }
        if ( equiBoundary + prolongBoundary > 0 ) {
            //! determine cell marker
            for ( size_t i = 0 ; i < mesh.cellCount() ; i++ ) {
                int num = 1;
                RVector3 midpoint = mesh.cell( i ).center();
                if ( is2dMesh ) {
                    if ( midpoint[ 0 ] > xmin && midpoint[ 0 ] < xmax && midpoint[ 1 ] > zmin ) {
                        if ( makeLayers ) { //! each layer another number
                            for ( int k = z.size()-1 ; k>=0 ; k-- ) {
                                if ( midpoint[ 1 ] < z[ k ] ) num++;
                            }
                        } else { //! normal 1/2 mesh
                            num = 2;
                        }
                    }
                } else { //! a 3D mesh                    
                    if ( midpoint[ 0 ] > xmin && midpoint[ 0 ] < xmax && midpoint[ 1 ] > ymin && midpoint[ 1 ] < ymax && midpoint[ 2 ] > zmin ) {
                        num = 2;
                    }
                }
                mesh.cell( i ).setMarker( num );

            }
        }
        if ( allNeumannBoundary ) {
            setAllNeumannBoundaryConditions( mesh );
        } else {
            setDefaultWorldBoundaryConditions( mesh );
        }
        saveBin = true;
    } else { // read bms file
        mesh.load( meshFile );
    }
    std::cout << "Mesh:\t"; mesh.showInfos();
    if ( refineH ) {
        mesh = mesh.createH2();
        std::cout << "MeshH:\t"; mesh.showInfos();
        saveBin = true;
    }
    if ( refineP ) {
        mesh = mesh.createP2( );
        std::cout << "MeshP:\t"; mesh.showInfos();
        saveBin = true;
    }
    if ( modelFileName != NOT_DEFINED ) {
        RVector model( modelFileName );
        mesh.addExportData( modelFileName, model );
    }
    if ( saveBin ) mesh.saveBinary( outFile );
    if ( exportVTK || !saveBin ) mesh.exportVTK( outFile );
    if ( exportBoundary ) mesh.exportBoundaryVTU( "boundary.vtu" );

    return EXIT_SUCCESS;
}
