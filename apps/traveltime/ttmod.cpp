/***************************************************************************
 *   Copyright (C) 2006-2008 by the resistivity.net development team       *
 *   Thomas Günther thomas@resistivity.net                                 *
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

#include <gimli.h>
#include <mesh.h>
#include <datacontainer.h>
#include <inversion.h>
#include <optionmap.h>
#include <modellingbase.h>

#include <string>
using namespace std;
using namespace GIMLI;

#define vcout if ( verbose ) cout
#define dcout if ( debug ) cout
#define DEBUG if ( debug )

#include <ttdijkstramodelling.h>

//** MAIN
int main( int argc, char *argv [] ) {

    bool isVelocity = true, isSlowness = false, doError = false, addNoise = false;
    int nSegments = 6, verboseCount = 0;
    double relativeInnerMaxEdgeLength = 0.01, errTime = 0.001, errPerc = 0;
    string meshFilename( NOT_DEFINED ), modelFilename( NOT_DEFINED ), outFilename( "out.dat" );
    string dataFileName(NOT_DEFINED ), attributemap ( NOT_DEFINED );

    OptionMap oMap;
    oMap.setDescription("Description. TTMod - Travel time modelling using Dijkstra algorithm.\n");
    oMap.addLastArg( dataFileName, "Data file" );
    oMap.add( verboseCount,       "v" , "verbose", "Verbose mode (2 times for debug mode)." );
    oMap.add( doError,            "E" , "doError", "Calculate error as well." );
    oMap.add( isVelocity,         "V" , "isVelocity", "Input is velocity" );
    oMap.add( isSlowness,         "S" , "isSlowness", "Input is slowness" );
    oMap.add( addNoise,           "N" , "addNoise", "Noisify data." );
    oMap.add( meshFilename,       "p:" , "meshFilename", "Mesh file." );
    oMap.add( outFilename,        "o:" , "outFilename", "Output file name [out.dat]." );
    oMap.add( errPerc,            "e:" , "errorPercent", "Percentage error (for -N/-E) [0]." );
    oMap.add( errTime,            "t:" , "errorTime", "Absolute traveltime error (for -N/-E) [1ms]." );
    oMap.add( attributemap,       "a:" , "attributeMap", "Map file associating velocity/slowness [all 1m/s]." );
    oMap.add( modelFilename,      "m:" , "modelFilename", "Model file instead of attribute map." );
    oMap.parse( argc, argv );

    bool verbose = ( verboseCount > 0 );//, debug = ( verboseCount > 1 );
    if ( isSlowness ) isVelocity = true;

    DataContainer data( dataFileName );
    if ( verbose ) data.showInfos();

    Mesh mesh;
    if ( meshFilename != NOT_DEFINED ) {
        mesh.load( meshFilename );
    } else {
        cerr << WHERE_AM_I << " no mesh given. Creating one. ..."  << endl;
        mesh.createClosedGeometryParaMesh( data.sensorPositions(),
                                           nSegments, relativeInnerMaxEdgeLength,
                                           data.additionalPoints() );
    }
    vcout << "Mesh: ";
    mesh.showInfos( );

    if ( attributemap != NOT_DEFINED ) {
        std::map < float, float > attMap( loadFloatMap( attributemap ) );
        if ( isVelocity ) {
            for (std::map < float, float >::iterator it = attMap.begin(); it!=attMap.end(); it ++ ) it->second = 1.0 / it->second;
        }
        mesh.mapCellAttributes( attMap );
    }

    /*! set up TT modeling class; */
    TravelTimeDijkstraModelling f( mesh, data );

//  vcout << "model size=" << model.size() << " min/max=" << min(model) << "/" << max(model) << endl;
    RVector traveltime( data.size() );
    if ( modelFilename != NOT_DEFINED ) {
        RVector model;
        load( model, modelFilename );
        if ( isVelocity ) model = 1.0 / model;
        SparseMapMatrix < double, size_t > W;
        f.createJacobian( W, model );
        traveltime = W * model;
        std::cout << "rms(data, modelResponse) = " << rms( data( "t" ), traveltime ) << "s." << std::endl;
        RVector coverage( model.size() );
        coverage = transMult( W, RVector( data.size(), 1.0 ) );
        save( coverage, "coverage.vec", Ascii );
//    } else {
//        traveltime = f.response( );
    }
    cout << "min/max traveltime = " << min( traveltime ) << "/" << max( traveltime ) << "s" << endl;
    data.set( "t", traveltime );

    if ( doError ) {
        /*! estimate error by error time and percentage value; */
        vcout << "Estimate error: " << errPerc << "% + " << errTime << "s" << endl;
        data.set( "err", errTime / data("t") + errPerc / 100.0 ); // always relative error
        vcout << "Data error:" << " min = " << min( data("err") ) * 100 << "%"
        << " max = " << max( data("err") ) * 100 << "%" << endl;
        /*! add Gaussian noise of error level to the data; */
        if ( addNoise ) {
            vcout << "Noisifying data by Gaussian error distribution" << endl;
            RVector noise( data.size() );
            randn( noise );  //! Standard normalized Gaussian distribution
            data.set( "t", data("t") * ( ( noise * data("err") ) + 1.0 ) ); //! a_noise = a * ( 1 + randn * da / a)
        }
        data.save( outFilename, std::string( "a m t err" ) ); // a=shot, m=geophone, t=time
    } else {
        data.save( outFilename, std::string( "a m t" ) ); // a=shot, m=geophone, t=time
    }

    return EXIT_SUCCESS;
}
