/***************************************************************************
 *   Copyright (C) 2006-2010 by the resistivity.net development team       *
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
#include <optionmap.h>
#include <mesh.h>
#include <datacontainer.h>
#include <node.h>

#include <inversion.h>
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

    //!** Initialisation of mesh and inversion options;
    int nSegments = 6;
    double relativeInnerMaxEdgeLength = 0.01;
    bool lambdaOpt = false, isBlocky = false, isRobust = false, createGradientModel = false;
    bool useAppPar = false, isBertMesh = false, isVelocity = false;//, useWater = false;
    double errTime = 0.001, errPerc = 0.1, lbound = 0.0, ubound = 0.0;
    double lambda = 100.0, zWeight = 1.0;
    int maxIter = 20, verboseCount = 0;
    string paraMeshFilename = NOT_DEFINED, startModelFilename = NOT_DEFINED, dataFileName;

    //! Parse command line
    OptionMap oMap;
    oMap.setDescription("Description. TTInv - Travel time inversion using Dijkstra algorithm.\n");
    oMap.addLastArg( dataFileName, "Data file" );
    oMap.add( verboseCount,       "v" , "verbose", "Verbose mode (2 times for debug mode)." );
    oMap.add( isBertMesh,         "1" , "isBertMesh", "The mesh is a BERT mesh with 1 background" );
    oMap.add( useAppPar,          "A" , "UseAppPar", "Use Apparent Parameters." );
    oMap.add( isBlocky,           "B" , "BlockyModel", "Blocky (L1) model constraints." );
    oMap.add( createGradientModel,"G" , "gradientModel", "Gradient starting model." );
    oMap.add( isRobust,           "R" , "RobustData", "Robust (L1) data weighting." );
    oMap.add( lambdaOpt,          "O" , "OptimizeLambda", "Optimize model smoothness using L-curve." );
    oMap.add( isVelocity,         "V" , "isVelocity", "boundaries are velocities (not slowness)." );
    oMap.add( lambda,             "l:", "lambda", "Regularization parameter lambda." );
    oMap.add( maxIter,            "i:", "iterations", "Maximum iteration number." );
    oMap.add( errPerc,            "e:", "errorPerc", "Percentage error part." );
    oMap.add( errTime,            "t:", "errorTime", "Absolute error for first arrival times." );
    oMap.add( startModelFilename, "s:", "startModel", "Starting model." );
    oMap.add( zWeight,            "z:", "zWeight", "z weight (vertical smoothness)[1=isotropic]" );
    oMap.add( paraMeshFilename,   "p:", "paraMeshFile", "Parameter mesh file." );
    oMap.add( lbound,             "b:", "lowerBound", "Lower velocity/slowness bound." );
    oMap.add( ubound,             "u:", "upperBound", "Upper velocity/slowness bound(0=inactive)." );
    oMap.parse( argc, argv );

    bool verbose = ( verboseCount > 0 ), debug = ( verboseCount > 1 );
    if ( verbose ) setbuf(stdout, NULL);

    if ( ( ubound > 0. ) && ( lbound > ubound ) ) isVelocity = true; // apparently velocities are meant
    if ( isVelocity ) { // velocity specified instead of slowness
        double dummy = lbound;
        if ( ubound > 0.0 ) { // both bounds are given
            lbound = 1.0 / ubound;
            ubound = 1.0 / dummy;
        } else if ( lbound > 0.0 ) {
            lbound = 1.0 / dummy;
        }
    }
    vcout << "Lower/upper slowness bound = " << lbound << "/" << ubound << std::endl;

    //!** load data file
    DataContainer dataIn( dataFileName );
    if ( verbose ) dataIn.showInfos();

    //!** apply error model if not defined;
    if ( !dataIn.nonZero( "err" ) ) {
        vcout << "Estimate error: " << errPerc << "% + " << errTime << "s" << endl;
        dataIn.set( "err", errTime / dataIn( "t" ) + errPerc / 100.0 ); // always relative error
    }
    vcout << "Data error:" << " min = " << min( dataIn( "err" ) ) * 100 << "%" 	<< " max = " << max( dataIn( "err" ) ) * 100 << "%" << endl;

    //!** Load mesh
    Mesh paraMesh;
    if ( paraMeshFilename != NOT_DEFINED ) {
        paraMesh.load( paraMeshFilename );
    } else {
        cerr << " no para mesh given. Creating one. ..."  << endl;
        paraMesh.createClosedGeometryParaMesh( dataIn.sensorPositions(), //** better not
                                               nSegments, relativeInnerMaxEdgeLength,
                                               dataIn.additionalPoints() );
        DEBUG paraMesh.save( "meshPara.bms" );
    }
    vcout << "Paramesh: ";
    paraMesh.showInfos( );

    //!** set up TT modeling class;
    TravelTimeDijkstraModelling f( paraMesh, dataIn, false ); //**better with sec?
    RVector appSlowness ( f.getApparentSlowness() );
    vcout << "min/max apparent velocity = " << 1.0 / max( appSlowness ) << "/" << 1.0 / min( appSlowness ) << "m/s" << endl;

    //!** get mesh from region manager (BERT convention: first/smallest region is background
    if ( isBertMesh && f.regionManager().regionCount() > 1 ) f.regionManager().regions()->begin()->second->setBackground( true );
    Mesh paraDomain( f.regionManager().paraDomain() );
    vcout << "ParaDomain:\t";
    paraDomain.showInfos( );
    DEBUG paraDomain.save( "meshParaDomain.bms" );
    DEBUG paraDomain.exportVTK( "meshParaDomain" );

    //f.createRefinedForwardMesh();
    //vcout << "Secmesh:\t"; f->mesh()->showInfos( );
    //DEBUG f->mesh()->save( "meshSec.bms" );
    //DEBUG f->mesh()->exportVTK( "meshSec" );

    size_t nModel = f.regionManager().parameterCount();
    vcout << "model size = " << nModel << endl;
    RVector startModel( nModel );

    if ( createGradientModel ) { // auslagern in ttdijkstramodel
        vcout << "Creating Gradient model ..." << endl;
        double smi = min( appSlowness ) / 2.0;
        if ( smi < lbound ) smi = lbound * 1.1;
        double sma = max( appSlowness ) / 2.0;
        if ( ubound > 0.0 && sma > ubound ) sma = ubound * 0.9;

        RVector zmid( nModel);
        int di = paraDomain.dim() - 1;
        for ( size_t i = 0; i < paraDomain.cellCount() ; i++ ) {
            for ( size_t j = 0; j < paraDomain.cell( i ).nodeCount(); j++ ) {
                zmid[ i ] += paraDomain.cell( i ).node( j ).pos()[ di ];
            }
            zmid[ i ] /= double( paraDomain.cell( i ).nodeCount() );
        }
        double zmi = min( zmid );
        double zma = max( zmid );
        for ( size_t i = 0; i < startModel.size(); i++ ) {
            startModel[i] = smi * std::exp( ( zmid[ i ] - zmi ) / ( zma - zmi ) * std::log( sma / smi) );
        }
        DEBUG save( startModel, "startModel.vector" );
    }
    if ( startModelFilename != NOT_DEFINED ) {
        load( startModel, startModelFilename );
    }
    f.setStartModel( startModel );
    double sloRef = max( appSlowness );
    if ( sloRef < lbound || ( ubound > 0.0 && sloRef > ubound ) ) {
        sloRef = std::max( std::sqrt( lbound * ubound ), ( lbound + ubound ) / 2 );
    }

    vcout << "Start model size=" << startModel.size() << " min/max=" << min(startModel) << "/" << max(startModel) << endl;

    RTransLogLU transModel( lbound, ubound );
    RTrans transData;

    //!** Model transformation: log slowness with lower and upper bound
    Inversion < double, DSparseMapMatrix > inv( dataIn( "t" ), f, transData, transModel, verbose, debug ); 
    inv.setModel( startModel );
    inv.setLambda( lambda );
    inv.setMaxIter( maxIter );
    inv.setError( dataIn( "err" ) );
    inv.setOptimizeLambda( lambdaOpt );
    inv.setRobustData( isRobust );
    inv.setBlockyModel( isBlocky );

    vcout << "setting z weight...";
    f.regionManager().setZWeight( zWeight );
    vcout << "ready. Start inversion." << std::endl;
    //!** Start Inversion
    RVector model( inv.run() );
    vcout << "inversion ready." << std::endl;
    //!** Save velocity, model response and coverage
    RVector velocity( 1.0 / ( model + TOLERANCE ) );
    save( velocity, "velocity.vector" );
    save( inv.response(), "response.vector" );
    RVector one( dataIn.size(), 1.0 );
    RVector coverage( transMult( *inv.jacobian(), one ) );
    save( coverage, "coverage.vector" );
    RVector scoverage( transMult( inv.constraintMatrix(), inv.constraintMatrix() * coverage ) );
    for ( Index i ; i < scoverage.size() ; i++ ) scoverage[ i ] = sign( std::abs( scoverage[ i ] ) );
    save( scoverage, "scoverage.vector" );
    
    //!** Save VTK file with the vectors
    std::map< std::string, RVector > resultsToShow;
    resultsToShow.insert( make_pair( "slowness" , model ) );
    resultsToShow.insert( make_pair( "velocity" , velocity ) );
    resultsToShow.insert( make_pair( "coverage" , coverage ) );
    resultsToShow.insert( make_pair( "scoverage" , scoverage ) );
    paraDomain.exportVTK( "ttinv.result", resultsToShow );

    //!** Cleanup and exit
    return EXIT_SUCCESS;
}
