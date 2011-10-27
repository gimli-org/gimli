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
#include <datamap.h>

#include <inversion.h>
#include <longoptions.h>
#include <modellingbase.h>
#include <node.h>
#include <shape.h>

#include <string>
using namespace std;
using namespace GIMLI;

#define vcout if ( verbose ) cout
#define dcout if ( debug ) cout
#define DEBUG if ( debug )

#include <ttdijkstramodelling.h>

//** MAIN
int main( int argc, char *argv [] ){

    //** Initialisation of mesh options;
    int nSegments = 6;
    double relativeInnerMaxEdgeLength = 0.01;

    //** Initialisation of inversion options;
    bool lambdaOpt = false, isBlocky = false, isRobust = false, createGradientModel = false, useAppPar = false;
    bool isBertMesh = false;
    double lambda = 100.0, zPower = 0.0, zWeight = 1.0, errTime = 0.001, errPerc = 0.1, lbound = 0.0, ubound = 0.0;
    int iModel = 0, verboseCount = 0, maxIter = 20;
    string paraMeshFilename = NOT_DEFINED, startModelFilename = NOT_DEFINED, dataFileName;

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


  int option_char = 0, option_index = 0;

  if ( ( lbound > ubound ) && ( ubound > 0 ) ) {
    double dummy = lbound;
    lbound = 1.0 / ubound;
    ubound = 1.0 / dummy;
  }

  DEBUG verbose = true;
  if ( verbose ) setbuf(stdout, NULL);

  DataContainer dataIn( dataFileName );
  if ( verbose ) dataIn.showInfos();

  Mesh paraMesh;
  if ( paraMeshFilename != NOT_DEFINED ){
    paraMesh.load( paraMeshFilename );
  } else {
    cerr << WHERE_AM_I << " no para mesh given. Creating one. ..."  << endl;
    paraMesh.createClosedGeometryParaMesh( dataIn.electrodePositions(),
				       nSegments, relativeInnerMaxEdgeLength,
				       dataIn.additionalPoints() );
    DEBUG paraMesh.save( "meshPara.bms" );

  }
  vcout << "Paramesh: "; paraMesh.showInfos( );

  //** count unique attributes greater 0;
  std::set < int > paraAttributes;
  for ( size_t i = 0; i < paraMesh.cellCount(); i ++ ){
    if ( paraMesh.cell( i ).attribute() > 0 ) paraAttributes.insert( (int)paraMesh.cell( i ).attribute() );
  }

  //** if there are just 2 different atts(1-boundarydomain, 2-paradomain) we assume the loaded mesh has no para-attributes;
  if ( paraAttributes.size() == 2 ){
    if ( verbose ) std::cout << " create paramesh attributes. " << std::endl;
    //** so we will create them;
    paraMesh.createNeighbourInfos();
    int count = 1;
    for ( size_t i = 0; i < paraMesh.cellCount(); i ++ ){
      if ( paraMesh.cell( i ).attribute() == 2 ){
        paraMesh.cell( i ).setAttribute( ++count );
      }
    }
  }

  if ( paraAttributes.size() < 3 ){
    if ( verbose ) std::cout << " create paradomain mesh attributes. " << std::endl;
    paraMesh.createNeighbourInfos();
    int count = 1;
    for ( size_t i = 0; i < paraMesh.cellCount(); i ++ ) paraMesh.cell( i ).setAttribute( ++count );
  }

  Mesh paraDomain;
  vcout << "ParaDomain:"; paraDomain.showInfos( );

  Mesh secDomain;
  secDomain.createH2Mesh( paraMesh );

  //** put data error to 1% if not defined;
  if ( min( abs( dataIn.err() ) ) <= 0.0 ){
    vcout << "Estimate error: " << errPerc << "% + " << errTime << "s" << endl;
    dataIn.setErr( errTime / dataIn.t() + errPerc / 100.0 ); // always relative error
  }
  vcout << "Data error:"
	<< " min = " << min( dataIn.err() ) * 100 << "%"
	<< " max = " << max( dataIn.err() ) * 100 << "%" << endl;

  //** set up TT modeling class;
#ifdef WITH_OFFSET
  TravelTimeDijkstraModellingOffset f( secDomain, dataIn );
#else
  TravelTimeDijkstraModelling f( secDomain, dataIn );
#endif
  RVector appSlowness ( f.getApparentSlowness() );
  vcout << "min/max apparent velocity = " << 1.0 / max( appSlowness )
        << "/" << 1.0 / min( appSlowness ) << "m/s" << endl;

  size_t nModel = paraDomain.cellCount();
  vcout << "model size = " << nModel << endl;
  RVector startModel( f.startModel() );
  if ( nModel != startModel.size() ){
    cerr << "crosscheck WARNING! nMode !=  startModel.size() " << nModel << " != " << startModel.size() << std::endl;
  }

  if ( createGradientModel ) {
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
       startModel[i] = smi * exp( ( zmid[ i ] - zmi ) / ( zma - zmi ) * log( sma / smi) );
     }
     DEBUG save( startModel, "startModel.vector" );
  }
  if ( startModelFilename != NOT_DEFINED ){
    load( startModel, startModelFilename );
  }

  f.setStartModel( startModel );

  vcout << "Start model size=" << startModel.size() << " min/max=" << min(startModel) << "/" << max(startModel) << endl;
  vcout << "Creating Jacobian...";// << endl;
  //** create Jacobian matrix;

  DSparseMapMatrix S;
  f.createJacobian( S, startModel );
  vcout << "Ready. size = " << S.rows() << " x " << S.cols() << " " << S.size() << " elements." << endl;
  DEBUG save( S * startModel, "dataJac.vec" );
  DEBUG save( f( startModel ), "dataFor.vec" );
  DEBUG S.save( "S.matrix");
  //  RVector tmp( S * startModel / dataIn.t() -1 );
  //  save( tmp, "deltaData.vec" );
  //  exit( 0 );

  //** Data transformation: traveltime or apparent slowness (with bounds?);
  Trans < RVector > * transData;
  Trans < RVector > * transModel = new TransLogLU< RVector >( lbound, ubound );
  //Trans < RVector > * transModel = new Trans< RVector >( );
  if ( useAppPar ) {
      transData = new TransNest< RVector >( *new TransLogLU< RVector >( lbound, ubound ),
                                          *new TransLinear< RVector >( RVector( appSlowness / dataIn.t() ) ) );
  } else {
      transData = new Trans< RVector >( );
  }

#ifdef WITH_OFFSET
  //Inversion < TravelTimeDijkstraModellingOffset, RVector, SensMatType > inv( dataIn.t(), paraDomain, f, *transData, *transModel, verbose, debug );
  Inversion < double, DSparseMapMatrix > inv( dataIn.t(), paraDomain, f, *transData, *transModel, verbose, debug );
#else
  Inversion < double, DSparseMapMatrix > inv( dataIn.t(), paraDomain, f, *transData, *transModel, verbose, debug );
  //Inversion < TravelTimeDijkstraModelling, RVector, SensMatType > inv( dataIn.t(), paraDomain, f, *transData, *transModel, verbose, debug );
#endif

  inv.setModel( startModel, true );
  inv.setLambda( lambda );
  inv.setLineSearch( true );
  inv.setJacobian( S );
  inv.setError( dataIn.err() );
  inv.setOptimizeLambda( lambdaOpt );
  inv.setRobustData( isRobust );
  inv.setBlockyModel( isBlocky );
  inv.echoStatus( );

  int nBounds = paraDomain.boundaryCount();
  RVector flatWeight( nBounds );
  //* das folgende geht irgendwie nicht, deshalb per Hand
  RVector normDir( getParaBoundaryNormDir( paraDomain, 1 ) );
  for ( size_t i = 0; i < normDir.size(); i++ ) flatWeight[ i ] = std::pow( 1.0 - std::fabs( normDir[ i ] ), zPower ) + 0.1;

  inv.setCWeight( flatWeight );

  //** Compute resolution by solution of one inverse substep with sens column
  RVector resolution( inv.modelCellResolution( iModel ) );

  save( resolution, "resolution.vec" );
  double area( paraDomain.cell( iModel ).shape().domainSize() );
  double pdim( (double) paraDomain.dim() );
  double resrad( pow( area / resolution[ iModel ] * 3 / PI / ( 1.0 + pdim ), 1 / pdim ));
  std::cout << "Resolution at cell " << iModel << " = " << resolution[ iModel ] << " radius = " << resrad
            << " (domain size = " << area << ")" << std::endl;

  std::map< std::string, RVector > resultsToShow;
  resultsToShow.insert( make_pair( "resolution" , resolution ) );
  paraDomain.exportVTK( "ttinv.resolution", resultsToShow );

  delete transData;
  delete transModel;
  return EXIT_SUCCESS;
}

