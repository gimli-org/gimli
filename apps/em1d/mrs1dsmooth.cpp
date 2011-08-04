/***************************************************************************
 *   Copyright (C) 2009-2010 by the resistivity.net development team       *
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
#include <optionmap.h>
#include <em1dmodelling.h>
#include <inversion.h>

#include <string>
using namespace GIMLI;

#define vcout if ( verbose ) std::cout
#define dcout if ( debug ) std::cout
#define DEBUG if ( debug )

int main( int argc, char *argv [] ){

    bool lambdaOpt = false, isRobust = false, isBlocky = false, doResolution = false;
    bool useAppPar = false, useLog = false, optimizeChi1 = false;
    double lambda = 1, lbound = 0, ubound = 0, errorLevel = 20, zWeight = 1;
    int maxIter = 20, verboseCount = 0;//, ctype = 1;
    std::string matBasename( "K" ), dataFile( "b.vector" ), meshFile = NOT_DEFINED;

    OptionMap oMap;
    oMap.setDescription("Description. MRS1dSmooth - Smooth 1d/2d MRS inversion\n");
    oMap.addLastArg( dataFile, "Data file" );
    oMap.add( verboseCount, "v" , "verbose"       , "Verbose/debug/dosave mode (multiple use)." );
    oMap.add( lambdaOpt   , "O" , "OptimizeLambda", "Optimize model smoothness using L-curve." );
    oMap.add( optimizeChi1, "C" , "OptimizeChi1"  , "Optimize lambda subject to chi^2=1." );
    oMap.add( doResolution, "D" , "doResolution"  , "Do resolution analysis." );
    oMap.add( isRobust    , "R" , "RobustData"    , "Robust (L1) data weighting." );
    oMap.add( isBlocky    , "B" , "BlockyModel"   , "Blocky (L1) model constraints." );
    oMap.add( useLog      , "L" , "useLog"        , "Use (co-)Tan instead of Log for LU." );
    oMap.add( useAppPar   , "A" , "useAppPar"     , "Use apparent parameter transformation." );
    oMap.add( matBasename , "m:", "matBasename"   , "Matrix base name [K -> KR/KI.bmat]" );
    oMap.add( lambda      , "l:", "lambda"        , "Regularization strength lambda[100]." );
    oMap.add( lbound      , "b:", "lbound"        , "Lower parameter bound[0]" );
    oMap.add( ubound      , "u:", "ubound"        , "Upper parameter bound[0-inactive]" );
    oMap.add( errorLevel  , "e:", "error"         , "Absolute error level" );
    oMap.add( maxIter     , "i:", "maxIter"       , "Maximum Iteration number" );
    oMap.add( meshFile    , "p:", "meshFile"      , "Parameter mesh file." );
    oMap.add( zWeight     , "z:", "zWeight"       , "Weight for vertical smoothness (1=isotrope)." );
    oMap.parse( argc, argv );
    bool verbose = ( verboseCount > 0 ), debug = ( verboseCount > 1 ), dosave = ( verboseCount > 2 );

    RMatrix dataerror;
    loadMatrixCol( dataerror, dataFile );
    RVector data( dataerror[ 0 ] ); 
    std::cout << "Loaded " << data.size() << " data points." << std::endl;

    RVector error;
    if ( dataerror.cols() > 1 ) { // data errors present in 2nd column
        error = dataerror[ 1 ] / data;
        vcout << "loaded errors from file" << std::endl;
    } else { // no errors given -> constant noise level (-e)
        error = errorLevel / data;
    }
    vcout << "error: min/max = " << min( error ) * 100 << "/" << max( error ) * 100 << "%" << std::endl;
    
    // compute kernelfunction if possible, now just load it
    RMatrix KR, KI;
    vcout << "Trying to load " << matBasename << "R.bmat and " << matBasename << "I.bmat" << std::endl;
    if ( ! loadMatrixSingleBin( KR, matBasename + "R.bmat" ) ) { 
        std::cerr << "Did not find " << matBasename + "R.bmat" << "!" << std::endl; 
        return EXIT_OPEN_FILE; 
    }
    if ( ! loadMatrixSingleBin( KI, matBasename + "I.bmat" ) ) { 
        std::cerr << "Did not find " << matBasename + "I.bmat" << "!" << std::endl; 
        return EXIT_OPEN_FILE; 
    }

    size_t nModel( KR.cols() );
    if ( KR.rows() != data.size() ) {
        std::cerr << "Matrix (" << KR.rows() << " rows) mismatches data( " << data.size() << "points" << std::endl;
        return EXIT_MATRIX_SIZE_INVALID;
    }
    vcout << "size(K) = " << KR.rows() << "x" << nModel << std::endl;
    //! the parameterization
    Mesh mesh;
    if ( meshFile != NOT_DEFINED ) { //! real mesh given
        mesh.load( meshFile );
    } else {
        mesh = createMesh1D( nModel );
        mesh.save( "mesh1d.bms" );
//        for ( size_t i = 0; i < mesh.cellCount(); i ++ ) mesh.cell( i ).setAttribute( 2.0 + i );
    }
    mesh.createNeighbourInfos();
    if ( verbose ) mesh.showInfos();

    //! the forward operator
    MRSModelling f( mesh, KR, KI, verbose );
    //f.region( 0 )->setConstraintType( ctype );
    if ( zWeight > 0.0 ) f.regionManager().setZWeight( zWeight );
    //! apparent water content
    RVector K1( f.response( RVector( nModel, 1.0 ) ) );
    RVector wca( data / K1 );
    save( wca, "wca.vec" );
    vcout << "apparent water content: min/max = " << min( wca ) << "/" << max( wca ) << std::endl;
    
    //! transformations (variable as pointer)
    Trans < RVector > * transData;
    Trans < RVector > * transModel;
    //! model transformation
    if ( useLog ) {
            transModel = new TransLogLU< RVector >( lbound, ubound );
    } else {
            if ( ubound <= lbound ) ubound = lbound + 1.0;
            transModel = new TransCotLU< RVector >( lbound, ubound );
    }
    //! data transformation
    if ( useAppPar ) {
        if ( useLog ) { //! logarithm of water content
            transData = new TransNest< RVector >( *new TransLogLU< RVector >( lbound, ubound ),
                                                  *new TransLinear< RVector >( RVector( wca / data ) ) );
        } else { //! cotangens bound trans of apparent water content
            transData = new TransNest< RVector >( *new TransCotLU< RVector >( lbound, ubound ),
                                                  *new TransLinear< RVector >( RVector( wca / data ) ) );
        }
    } else { //! linear (default)
        transData = new Trans< RVector >( );
    }
    
    //! the inversion class
    Inversion< double, RMatrix > inv( data, f, *transData, *transModel, verbose, dosave );
    inv.setLambda( lambda );
    inv.setOptimizeLambda( lambdaOpt );
    inv.setRobustData( isRobust );
    inv.setBlockyModel( isBlocky );
    inv.setMaxIter( maxIter );
    inv.setDeltaPhiAbortPercent( 0.1 );
    inv.setError( error );
    //! starting model
    RVector wc( nModel, mean( wca ) );
    inv.setModel( wc );
    inv.setReferenceModel( wc );
    
    //! inversion run
    if ( optimizeChi1 ) { 
        wc = inv.runChi1(); 
    } else { 
        wc = inv.run(); 
    }
    if ( verbose && wc.size() < 200 ) std::cout << "wc = " << wc << std::endl;
    save( wc, "wc.vec" );
    save( inv.response(), "response.vec" );
    
    //! resolution analysis
    if ( doResolution ) {
        RVector resolution( nModel );
        RVector resMDiag ( nModel );
        RMatrix resM;
        for ( size_t iModel = 0; iModel < nModel; iModel++ ) {
            resolution = inv.modelCellResolution( iModel );
            resM.push_back( resolution );
            resMDiag[ iModel ] = resolution[ iModel ];
        }
        save( resMDiag, "resMDiag.vec" );
        save( resM, "resolution.matrix" );
    }

    delete transData;
    delete transModel;

    return EXIT_SUCCESS;
}

