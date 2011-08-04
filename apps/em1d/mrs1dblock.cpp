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

    bool lambdaOpt = false, doResolution = false, optimizeChi1 = false;
    double lambda = 1, lbound = 0, ubound = 1.0, errorLevel = 20;
    int maxIter = 20, verboseCount = 0, nlay = 3;
    std::string matBasename( "K" ), dataFile( "b.vector" ), startModelFile = NOT_DEFINED;

    OptionMap oMap;
    oMap.setDescription("Description. MRS1dBlock - MRS 1d block inversion\n");
    oMap.addLastArg( dataFile, "Data file" );
    oMap.add( verboseCount, "v" , "verbose"       , "Verbose/debug/dosave mode (multiple use)." );
    oMap.add( lambdaOpt   , "O" , "OptimizeLambda", "Optimize model smoothness using L-curve." );
    oMap.add( optimizeChi1, "C" , "OptimizeChi1"  , "Optimize lambda subject to chi^2=1." );
    oMap.add( doResolution, "D" , "doResolution"  , "Do resolution analysis." );
    oMap.add( matBasename , "m:", "matBasename"   , "Matrix base name [K -> KR/KI.bmat]" );
    oMap.add( nlay,         "n:", "nlay", "Number of layers" );
    oMap.add( lambda      , "l:", "lambda"        , "Regularization strength lambda[100]." );
    oMap.add( lbound      , "b:", "lbound"        , "Lower parameter bound[0]" );
    oMap.add( ubound      , "u:", "ubound"        , "Upper parameter bound[0-inactive]" );
    oMap.add( errorLevel  , "e:", "error"         , "Absolute error level" );
    oMap.add( maxIter     , "i:", "maxIter"       , "Maximum Iteration number" );
    oMap.add( startModelFile, "s:", "startModel"  , "Start model file" );
    oMap.parse( argc, argv );
    bool verbose = ( verboseCount > 0 ), debug = ( verboseCount > 1 ), dosave = ( verboseCount > 2 );

    RMatrix dataerror;
    loadMatrixCol( dataerror, dataFile );
    RVector data( dataerror[ 0 ] ); 
    std::cout << "Loaded " << data.size() << " data points." << std::endl;
        
    RVector error;
    if ( dataerror.rows() > 1 ) { // data errors present in 2nd column
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
    RVector zvec( "zkernel.vec" );

    size_t nModel( KR.cols() );
    if ( KR.rows() != data.size() ) {
        std::cerr << "Matrix (" << KR.rows() << " rows) mismatches data( " << data.size() << "points" << std::endl;
        return EXIT_MATRIX_SIZE_INVALID;
    }
    vcout << "size(K) = " << KR.rows() << "x" << nModel << std::endl;
    vcout << "size(z) = " << zvec.size() << std::endl;

    //! start model
    RVector model;
    if ( startModelFile != NOT_DEFINED ) {
        vcout << "Loading starting model from file " << startModelFile << std::endl;
        load( model, startModelFile );
        nlay = size_t( ( model.size() + 1.0 ) / 2 );
        vcout << nlay << " layers" << std::endl;
    }

    //! the parameterization
    Mesh mesh = createMesh1DBlock( nlay );
    if ( verbose ) mesh.showInfos();

    //! the forward operator
    MRS1dBlockModelling f( mesh, KR, KI, zvec, debug );

    //! apparent water content
    RVector K1( f.response( RVector( nlay * 2 - 1, 1.0 ) ) );
    save( K1, "K1.vec" );
    RVector wca( data / K1 );
    save( wca, "wca.vec" );
    vcout << "apparent water content: min/max = " << min( wca ) << "/" << max( wca ) << std::endl;
    
    //! transformations
    RTrans transData;
    RTransLogLU transWC( lbound, ubound );
    RTransLog transThk;
    f.region(0)->setTransModel( transThk );
    f.region(1)->setTransModel( transWC );
    
    //! zeroth order constraints, starting model
    f.region(0)->setStartValue( 5.0 );
    f.region(1)->setStartValue( median( wca ) );

    if ( startModelFile == NOT_DEFINED ) model = f.createStartVector();
    if ( debug ) {
        save( model, "startmodel.vec" );
        save( f( model ), "startresponse.vec" );
    }
    

    //! the inversion class
    Inversion< double, RMatrix > inv( data, f, verbose, dosave );
    inv.setTransData( transData );
    inv.setMarquardtScheme( 0.8 ); //! 0th order local decreasing regularization
    inv.setLambda( lambda );
    inv.setOptimizeLambda( lambdaOpt );
    inv.setMaxIter( maxIter );
    inv.setDeltaPhiAbortPercent( 0.2 );
    inv.setError( error );
    inv.setModel( model );
    if ( startModelFile != NOT_DEFINED ) inv.setReferenceModel( model );
    
    //! inversion run
    if ( optimizeChi1 ) { 
        model = inv.runChi1(); 
    } else { 
        model = inv.run(); 
    }
    save( model, "model.vec" );
    save( inv.response(), "response.vec" );
    vcout << "wc=";
    vcout << " " << model( nlay - 1, nlay * 2 - 1 ) << std::endl;
    vcout << "thk =  ";
    vcout << " " << model( 0 , nlay - 1 ) << std::endl;
    
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

    return EXIT_SUCCESS;
}

