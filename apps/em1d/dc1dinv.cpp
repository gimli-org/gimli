/***************************************************************************
 *   Copyright (C) 2008-2011 by the resistivity.net development team       *
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
#include <dc1dmodelling.h>
#include <meshgenerators.h>
#include <mesh.h>
#include <inversion.h>
#include <modellingbase.h>
#include <solver.h>

#include <string>
//using namespace std;
using namespace GIMLI;

#define vcout if ( verbose ) std::cout
#define dcout if ( debug ) std::cout
#define DEBUG if ( debug )

int main( int argc, char *argv [] )
{
    bool lambdaOpt = false, doResolution = false, useTan = false, optimizeChi1 = false;
    double lambda = 10.0, lambdaIP = 1.0, lbound = 0.0, ubound = 0.0, errPerc = 3.0;
    int maxIter = 10, nlay = 3, verboseCount = 0, linCount = 0;
    std::string modelFile( NOT_DEFINED ), dataFileName;

    OptionMap oMap;
    oMap.setDescription("Description. DC1dInv - 1D block inversion of dc resistivity data\n");
    oMap.addLastArg( dataFileName, "Data file" );
    oMap.add( verboseCount, "v" , "verbose", "Verbose/debug/dosave mode (multiple use)." );
    oMap.add( lambdaOpt,    "O" , "OptimizeLambda", "Optimize model smoothness using L-curve." );
    oMap.add( doResolution, "D" , "doResolution", "Do resolution analysis." );
    oMap.add( optimizeChi1, "C" , "OptimizeChi1", "Optimize lambda subject to chi^2=1." );
    oMap.add( useTan,       "T" , "useTan", "Use (co-)Tan instead of Log for LU." );
    oMap.add( linCount,     "L" , "dataLin", "Use linear trafo for data, 2x also for model." );
    oMap.add( maxIter,      "i:", "maxIter", "Maximum iteration number." );
    oMap.add( lambda,       "l:", "lambda", "Regularization strength lambda." );
    oMap.add( lambdaIP,     "r:", "lambdaIP", "Regularization strength lambda for IP." );
    oMap.add( errPerc,      "e:", "error", "Error percentage" );
    oMap.add( nlay,         "n:", "nlay", "Number of layers" );
    oMap.add( lbound,       "b:", "lbound", "Lower Resistivity bound" );
    oMap.add( ubound,       "u:", "ubound", "Upper Resistivity bound" );
    oMap.add( modelFile,    "m:", "modelFile", "Model file for pure modelling. file: thick_i rho_i\\n" );
    oMap.parse( argc, argv );

    bool verbose = ( verboseCount > 0 ), debug = ( verboseCount > 1 ), dosave = ( verboseCount > 2 );
    bool dataLin = ( linCount > 0 ), modelLin = ( linCount > 1 );
    
    DataContainer data( dataFileName ); //! read data
    if ( verbose ) data.showInfos();

    if ( modelFile != NOT_DEFINED ) { //! pure forward modeling
        calculateDC1D( data, modelFile, dataFileName + ".sim" );
        return EXIT_SUCCESS;
    }

    if ( data.nonZero( "rhoa" ) ) {
        data.set( "rhoa", data("r") * geometricFactor( data ) );
    }

    size_t nModel( 2 * nlay - 1 );
    Mesh mesh( createMesh1DBlock( nlay ) );
    if ( debug ) mesh.showInfos();

    DC1dModelling f( mesh, data, nlay, false );

    /*! Error estimation */
    RVector error( data( "err" ) );
        if ( std::fabs( min( error ) ) < TOLERANCE ) {
        double current = 0.1, errVolt = 0;//1e-4;
        RVector voltage( data( "rhoa" ) / data( "k" ) * current );
        error = errVolt / voltage + errPerc / 100.0;
    }
    vcout << "error min/max = " << min( error ) << "/" << max( error ) << std::endl;

    /*! Transformations: log for app. resisivity and thickness, logLU for resistivity */
//    TransLog< RVector > transThk;
//    TransLogLU< RVector > transRho( lbound, ubound );
//    TransLog< RVector > transRhoa;
    Trans< RVector > *transThk, *transRho, *transRhoa;
    if ( modelLin ) {
         transThk = new Trans< RVector >;
         transRho = new Trans< RVector >;
    } else {
        transThk = new TransLog< RVector >;
        if ( useTan ) {
            transRho = new TransCotLU< RVector >( lbound, ubound );
        } else {
            transRho = new TransLogLU< RVector >( lbound, ubound );
        }
    }
    if ( dataLin ) {
        transRhoa = new Trans< RVector >;
    } else {
        transRhoa = new TransLog< RVector >;
    }

    f.region( 0 )->setTransModel( *transThk );
    f.region( 1 )->setTransModel( *transRho );

    double paraDepth = DCParaDepth( data );
    vcout << "Paradepth = " << paraDepth << std::endl;
    f.region( 0 )->setStartValue( paraDepth / nlay / 4.0 );
    f.region( 1 )->setStartValue( median( data( "rhoa" ) ) );

    RVector model;//( f.createDefaultStartModel() );
    model = f.createStartVector();
//    model[ nlay ] *= 1.5; //! in order to make thicknesses sensitive
    DEBUG save( model, "start.vec" );

    /*! Set up inversion with full matrix, data and forward operator */
    Inversion< double, RMatrix > inv( data( "rhoa" ), f, verbose, dosave );    
    inv.setTransData( *transRhoa );
    inv.setMarquardtScheme( 0.8 ); //! 0th order local decreasing regularization
    inv.setLambda( lambda );
    inv.setOptimizeLambda( lambdaOpt );
    inv.setMaxIter( maxIter );
    inv.setError( error );              //! error model
    inv.setModel( model );       //! starting model
    inv.setDeltaPhiAbortPercent( 0.5 );

    /*! actual computation: run the inversion */
    if ( optimizeChi1 ) {
        model = inv.runChi1( 0.1 );
    } else {
        model = inv.run();
    }

    RVector thk( nlay - 1 );
    RVector res( nlay );
    for ( int i = 0 ; i < nlay - 1 ; i++ ) thk[ i ] = model[ i ];
    for ( int i = 0 ; i < nlay ; i++ ) res[ i ] = model[ nlay - 1 + i ];
    save( res, "resistivity.vec" );
    save( thk, "thickness.vec" );
    save( inv.response(), "response.vec" );
    
    if ( verbose ) {
        RVector cumz( thk );
        for ( size_t i = 1 ; i < thk.size() ; i++ ) cumz[i] = cumz[ i-1 ] + thk[i];
        for ( size_t i = 0 ; i < cumz.size() ; i++ ) cumz[i] = round( cumz[i]*10.0 ) / 10.0;        
        for ( size_t i = 0 ; i < res.size() ; i++ ) res[i] = round( res[i] );        
        std::cout << "Res = " << res << std::endl;
        std::cout << "  z =  " << cumz << std::endl;
    }

    if ( data.nonZero( "ip" ) ) {
//        if ( verbose ) std::cout << "Found ip values, doing ip inversion" << std::endl;
        //! imaginary apparent resistivity
        RVector rhoai( data( "rhoa" ) * sin( data( "ip" ) / 1000 ) );
        //! modelling operator from fixed layer 
        Mesh ipMesh( createMesh1D( nlay ) );
        DC1dRhoModelling fIP( ipMesh, data, thk, verbose );
        fIP.region( 0 )->setTransModel( *transRho );
        //! IP (imaginary resistivity) inversion using fIP
        RInversion invIP( rhoai, fIP, verbose );
        if ( min( data.get( "iperr" ) ) > 0.0 ) { //! phase error present
            invIP.setRelativeError( data.get( "iperr" ) / data( "ip" ) );
        } else { //! default: 1 mrad
            invIP.setAbsoluteError( rhoai / data( "ip" ) * 1.0 );
        }
        //! take jacobian from real valued problems
        invIP.setModel( res );
        invIP.checkJacobian(); //! force jacobian calculation
        invIP.setRecalcJacobian( false );
        
        invIP.setLambda( lambdaIP );
        invIP.setLocalRegularization( true ); //! Marquardt method
        invIP.setModel( model );       //! starting model
        invIP.stopAtChi1( false );
        RVector ipModel( nlay, median( rhoai ) );
        invIP.setModel( ipModel );
        ipModel = invIP.run();
        RVector phase = atan( ipModel / res ) * 1000;
        save( phase, "phase.vec" );
        RVector aphase = atan( invIP.response() / data( "rhoa" ) ) * 1000;
        save( aphase, "aphase.vec" );         
        std::cout << " ip =  " << phase << std::endl;
    }
    
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
        save( resM,     "resM" );
        vcout << "diagRM = " << resMDiag << std::endl;
    }

    return EXIT_SUCCESS;
}

