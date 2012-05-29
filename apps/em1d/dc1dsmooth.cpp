/***************************************************************************
 *   Copyright (C) 2008-2012 by the resistivity.net development team       *
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

#include <optionmap.h>
#include <inversion.h>
#include <dc1dmodelling.h>
//#include <expressions.h>
#include <string>

#define vcout if ( verbose ) std::cout

using namespace GIMLI;

int main( int argc, char *argv [] )
{
    std::string dataFile( NOT_DEFINED );
    double errPerc = 3.0, lambda = 20.0, lambdaIP = 1.0, maxDepth = 0.0;
    int nlay = 30, maxIter = 20;
    bool verbose = false, lambdaOpt = false, isBlocky = false, isRobust = false, optimizeChi1 = false;

    OptionMap oMap;
    oMap.setDescription("DC1dsmooth - smooth 1d dc resistivity inversion");
    oMap.add( verbose,      "v" , "verbose"       , "Verbose output." );
    oMap.add( lambdaOpt,    "O" , "OptimizeLambda", "Optimize model smoothness using L-curve." );
    oMap.add( isRobust,     "R" , "RobustData"    , "Robust (L1) data weighting." );
    oMap.add( isBlocky,     "B" , "BlockyModel"   , "Blocky (L1) model constraints." );
    oMap.add( optimizeChi1, "C" , "OptimizeChi1"  , "Optimize lambda subject to chi^2=1." );
    oMap.add( maxIter,      "i:", "maxIter"       , "Maximum iteration number." );
    oMap.add( lambda,       "l:", "lambda"        , "Regularization strength lambda." );
    oMap.add( lambdaIP,     "r:", "lambdaIP"      , "Regularization strength lambda for IP." );
    oMap.add( errPerc,      "e:", "error"         , "Error percentage" );
    oMap.add( nlay,         "n:", "nLay"          , "Number of layers" );
    oMap.add( maxDepth,     "d:", "maxDepth"      , "Maximum depth" );
    oMap.addLastArg( dataFile, "Datafile" );
    oMap.parse( argc, argv );

    /*! Data: read data file from column file */
    RMatrix abmnr; 
    loadMatrixCol( abmnr, dataFile ); //! read data
    RVector ab2(  abmnr[ 0 ] );        //! first column
    RVector mn2(  abmnr[ 1 ] );        //! second column
    RVector rhoa( abmnr[ 2 ] );        //! third column

    /*! Define discretization according to AB/2 */
    if ( maxDepth == 0.0 ) {
        maxDepth = max( ab2 ) / 2;   //! rule of thumb
        std::cout << "Maximum depth estimated to " << maxDepth << std::endl;
    }
    RVector thk( nlay - 1, maxDepth / ( nlay - 1 ) );
    thk *= ( maxDepth / sum( thk ) );
    RVector model( nlay, median( rhoa ) );

    /*! Transformations: log for app. resisivity and thickness, logLU for resistivity */
//    TransLogLU< RVector > transRho( lbound, ubound );
    TransLog< RVector > transRho;
    TransLog< RVector > transRhoa;

    /*! Forward operator, transformations and constraints */
    DC1dRhoModelling f( thk, ab2, mn2, false );
    f.region( 0 )->setTransModel( transRho );

    /*! Error estimation */
    double current = 0.1, errVolt = 0;//1e-4;
    RVector voltage( rhoa * f( RVector( 3, 1.0 ) ) * current );
    RVector error = errVolt / voltage + errPerc / 100.0;
    vcout << "error min/max = " << min( error ) << "/" << max( error ) << std::endl;

    /*! Set up inversion with full matrix, data and forward operator */
    RInversion inv( rhoa, f, verbose );
    inv.setTransData( transRhoa );
    inv.setRelativeError( error );
    inv.setLambda( lambda );
    inv.setModel( model );
    inv.setBlockyModel( isBlocky );
    inv.setRobustData( isRobust );
    inv.setOptimizeLambda( lambdaOpt );
    inv.setMaxIter( maxIter );
    inv.setDeltaPhiAbortPercent( 0.5 );
    inv.stopAtChi1( false );

    if ( optimizeChi1 ) {
        model = inv.runChi1( 0.1 );
    } else {
        model = inv.run();
    }
    if ( verbose ) std::cout << "model = " << model << std::endl;
    save( model, "resistivity.vec" );
    save( thk, "thickness.vec" );
    save( inv.response(), "response.vec" );
    
    if ( abmnr.cols() > 3 ) {
        if ( verbose ) std::cout << "Found ip values, doing ip inversion" << std::endl;
        //! imaginary apparent resistivity
        RVector rhoai( rhoa * sin( abmnr[ 3 ] / 1000 ) );
        //! linear modelling operator using the amplitude jacobian
        Mesh mesh( createMesh1D( thk.size() ) );
        LinearModelling fIP( mesh, f.jacobian(), verbose );
        fIP.region( 0 )->setTransModel( transRho );
        //! IP (imaginary resistivity) inversion using fIP
        RInversion invIP( rhoai, fIP, verbose );
        invIP.setAbsoluteError( rhoai / abmnr[ 3 ] * 1.0 );
        invIP.setBlockyModel( isBlocky );
        invIP.setRobustData( isRobust );
        invIP.setLambda( lambdaIP );
        invIP.setRecalcJacobian( false );
        invIP.stopAtChi1( false );
        RVector ipModel( nlay, median( rhoai ) );
        invIP.setModel( ipModel );
        ipModel = invIP.run();
        RVector phase = atan( ipModel / model ) * 1000.;
        save( phase, "phase.vec" );
        RVector aphase = atan( invIP.response() / rhoa ) * 1000.;
        save( aphase, "aphase.vec" ); 
    }

    return EXIT_SUCCESS;
}
