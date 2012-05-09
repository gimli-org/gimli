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
#include <meshgenerators.h>
#include <mesh.h>
#include <inversion.h>
#include <optionmap.h>
#include <modellingbase.h>
//#include <trans.h>

#include <string>
//using namespace std;
using namespace GIMLI;

#define vcout if ( verbose ) std::cout
#define dcout if ( debug ) std::cout
#define DEBUG if ( debug )

int main( int argc, char *argv [] ){

    bool lambdaOpt = false, isRobust = false, isBlocky = false, doResolution = false;
    bool useAppPar = false, useTan = false, useWater = false, optimizeChi1 = false;
    double lambda = 1, lbound = 0, ubound = 0, errbabs = 20;
    int maxIter = 10, verboseCount = 0, ctype = 0;
    std::string matFile( "A.mat" ), bFile( "b.vec" );

    OptionMap oMap;
    oMap.setDescription("Description. InvLinearMat - Linear Inversion with given matrix and vector\n");
    oMap.addLastArg( bFile, "Data file" );
    oMap.add( verboseCount, "v" , "verbose", "Verbose/debug/dosave mode (multiple use)." );
    oMap.add( lambdaOpt,    "O" , "OptimizeLambda", "Optimize model smoothness using L-curve." );
    oMap.add( optimizeChi1, "C" , "OptimizeChi1", "Optimize lambda subject to chi^2=1." );
    oMap.add( doResolution, "D" , "doResolution", "Do resolution analysis." );
    oMap.add( isRobust,     "R" , "RobustData", "Robust (L1) data weighting." );
    oMap.add( isBlocky,     "B" , "BlockyModel", "Blocky (L1) model constraints." );
    oMap.add( useTan,       "T" , "useTan", "Use (co-)Tan instead of Log for LU." );
    oMap.add( useAppPar,    "A" , "useAppPar", "Use apparent parameter transformation." );
    oMap.add( useWater,     "W" , "useWater", "Use water content transformation." );
    oMap.add( matFile,      "m:", "matFile", "Matrix file [A.mat]" );
    oMap.add( lambda,       "l:", "lambda", "Regularization strength lambda[100]." );
    oMap.add( lbound,       "b:", "lbound", "Lower parameter bound[0]" );
    oMap.add( ubound,       "u:", "ubound", "Upper parameter bound[0-inactive]" );
    oMap.add( errbabs,      "e:", "error", "Absolute error level" );
    oMap.add( maxIter,      "i:", "maxIter", "Maximum Iteration number" );
    oMap.parse( argc, argv );
    bool verbose = ( verboseCount > 0 ), debug = ( verboseCount > 1 ), dosave = ( verboseCount > 2 );

    RMatrix A;
    if ( ! loadMatrixSingleBin( A, matFile ) ) { std::cerr << "Did not find A.mat!" << std::endl; return EXIT_OPEN_FILE; }
    RVector b; load( b, bFile );
    size_t nModel( A.cols() );
    dcout << "size(A) = " << A.rows() << "x" << nModel << "size(b) = " << b.size() << std::endl;

    RVector Asum( A * RVector( nModel, 1.0 ) );
    RVector xapp( b / Asum );
    DEBUG save( xapp, "xapp.vec" );
    dcout << "apparent x: min/max = " << min( xapp ) << "/" << max( xapp ) << std::endl;

    Mesh mesh( createMesh1D( nModel ) );
    DEBUG mesh.save( "mesh1d.bms" );
    mesh.showInfos();
    for ( size_t i = 0; i < mesh.cellCount(); i ++ ) mesh.cell( i ).setAttribute( 2.0 + i );
    mesh.createNeighbourInfos();

    LinearModelling f( mesh, A, verbose );
    f.regionManager().region( 0 )->setConstraintType( ctype );

    Trans < RVector > * transData;
    Trans < RVector > * transModel;
    if ( useAppPar ) {
        if ( useTan ) {
            if ( ubound <= lbound ) ubound = lbound + 1.0;
            transData = new TransNest< RVector >( *new TransCotLU< RVector >( lbound, ubound ),
                                                  *new TransLinear< RVector >( RVector( xapp / b ) ) );
            transModel = new TransCotLU< RVector >( lbound, ubound );
        } else {
            transData = new TransNest< RVector >( *new TransLogLU< RVector >( lbound, ubound ),
                                                  *new TransLinear< RVector >( RVector( xapp / b ) ) );
            transModel = new TransLogLU< RVector >( lbound, ubound );
        }
    } else {
        transData = new Trans< RVector >( );
        transModel = new Trans< RVector >( );
    }
    /*! set up inversion */
    RInversion inv( b, f, *transData, *transModel, verbose, dosave );

    inv.setRecalcJacobian( false ); //! no need since it is linear
    inv.setLambda( lambda );
    inv.setOptimizeLambda( lambdaOpt );
    inv.setRobustData( isRobust );
    inv.setBlockyModel( isBlocky );
    inv.setMaxIter( maxIter );
    inv.setDeltaPhiAbortPercent( 1.0 );
    RVector error( errbabs / b);
    vcout << "error: min/max = " << min( error ) << "/" << max( error ) << std::endl;
    inv.setError( error );

    RVector x( nModel, mean( xapp ) );
    inv.setModel( x );
    inv.setReferenceModel( x );

    if ( optimizeChi1 ) { x = inv.runChi1(); } else { x = inv.run(); }
    vcout << "x = " << x << std::endl;
    save( x, "x.vec" );

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

