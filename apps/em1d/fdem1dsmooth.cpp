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
#include <em1dmodelling.h>
#include <inversion.h>

#include <string>
//using namespace std;
using namespace GIMLI;

#define vcout if ( verbose ) std::cout
#define dcout if ( debug ) std::cout
#define DEBUG if ( debug )

int main( int argc, char *argv [] )
{
    bool lambdaOpt = false, doResolution = false, useTan = false, optimizeChi1 = false, delinhigh = true;
    double lambda = 10.0, lbound = 0.0, ubound = 0.0, errPerc = 3.0, coilspacing = 50.0, depth = 30.0;
    int maxIter = 10, nlay = 30, verboseCount = 0;
    std::string modelFile( NOT_DEFINED ), dataFileName;

    OptionMap oMap;
    oMap.setDescription("Description. DC1dInv - 1D block inversion of dc resistivity data\n");
    oMap.addLastArg( dataFileName, "Data file" );
    oMap.add( verboseCount, "v" , "verbose", "Verbose/debug/dosave mode (multiple use)." );
    oMap.add( lambdaOpt,    "O" , "OptimizeLambda", "Optimize model smoothness using L-curve." );
    oMap.add( doResolution, "D" , "doResolution", "Do resolution analysis." );
    oMap.add( optimizeChi1, "C" , "OptimizeChi1", "Optimize lambda subject to chi^2=1." );
    oMap.add( useTan,       "T" , "useTan", "Use (co-)Tan instead of Log for LU." );
    oMap.add( maxIter,      "i:", "maxIter", "Maximum iteration number." );
    oMap.add( lambda,       "l:", "lambda", "Regularization strength lambda." );
    oMap.add( coilspacing,  "c:", "coilSpacing", "Coil spacing in m." );
    oMap.add( depth,        "d:", "depth", "Depth of model." );
    oMap.add( errPerc,      "e:", "error", "Error in percent." );
    oMap.add( nlay,         "n:", "nlay", "Number of layers." );
    oMap.add( lbound,       "b:", "lbound", "Lower Resistivity bound." );
    oMap.add( ubound,       "u:", "ubound", "Upper Resistivity bound." );
    oMap.add( modelFile,    "m:", "modelFile", "Model file for pure modelling. file: thick_i rho_i\\n" );
    oMap.parse( argc, argv );

    bool verbose = ( verboseCount > 0 ), debug = ( verboseCount > 1 ), dosave = ( verboseCount > 2 );
    
    RMatrix DATA;
    loadMatrixCol( DATA, dataFileName );
    RVector freq( DATA[0] );
    size_t nfreq = freq.size(), ndata = nfreq * 2;
    RVector data( cat(DATA[1], DATA[2]) ); //** inphase and quadrature
    RVector thk( nlay, depth / nlay );

    FDEM1dRhoModelling f( thk, freq, coilspacing, 1.0, false );
    DEBUG save( f.freeAirSolution(), "freeairsolution.vec" );
    
    /*! Error estimation */
    RVector error( ndata, errPerc );
    if ( delinhigh ) error[ nfreq - 1 ] = errPerc * 10.0;   //** strange reading for highest f inphase

    /*! Transformations: log for app. resisivity and thickness, logLU for resistivity */
//    RTransLogLU transThk( 0.1, 1e2 );
    RTransLog transThk;
    RTransLogLU transRho( lbound, ubound );    
//    RTransLog transRho;
    f.region( 0 )->setTransModel( transThk );
    f.region( 1 )->setTransModel( transRho );
    f.region( 0 )->setStartValue( 10.0 );
    f.region( 1 )->setStartValue( 10.0 );

    RVector model, response;
    model = f.createStartVector();
    DEBUG save( model, "startModel.vec" );
    response = f(model);
    DEBUG save( response, "startResponse.vec" );

    /*! Set up inversion with full matrix, data and forward operator */
    RInversion inv( data, f, verbose, dosave );    
    inv.setMarquardtScheme( 0.8 ); //! 0th order local decreasing regularization
    inv.setLambda( lambda );
    inv.setOptimizeLambda( lambdaOpt );
    inv.setMaxIter( maxIter );
    inv.setAbsoluteError( error );              //! error model
    inv.setModel( model );        //! starting model
    inv.setDeltaPhiAbortPercent( 0.5 );

    /*! actual computation: run the inversion */
    model = inv.run();

    save( model, "resistivity.vec" );
    save( thk, "thickness.vec" );
    save( inv.response(), "response.vec" );
    
    if ( verbose ) {
        RVector cumz( thk );
        for ( size_t i = 1 ; i < thk.size() ; i++ ) cumz[i] = cumz[ i-1 ] + thk[i];
		cumz.round(0.1);        
		model.round(0.1);        
        std::cout << "Res = " << model << std::endl;
        std::cout << "  z =  " << cumz << std::endl;
    }

    if ( doResolution ) {
        RVector resolution( nlay );
        RVector resMDiag ( nlay );
        RMatrix resM;
        for ( int iModel = 0; iModel < nlay; iModel++ ) {
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

