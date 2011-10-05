/***************************************************************************
 * This file belongs to GIMLi (Geophysical Inversion and Modelling) library*
 * and was created for documentation purposes only.                        *
 * DC1dBLOCK is a block inversion of resistivity soundings               *
 * Run with example file using: dc1dblock sond1-100dat                    *
***************************************************************************/

#include <gimli.h>
#include <dc1dmodelling.h>
#include <inversion.h>
#include <optionmap.h>
#include <modellingbase.h>

#include <string>
//using namespace std;
using namespace GIMLI;

#define vcout if ( verbose ) std::cout
#define dcout if ( debug ) std::cout
#define DEBUG if ( debug )

int main( int argc, char *argv [] )
{
    std::string dataFile( NOT_DEFINED );
    double errPerc = 3.0, lambda = 20.0, lbound = 0.0, ubound = 0.0;
    int nlay = 4;
    bool verbose = false, optimizeChi1 = false, doResolution = false, lambdaOpt = false;

    OptionMap oMap;
    oMap.setDescription("DC1dblock - block 1d dc resistivity inversion");
    oMap.add( verbose,      "v" , "verbose", "Be verbose." );
    oMap.add( lambdaOpt,    "O" , "OptimizeLambda", "Optimize model smoothness using L-curve." );
    oMap.add( optimizeChi1, "C" , "OptimizeChi1"  , "Optimize lambda subject to chi^2=1." );
    oMap.add( doResolution, "D" , "DoResolutionAnalysis", "Calculate resolution matrix." );
    oMap.add( lambda,       "l:", "lambda"        , "Regularization strength lambda." );
    oMap.add( nlay,         "n:", "nLayers"       , "Number of layers." );
    oMap.addLastArg( dataFile, "Datafile" );
    oMap.parse( argc, argv );

    /*! Data: read data file from column file */
    RMatrix abmnr; 
    loadMatrixCol( abmnr, dataFile ); //! read data
    RVector ab2(  abmnr[ 0 ] );        //! first column
    RVector mn2(  abmnr[ 1 ] );        //! second column
    RVector rhoa( abmnr[ 2 ] );        //! third column

    DC1dModelling f( nlay, ab2, mn2 );

    /*! Transformations: log for app. resisivity and thickness, logLU for resistivity */
    TransLog< RVector > transThk;
    TransLogLU< RVector > transRho( lbound, ubound );
    RTransLog transRhoa;
    
    f.region( 0 )->setTransModel( transThk );
    f.region( 1 )->setTransModel( transRho );

    double paraDepth = max( ab2 ) / 3;
    vcout << "Paradepth = " << paraDepth << std::endl;
    f.region( 0 )->setStartValue( paraDepth / nlay / 2.0 );
    f.region( 1 )->setStartValue( median( rhoa ) );

    RVector model;//( f.createDefaultStartModel() );
    model = f.createStartVector();
    model[ nlay ] *= 1.5;

    /*! Set up inversion with full matrix, data and forward operator */
    Inversion< double, RMatrix > inv( rhoa, f, verbose );
    inv.setTransData( transRhoa );
    inv.setLambda( lambda );
    inv.setRelativeError( errPerc / 100 );      //! error model
    inv.setModel( model );       //! starting model
    inv.setLocalRegularization( true ); //! Marquardt method
    inv.setOptimizeLambda( lambdaOpt );
    inv.stopAtChi1( false );
    inv.setMarquardtScheme( 0.9 );

    /*! actual computation: run the inversion */
    if ( optimizeChi1 ) {
        model = inv.runChi1( 0.1 );
    } else {
        model = inv.run();
    }
    std::cout << "model = " << model << std::endl;
    save( model, "model.vec" );

    if ( doResolution ) {
        size_t nModel( 2 * nlay - 1 );
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

