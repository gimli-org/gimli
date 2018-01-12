/***************************************************************************
 * This file belongs to GIMLi (Geophysical Inversion and Modelling) library*
 * and was created for documentation purposes only.                        *
 * DC1dSMOOTH is a smooth inversion of resistivity soundings               *
 * Run with example file using: dc1dsmooth sond1-100dat                    *
***************************************************************************/

#include <optionmap.h>
#include <inversion.h>
#include <dc1dmodelling.h>
#include <string>

using namespace GIMLI;

int main( int argc, char *argv [] )
{
    std::string dataFile( NOT_DEFINED );
    double errPerc = 3.0, lambda = 20.0;
    int nlay = 30, maxIter = 20;
    bool verbose = true, lambdaOpt = false, isBlocky = false, optimizeChi1 = false;

    OptionMap oMap;
    oMap.setDescription("DC1dsmooth - smooth 1d dc resistivity inversion");
    oMap.add( lambdaOpt,    "O" , "OptimizeLambda", "Optimize model smoothness using L-curve." );
    oMap.add( isBlocky,     "B" , "BlockyModel"   , "Blocky (L1) model constraints." );
    oMap.add( optimizeChi1, "C" , "OptimizeChi1"  , "Optimize lambda subject to chi^2=1." );
    oMap.add( lambda,       "l:", "lambda"        , "Regularization strength lambda." );
    oMap.add( maxIter,      "i:", "maxIter"       , "Maximum iteration number." );
    oMap.addLastArg( dataFile, "Datafile" );
    oMap.parse( argc, argv );

    /*! Data: read data file from column file */
    RMatrix abmnr; 
    loadMatrixCol( abmnr, dataFile ); //! read data
    RVector ab2(  abmnr[ 0 ] );        //! first column
    RVector mn2(  abmnr[ 1 ] );        //! second column
    RVector rhoa( abmnr[ 2 ] );        //! third column

    /*! Transformations: log for app. resisivity and thickness, logLU for resistivity */
//    RTransLogLU transRho( lbound, ubound );
    RTransLog transRho;
    RTransLog transRhoa;

    /*! Forward operator, transformations and constraints */
    RVector thk( nlay - 1, 0.5 );
    DC1dRhoModelling f( thk, ab2, mn2 );
    f.region( 0 )->setTransModel( transRho );
    //! variant 1: set marker in mesh
    f.mesh()->boundary( 8 ).setMarker( 1 );
    //! variant 2: set region boundary control
//    RVector bc( f.regionManager().constraintCount(), 1.0 );
//    bc[ 7 ] = 0.0;
//    f.region( 0 )->setConstraintsWeight( bc );

    /*! Set up inversion with full matrix, data and forward operator */
    RInversion inv( rhoa, f, transRhoa, transRho, verbose );
    inv.setRelativeError( errPerc / 100.0 );
    RVector model( nlay, median( rhoa ) ); //! starting model    
    inv.setModel( model );
    inv.setBlockyModel( isBlocky );
    inv.setLambda( lambda );
    inv.setOptimizeLambda( lambdaOpt );
    inv.setMaxIter( maxIter );

    //! variant 3: set inversion property
//    RVector bc( inv.constraintsCount(), 1.0 );
//    bc[ 6 ] = 0.0;
//    inv.setCWeight( bc ); 
    
    if ( optimizeChi1 ) {
        model = inv.runChi1( 0.1 );
    } else {
        model = inv.run();
    }
    if ( verbose ) std::cout << "model = " << model << std::endl;
    save( model, "resistivity.vec" );
    save( thk, "thickness.vec" );

    return EXIT_SUCCESS;
}
