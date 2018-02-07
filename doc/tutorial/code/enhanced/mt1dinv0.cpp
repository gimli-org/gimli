/***************************************************************************
 *   Copyright (C) 2008-2009 by the resistivity.net development team       *
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
#include <gimli.h>
#include <em1dmodelling.h>
#include <inversion.h>
#include <complex.h>
using namespace GIMLI;

#define vcout if ( verbose ) std::cout
#define dcout if ( debug ) std::cout
#define DEBUG if ( debug )

int main( int argc, char *argv [] ) {
    bool verbose = true, lambdaOpt = false, doResolution = false;
    double lambda = 10.0, lbound = 0.0, ubound = 0.0, errRhoa = 5.0*1000, errPhase = 1.0*PI/180;
    int maxIter = 10, nlay = 3;
    std::string modelFile( NOT_DEFINED ), dataFileName;

    OptionMap oMap;
    oMap.setDescription("Description. MT1dInv - 1D block or smooth inversion of dc resistivity data\n");
    oMap.addLastArg( dataFileName, "Data file" );
    oMap.add( verbose,      "v" , "verbose", "Verbose output." );
    oMap.add( lambdaOpt,    "O" , "OptimizeLambda", "Optimize model smoothness using L-curve." );
    oMap.add( doResolution, "D" , "doResolution", "Do resolution analysis." );
    oMap.add( lambda,       "l:", "lambda", "Regularization strength lambda." );
    oMap.add( errRhoa,      "e:", "error", "Relative error for rhoa in percent" );
    oMap.add( nlay,         "n:", "nlay", "Number of layers" );
    oMap.add( lbound,       "b:", "lbound", "Lower Resistivity bound" );
    oMap.add( ubound,       "u:", "ubound", "Upper Resistivity bound" );
    oMap.parse( argc, argv );

    RMatrix TRP; loadMatrixCol( TRP, dataFileName );
    size_t nP = TRP.cols();
    RVector T( TRP[ 0 ] ), rhoaMT( TRP[ 1 ] ), phiMT( TRP[ 2 ] ); //! periods/app.res./phases
    double medrhoa = median( rhoaMT );

    /*! Model transformations: log for resistivity and thickness */
    TransLog< RVector > transThk;
    TransLogLU< RVector > transRho( lbound, ubound );
    /*! Data transformations: log apparent resistivity, linear phases */
    TransLog< RVector > transRhoa;
    Trans< RVector > transPhi;
    CumulativeTrans< RVector > transData; //! combination of two trans functions
    transData.push_back( transRhoa, nP ); //! append rhoa trans
    transData.push_back( transPhi, nP );  //! append phi trans

    /*! Modelling operator and constraints */
    MT1dModelling f( T, nlay, false );
    f.regionManager().setConstraintType( 0 ); //! minimum length (no smoothing) for all
    f.region( 0 )->setTransModel( transThk );
    f.region( 1 )->setTransModel( transRho );
    double medskindepth = sqrt( median( T ) * medrhoa ) *503.0;
    f.region( 0 )->setStartValue( medskindepth / nlay );
    f.region( 1 )->setStartValue( medrhoa );

    /*! Set up inversion with full matrix, data and forward operator */
    RInversion inv( cat( rhoaMT, phiMT ), f, verbose ); //! only rhoa
    inv.setTransData( transData );
    inv.setLambda( lambda );
    inv.setOptimizeLambda( lambdaOpt );
    inv.setMaxIter( maxIter );
    inv.setLocalRegularization( true ); //! Marquardt method

    /*! Error model combined of rhoa error and phi error */
    RVector error( cat( RVector( nP, errRhoa / 100.0 ), RVector( errPhase / phiMT ) ) );
    inv.setRelativeError( error );              //! error model

    /*! actual computation: run the inversion */
    RVector model = inv.run();
    std::cout << "model = " << model << std::endl;
    save( model, "model.vec" );
    save( inv.response(), "response.vec" );

    if ( doResolution ) {
        RVector resolution( model.size() );
        RVector resMDiag ( model.size() );
        RMatrix resM;
        for ( size_t iModel = 0; iModel < model.size(); iModel++ ) {
            resolution = inv.modelCellResolution( iModel );
            resM.push_back( resolution );
            resMDiag[ iModel ] = resolution[ iModel ];
        }
        save( resMDiag, "resMDiag.vec" );
        save( resM,     "resM" );
    }

    return EXIT_SUCCESS;
}
//    RVector rhoa2( RVector & rho, RVector & thk ) { //after WAIT.m by A. Junge
//        int nperiods = periods_.size();
//        static double mu0 = PI * 4e-7;
//        RVector rhoa( nperiods), phi( nperiods );
//        RVector freq( 2 * PI / periods_ );
//        Complex i_unit( 0.0 , 1.0);
//        CVector k1( nperiods), k2( nperiods);
//        CVector g( nperiods, Complex( 1.0, 0.0 ) );
//        for ( int k = nlay_ - 2 ; k >= 0 ; k-- ) {
//            k1 = sqrt( toComplex( freq * mu0 / rho[ k ] ) * i_unit );
//            save( k1, "k1.vec" );
////            k2 = sqrt( toComplex( freq * mu0 / rho[ k+1 ] ) * i_unit );
//            std::cout << k << " " << rho[ k + 1 ] << std::endl;
//            k2 = toComplex( freq * mu0 / rho[ k+1 ] );
//            save( k2, "k2.vec" );
//            g = ( g * k2 + k1 * tanh( k1 * thk[ k ] ) ) / ( k1 + g * k2 * tanh( k1 * thk[ k ] ) );
//            save( g, "g.vec" );
//        }
//        CVector z( toComplex( 0.0, freq ) / ( k1 * g ) );
//        save( z, "z.vec" );
//        rhoa = ( real(z) * real(z) + imag(z) * imag(z) ) / freq * mu0;
//        phi = angle( z );
////        phi = atan( imag( z ) / real( z ) );  //angle(z);
//        return rhoa;
//    }


