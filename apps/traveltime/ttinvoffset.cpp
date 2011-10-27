/***************************************************************************
 *   Copyright (C) 2008-2010 by the resistivity.net development team       *
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
#include <ttdijkstramodelling.h>
#include <datacontainer.h>
#include <mesh.h>
#include <meshgenerators.h>
#include <inversion.h>
#include <modellingbase.h>
#include <complex.h>
#include <string>
#include <vector.h>
#include <blockmatrix.h>

using namespace GIMLI;

#define vcout if ( verbose ) std::cout
#define dcout if ( debug ) std::cout
#define DEBUG if ( debug )
#define TTMod TravelTimeDijkstraModelling
#define NEWREGION 33333

/*! New Class derived from standard travel time modelling */
class TTOffsetModelling: public TravelTimeDijkstraModelling{
public:
    TTOffsetModelling( Mesh & mesh, DataContainer & dataContainer, bool verbose) :
        TravelTimeDijkstraModelling( mesh, dataContainer, verbose ) {
        //! find occuring shots, and maps them to indices starting from zero
        shots_ = unique( sort( dataContainer.get("s") ) );
        std::cout << "found " << shots_.size() << " shots." << std::endl;
        for ( size_t i = 0 ; i < shots_.size() ; i++ ) {
            shotMap_.insert( std::pair< int, int >( shots_[ i ], i ) );
        }
        //! create new region containing offsets and extend forward operator with it
        offsetMesh_ = createMesh1D( shots_.size() );
        for ( size_t i = 0 ; i < offsetMesh_.cellCount() ; i++ ) offsetMesh_.cell( i ).setMarker( NEWREGION );
        regionManager().createRegion( NEWREGION, offsetMesh_ );
        regionManager().recountParaMarker_();
    }

    virtual ~TTOffsetModelling() { }

    RVector createDefaultStartModel() { 
        return cat( TravelTimeDijkstraModelling::createDefaultStartModel(), RVector( shots_.size() ) );
    }
    
    RVector response( const RVector & model ) {
        RVector slowness( model, 0, model.size() - shots_.size() );
        RVector offsets( model, model.size() - shots_.size(), model.size() ); 
        RVector resp = TravelTimeDijkstraModelling::response( slowness ); //! normal response
        RVector shotpos = dataContainer_->get( "s" );
        for ( size_t i = 0; i < resp.size() ; i++ ){
            resp[ i ] += offsets[ shotMap_[ shotpos[ i ] ] ];
        }
        return resp;
    }

    void createJacobian( H2SparseMapMatrix & jacobian, const RVector & model ) { 
        RVector slowness( model, 0, model.size() - shots_.size() );
        RVector offsets( model, model.size() - shots_.size(), model.size() );         
        TravelTimeDijkstraModelling::createJacobian( jacobian.H1(), slowness );
        jacobian.H2().setRows( dataContainer_->size() );
        jacobian.H2().setCols( offsets.size() );
        //** set 1 entries for the used shot
        RVector shotpos = dataContainer_->get( "s" ); // shot=C1/A
        for ( size_t i = 0; i < dataContainer_->size(); i++ ) {
            jacobian.H2().setVal( i, shotMap_[ shotpos[ i ] ], 1.0 ); 
        }
    }
    
    size_t nShots(){ return shots_.size(); }
protected:
    RVector shots_;
    std::map< int, int > shotMap_;
    Mesh offsetMesh_;
};

int main( int argc, char *argv [] ) {
    bool verbose = true, createGradientModel = false;
    double errTime = 0.001;
    std::string dataFileName( NOT_DEFINED ), paraMeshFilename( "mesh.bms" );

    OptionMap oMap;
    oMap.setDescription("Description. TTInvOffset - travel time inversion with shot offsets\n");
    oMap.addLastArg( dataFileName, "Data file" );
    oMap.add( verbose,            "v" , "verbose"      , "Verbose output." );
    oMap.add( paraMeshFilename,   "p:", "paraMeshFile" , "Parameter mesh file." );
    oMap.add( createGradientModel,"G" , "gradientModel", "Gradient starting model." );
    oMap.parse( argc, argv );

    //!** load data file
    DataContainer dataIn( dataFileName );
    if ( verbose ) dataIn.showInfos();
    
    //!** load mesh for slownesses
    Mesh paraMesh;
    paraMesh.load( paraMeshFilename );
    
    /*! Error model combined of rhoa error and phi error */
    TTOffsetModelling f( paraMesh, dataIn, false );
    RVector appSlowness ( f.getApparentSlowness() );

    RTransLog transSlo;
    Region * mainRegion = f.regionManager().regions()->begin()->second;
    mainRegion->setConstraintType( 1 ); //! smoothness
    mainRegion->setStartValue( median( appSlowness ) );
    mainRegion->setTransModel( transSlo );
    f.region( NEWREGION )->setConstraintType( 0 ); //! minimum length (no smoothing)
    f.region( NEWREGION )->setStartValue( 0.0 );
    RVector model;

    /*! Set up inversion with full matrix, data and forward operator */
    Inversion< double, H2SparseMapMatrix > inv( dataIn.get( "t" ), f, verbose );
    inv.setAbsoluteError( errTime );

    /*! actual computation: run the inversion */
    model = inv.run();
    RVector slowness( model, 0, model.size() - f.nShots() );
    RVector velocity( 1 / slowness );
    save( velocity, "velocity.vec" );
    RVector offsets( model, model.size() - f.nShots(), model.size() ); 
    save( offsets, "offsets.vec" );
    
    vcout << "min/max( velocity ) = " << 1 / max( slowness ) << "/" << 1 / min( slowness ) << std::endl;
    vcout << "offsets: " << offsets << std::endl;

    return EXIT_SUCCESS;
}
