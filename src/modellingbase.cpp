/***************************************************************************
 *   Copyright (C) 2005-2011 by the resistivity.net development team       *
 *   Carsten Rücker carsten@resistivity.net                                *
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

#include "modellingbase.h"

#include "datacontainer.h"
#include "mesh.h"
#include "regionManager.h"
#include "stopwatch.h"
#include "vector.h"
#include "vectortemplates.h"

namespace GIMLI{

ModellingBase::ModellingBase( bool verbose )
    : dataContainer_( NULL ), verbose_( verbose ){
    init_();
}

ModellingBase::ModellingBase( const DataContainer & data, bool verbose )
    : dataContainer_( NULL ), verbose_( verbose ){
    init_();
    setData( data );
}

ModellingBase::ModellingBase( Mesh & mesh, bool verbose )
    : dataContainer_( NULL ), verbose_( verbose ){
    init_();
    setMesh( mesh );
}

ModellingBase::ModellingBase( Mesh & mesh, const DataContainer & data, bool verbose )
    : dataContainer_( NULL ), verbose_( verbose ){
    init_();
    setData( data );
    setMesh( mesh );
}

ModellingBase::~ModellingBase( ) {
    delete regionManager_;
    if ( mesh_ ) delete mesh_;
    if ( dataContainer_ ) delete dataContainer_;
}


void ModellingBase::init_() {
    regionManager_      = new RegionManager( verbose_ );
    mesh_               = NULL;
    regionManagerInUse_  = false;
}

void ModellingBase::setData( const DataContainer & data ){
    if ( dataContainer_ ) {
        *dataContainer_ = data;
    } else {
        dataContainer_ = new DataContainer( data );
    }
    updateDataDependency_();
}

RVector ModellingBase::startModel( ) {
    if ( startModel_.size() > 0 ) return startModel_;

    setStartModel( regionManager_->createStartVector( ) );

    if ( startModel_.size() == 0 ){
/*            std::cerr << WHERE_AM_I << "Warning! no mesh defined. "
                "You should either define a mesh or override this function." << std::endl;*/

        setStartModel( createDefaultStartModel() );
    }
    return startModel_;
}

void ModellingBase::createRefinedForwardMesh( bool refine, bool pRefine ){
    if ( !mesh_ ) {
        mesh_ = new Mesh();
    } else {
        mesh_->clear();
    }
    
    if ( refine ){
        if ( pRefine ){
            mesh_->createP2Mesh( regionManager_->mesh() );
        } else {
            mesh_->createH2Mesh( regionManager_->mesh() );
        }
    } else {
        setMesh_( regionManager_->mesh() );
    }
    updateMeshDependency_();
}
 
void ModellingBase::setRefinedMesh( const Mesh & mesh ){
    if ( verbose_ ) {
        std::cout << "set extenal secondary mesh:" << std::endl;
    }
    setMesh_( mesh );
    if ( verbose_ ) {
        std::cout << "nModel = " << regionManager_->parameterCount() << std::endl;
        std::vector < int > m( unique( sort( mesh_->cellMarker() ) ) );
        std::cout << "secMesh marker = [ " << m[0] <<", " << m[1] << ", " << m[2]  
         << ", ... ,  " << m.back() << " ]" << std::endl;
    }
    updateMeshDependency_();
}
 
void ModellingBase::setMesh( const Mesh & mesh, bool holdRegionInfos ) {
    deleteMeshDependency_();

    Stopwatch swatch( true );
    if ( regionManagerInUse_ ){
        regionManager_->setMesh( mesh, holdRegionInfos );
        if ( verbose_ ) std::cout << "ModellingBase::setMesh() switch to regionmanager mesh" << std::endl;
        this->setMesh_( regionManager_->mesh() );
    } else {
        if ( verbose_ ) std::cout << "ModellingBase::setMesh() copying new mesh ... ";
        this->setMesh_( mesh );
        if ( verbose_ ) std::cout << swatch.duration( true ) << " s" << std::endl;
    }

    if ( verbose_ ) std::cout << "FOP updating mesh dependencies ... ";
    startModel_.clear();
    updateMeshDependency_();
    if ( verbose_ ) std::cout << swatch.duration( true ) << " s" << std::endl;
}

void ModellingBase::setMesh_( const Mesh & mesh ){
    if ( !mesh_ ) mesh_ = new Mesh();
    (*mesh_) = mesh;
}

/*! create Jacobian matrix using brute force approach */
void ModellingBase::createJacobian( RMatrix & jacobian, const RVector & model ){
    if ( verbose_ ) std::cout << "Create Jacobian matrix (brute force) ...";

    Stopwatch swatch( true );
    RVector resp( response( model ) );
    if ( jacobian.rows() != resp.size() ){ jacobian.resize( resp.size(), model.size() ); }

    double fak = 1.05;
    for ( size_t i = 0; i < model.size(); i++ ) {
        RVector modelChange( model );
        modelChange[ i ] *= fak;
        RVector respChange( response( modelChange ) );

        for ( size_t j = 0; j < resp.size(); j++ ){
            jacobian[ j ][ i ] = ( respChange[ j ] - resp[ j ] ) / ( modelChange[ i ] - model[ i ] );
        }
    }

    swatch.stop();
    if ( verbose_ ) std::cout << " ... " << swatch.duration() << " s." << std::endl;
}

void ModellingBase::createJacobian( DSparseMapMatrix & jacobian, const RVector & model ){
    CERR_TO_IMPL
}

void ModellingBase::mapModel( const RVector & model, double background ){
    int marker = -1;
    std::vector< Cell * > emptyList;
    mesh_->createNeighbourInfos();


    for ( uint i = 0, imax = mesh_->cellCount(); i < imax; i ++ ){
        marker = mesh_->cell( i ).marker();
        if ( marker >= 0 ) {
            if ( (size_t)marker >= model.size() ){
                mesh_->exportVTK( "mapModelfail" );
                throwLengthError( 1, WHERE_AM_I + " marker greater= then model.size()" + toStr( marker )
                       + " != " + toStr( model.size() ) );
            }
            if ( model[ marker ] < TOLERANCE ){
                emptyList.push_back( &mesh_->cell( i ) );
            }
            mesh_->cell( i ).setAttribute( model[ marker ] );

        } else {
            mesh_->cell( i ).setAttribute( 0.0 );
            emptyList.push_back( &mesh_->cell( i ) );
        }
    }

    if ( emptyList.size() == mesh_->cellCount() ){
        throwLengthError( 1, WHERE_AM_I + " to many empty cells" + toStr( emptyList.size() )
                       + " == " + toStr( mesh_->cellCount() ) );
    }

    mesh_->fillEmptyCells( emptyList, background );
}

void ModellingBase::initRegionManager() {
    if ( !regionManagerInUse_ ){
        if ( mesh_ ){
            regionManager_->setMesh( *mesh_ );
            this->setMesh_( regionManager_->mesh() );
        }
        regionManagerInUse_   = true;
    }
}

const RegionManager & ModellingBase::regionManager() const {
    return *regionManager_;
}

RegionManager & ModellingBase::regionManager(){
    initRegionManager( );
    return *regionManager_;
}

Region * ModellingBase::region( int marker ){
    return regionManager().region( marker );
}

RVector ModellingBase::createStartVector( ) {
    return regionManager().createStartVector();
}


LinearModelling::LinearModelling( const RMatrix * A, bool verbose )
    : ModellingBase( verbose ), A_( A ) {
    this->regionManager().setParameterCount( A_->cols() );
}

RVector LinearModelling::response( const RVector & model ) {
    if ( A_->cols() != model.size() ){
        throwLengthError( 1, WHERE_AM_I + " Jacobian col size != model.size()"
                                + toStr( A_->cols() ) + " != " + toStr( model.size() ) );
    }
    return *A_ * model;
}

RVector LinearModelling::createDefaultStartModel( ) {
    return RVector( A_->cols(), 1.0 );
}

} // namespace GIMLI

