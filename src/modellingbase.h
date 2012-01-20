/***************************************************************************
 *   Copyright (C) 2005-2012 by the resistivity.net development team       *
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

#ifndef _GIMLI_MODELLINGBASE__H
#define _GIMLI_MODELLINGBASE__H

#include "gimli.h"
#include "matrix.h"
#include "blockmatrix.h"

namespace GIMLI{
    
//class H2SparseMapMatrix;    

/*! General modelling interface class.*/
class DLLEXPORT ModellingBase{
public:

    ModellingBase( bool verbose = false ) ;

    ModellingBase( DataContainer & dataContainer, bool verbose = false );

    ModellingBase( const Mesh & mesh, bool verbose = false );

    ModellingBase( const Mesh & mesh, DataContainer & dataContainer, bool verbose = false );

    virtual ~ModellingBase( );

    virtual RVector response( const RVector & model ) = 0;

    inline RVector operator() ( const RVector & model ){ return response( model ); }

    virtual RVector createDefaultStartModel( ) { return RVector( 0 ); }

    /*! Change the associated data container */
    void setData( DataContainer & data );
    
    /*! */
    virtual RVector startModel( );

    virtual void setStartModel( const RVector & startModel ){ startModel_ = startModel; }

/*! Set new mesh to the forward operator, optionaly hold region parameter for the new mesh (i.e. for rollalong)*/
    void setMesh( const Mesh & mesh, bool holdRegionInfos = false );

    inline Mesh * mesh() { return mesh_; }

    /*! Set external jacobian */
    virtual void setJacobian( MatrixBase * J );
    
    /*! Create and fill the Jacobian matrix with a given model vector. */
    virtual void createJacobian( const RVector & model );
    
    /*! Here you should initialize your Jacobian matrix. Default is RMatrix()*/
    virtual void initJacobian( );
    
    /*! Return the pointer to the Jacobian matrix associated with this forward operator. */
    virtual MatrixBase * jacobian( ) { return jacobian_; }
    
    /*! Return the pointer to the Jacobian matrix associated with this forward operator. */
    virtual MatrixBase * jacobian( ) const { return jacobian_; }
    
    /*! Return the Jacobian Matrix (read only) associated with this forward operator. 
     *  Throws an exception if the jacobian is not initialized. Cannot yet be overloaded py pyplusplus (Warning 1049). */
    virtual RMatrix & jacobianRef( ) const { 
        if ( ! jacobian_ ) {
            throwError( 1, WHERE_AM_I + " Jacobian matrix is not initialized." );
        }
        return *dynamic_cast< RMatrix * >( jacobian_ ); 
    }
    
    /*! Clear Jacobian matrix. */
    virtual void clearJacobian( ){ jacobian_->clear(); }
    
    const RMatrix & solution() const { return solutions_; }

    void createRefinedForwardMesh( bool refine = true, bool pRefine = false );
    
    /*! Set refined mesh for forward calculation. */
    void setRefinedMesh( const Mesh & mesh );

    void mapModel( const RVector & model, double background = 0 );

    const RegionManager & regionManager() const;

    RegionManager & regionManager();

    void setVerbose( bool verbose ) { verbose_ = verbose; }

    Region * region( int marker );

    RVector createStartVector( );

    void initRegionManager( );
     
    /*! Set the maximum number of allowed threads for MT calculation. Have to be greater than 0 */
    void setThreadCount( uint nThreads ) { nThreads_ = max( 1, (int)nThreads ); }
    
    /*! Return the maximum number of allowed threads for MT calculation */
    uint threadCount( ) const { return nThreads_; }

protected:
    virtual void init_();
    
    virtual void deleteMeshDependency_(){}

    virtual void updateMeshDependency_(){}

    virtual void updateDataDependency_(){}

    void setMesh_( const Mesh & mesh );
    
    Mesh                    * mesh_;

    DataContainer           * dataContainer_;

    MatrixBase              * jacobian_;
    
    bool                    ownJacobian_;
    
    RMatrix                     solutions_;

    RVector                 startModel_;
    
    bool                       verbose_;
    
    bool                       regionManagerInUse_;
    
    uint nThreads_;
    
private:
    RegionManager            * regionManager_;

};

//! Simple linear modelling class
/*! Linear modelling class.
    Calculate response from derivative matrix (Jacobian,Sensitivity) * model
*/
class DLLEXPORT LinearModelling : public ModellingBase {
public:
    LinearModelling( Mesh & mesh, const MatrixBase & A, bool verbose = false )
        : ModellingBase( mesh, verbose ), A_( &A ) { }

    LinearModelling( Mesh & mesh, const MatrixBase * A, bool verbose = false )
        : ModellingBase( mesh, verbose ), A_( A ) { }

    LinearModelling( const MatrixBase & A, bool verbose = false );

    virtual ~LinearModelling() { }

    /*! Calculate */
    RVector response( const RVector & model );

    RVector createDefaultStartModel( );

protected:

    const MatrixBase * A_;
};



} // namespace GIMLI

#endif // MODELLINGBASE__H
