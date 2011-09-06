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

#ifndef _GIMLI_MODELLINGBASE__H
#define _GIMLI_MODELLINGBASE__H

#include "gimli.h"
#include "matrix.h"

namespace GIMLI{
    
//class H2SparseMapMatrix;    

/*! General modelling interface class.*/
class DLLEXPORT ModellingBase{
public:

    ModellingBase( bool verbose = false ) ;

    ModellingBase( const DataContainer & dataContainer, bool verbose = false );

    ModellingBase( Mesh & mesh, bool verbose = false );

    ModellingBase( Mesh & mesh, const DataContainer & dataContainer, bool verbose = false );

    virtual ~ModellingBase( );

    virtual RVector response( const RVector & model ) = 0;

    inline RVector operator() ( const RVector & model ){ return response( model ); }

    virtual RVector createDefaultStartModel( ) { return RVector( 0 ); }

    /*! Change the associated data container */
    void setData( const DataContainer & data );
    
    /*! */
    virtual RVector startModel( );

    virtual void setStartModel( const RVector & startModel ){ startModel_ = startModel; }

/*! Set new mesh to the forward operator, optionaly hold region parameter for the new mesh (i.e. for rollalong)*/
    void setMesh( const Mesh & mesh, bool holdRegionInfos = false );

    inline Mesh * mesh() { return mesh_; }

#ifndef PYGIMLI
    // temporary exclusion hack for python bindings, think about abstract base class for matrices
    template < class Matrix > void createJacobian( Matrix & jacobian, const RVector & model );
#endif
    virtual void createJacobian( RMatrix & jacobian, const RVector & model );

    virtual void createJacobian( DSparseMapMatrix & jacobian, const RVector & model );

//     virtual void createJacobian( H2SparseMapMatrix & jacobian, const RVector & model ){}

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
    virtual void deleteMeshDependency_(){}

    virtual void updateMeshDependency_(){}

    virtual void updateDataDependency_(){}

    void setMesh_( const Mesh & mesh );

    virtual void init_();

    Mesh                     * mesh_;

    DataContainer            * dataContainer_;

    bool                       verbose_;

    //region std::vector < int >        cellMapIndex_;
    RMatrix                    solutions_;

    RegionManager            * regionManager_;
    bool                       regionManagerInUse_;

    RVector                    startModel_;
    
    uint nThreads_;
};

//! Simple linear modelling class
/*! Linear modelling class.
    Calculate response from derivative matrix (Jacobian,Sensitivity) * model
*/
class DLLEXPORT LinearModelling : public ModellingBase {
public:
    LinearModelling( Mesh & mesh, const RMatrix & A, bool verbose = false )
        : ModellingBase( mesh, verbose ), A_( &A ) { }

    LinearModelling( Mesh & mesh, const RMatrix * A, bool verbose = false )
        : ModellingBase( mesh, verbose ), A_( A ) { }

    LinearModelling( const RMatrix * A, bool verbose = false );

    virtual ~LinearModelling() { }

    /*! Calculate */
    RVector response( const RVector & model );

    RVector createDefaultStartModel( );

    //void createJacobian( RMatrix & jacobian, const RVector & model ) { return *A_; }
protected:

    const RMatrix * A_;
};



} // namespace GIMLI

#endif // MODELLINGBASE__H
