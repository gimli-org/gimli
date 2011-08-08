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

#ifndef _GIMLI_REGIONMANAGER__H
#define _GIMLI_REGIONMANAGER__H

#include "gimli.h"
#include "vector.h"
#include "trans.h"

#include <set>
#include <list>

namespace GIMLI{

class PosFunctor{
public:
    virtual ~PosFunctor(){}
    virtual const double operator()( const RVector3 & pos ) const = 0;
};

class DLLEXPORT Region{
public:
    Region( int marker, RegionManager * parent, bool single = false );

    Region( int marker, const Mesh & mesh, RegionManager * parent );

    Region( const Region & region );

    Region & operator = ( const Region & region );

    ~Region( );

    /*! Set the region marker id */
    inline void setMarker( int marker ) { marker_ = marker; }

    /*! Return the region marker id. */
    inline int marker() const { return marker_; }

    /*! Returns read_only acccess to the boundaries of this region */
    const std::vector < Boundary * > & boundaries() const { return bounds_; }

    /*! Set new parameter cells, i.e. update the related mesh and all sizes. */
    void resize( const Mesh & mesh );

    /*! Set new parameter cells, i.e. update the related mesh and all sizes. Only for single region.*/
    void resize( const std::vector < Cell * > & cells );

    /*! Mark this region to be a background region, need RegionManger::recount */
    inline void markBackground( bool background ){ isBackground_ = background; }

    /*! Set this region to be a background region. Forces the regionmanager to recount all regions. */
    void setBackground( bool background = true );

    /*! Return true if this region is a parameter region */
    inline bool isBackground() const { return isBackground_ ; }

    /*! Mark this region to be a single parameter region, need RegionManger::recount */
    inline void markSingle( bool issingle ){ isSingle_ = issingle; }

    /*! Set this region to be a single parameter region. Forces the regionmanager to recount all regions. */
    void setSingle( bool background = true );

    /*! Return true if this region is a single parameter region */
    inline bool isSingle() const { return isSingle_ ; }

    /*! Return amount of parameter for this region, 1 on single region */
    inline uint parameterCount() const { return parameterCount_; }

    /*! Returns the first parameter id by means of a global count, defined in countParameter called from RegionManager */
    inline uint startParameter() const { return startParameter_; }
    
    /*! Returns the last parameter id by means of a global count, defined in countParameter called from RegionManager */
    inline uint endParameter() const { return endParameter_; }


    /*! Set the constrainst-type for this region available( 0, 1 ) */
    inline void setConstraintType( uint type );
    
    /*! Return constrainst-type for this region. */
    inline uint constraintType() const { return constraintType_; }

    /*! Returns amount if constraints defined for this region,
        For single region return 1,\n
        for constraintstype == 0 return amount of parameter cells,\n
        else return amount of inner boundaries */
    uint constraintCount() const;
    
    /*! Fill given constraints matrix with local constraints. 
        For single region (startConstraintsID, startParameter_ ) = 1, \n
        for constraintstype == 0 fill  (startConstraintsID + i, startParameter_ + i ) = 1, i=1..constraintCount()\n
        else fill (startConstraintsID + i, Boundary_i_leftNeightbourParameterID) = 1, 
                  (startConstraintsID + i, Boundary_i_rightNeightbourParameterID) = -1, i = 1..nBoundaries.*/
    void fillConstraints( DSparseMapMatrix & C, uint startConstraintsID );

    /*! Set region wide constant constraints weight, (default = 1 ). If this method is called background is forced to false. */
    void setConstraintsWeight( double bc );
    
    /*! Set region wide variable constraints weight from RVector cw. If this method is called background is forced to false. */
    void setConstraintsWeight( const RVector & cw );
    
    /*! Return Read-only vector of given constraints weights as RVector (default RVector( constraintCount, 1.0) */
    inline const RVector & constraintsWeight(  ) const { return constraintsWeight_; }
    
    /*! Fill global constraints weight vector started at constraintStart. */ 
    void fillConstraintsWeight( RVector & vec, uint constraintStart  );
    
    /*! Set Region-Wide horizontal(z) weighting parameter for anisotropic smoothing \n
        1 - isotrope, 0 -- no vertical smoothing
    */
    inline void setZWeight( double zw ){ 
        zWeight_ = zw; 
        this->fillConstraintsWeightWithFlatWeight(); 
    }
    /*! Return Region-Wide horizontal(z)-weighting parameter*/
    inline double zWeight( ) const { return zWeight_; }
    
    /*! Possible DEPRECATED ?? */
    inline void setZPower( double zp ){ 
        zPower_ = zp; zWeight_ = 0.01; 
        fillConstraintsWeightWithFlatWeight( );
    }
    inline double zPower( ) const { return zPower_; }
    
    /*! Helper method that convert cweight parameter into individual constraintsWeights depending on the associated boundary norm. At the moment only zWeight is considered. */
    void fillConstraintsWeightWithFlatWeight( );
    
//DEPRECATED ??
//      inline RVector * constraintsWeight(  ) { return & boundaryControl_; }

    void fillBoundaryNorm( std::vector< RVector3 > & vnorm, uint boundCount );
 
    void fillBoundarySize( RVector & vec, uint boundStart );
    
    void fillStartVector( RVector & vec );

    void fillModelControl( RVector & vec );

    const std::vector < Cell * > & cells() const { return cells_; }

    std::vector < Cell * > & cells() { return cells_; }

    void countParameter( uint start );

    void setStartVector( const RVector & start );
    void setStartValue( double start );

    void setModelControl( double mc );
    void setModelControl( const RVector & mc );
    void setModelControl( PosFunctor * mcF );

    inline const RVector & modelControl(  ) const { return modelControl_; }
    inline RVector * modelControl(  ) { return & modelControl_; }

    void setTransModel( Trans< RVector > & tM );

    /*! Return a cumulative transform function based on transform functions for each region.
        Return NULL if no region is defined. */
    Trans< RVector > * transModel() { return tM_; }

    void setLowerBound( double lb );
    void setUpperBound( double ub );

    void setModelTransStr_( const std::string & val );

    void setModelControlStr_(   const std::string & val ){ setModelControl( toDouble( val ) ); }
    void setStartValueStr_(     const std::string & val ){ setStartValue( toDouble( val ) ); }
    void setZPowerStr_(         const std::string & val ){ setZPower( toDouble( val ) ); }
    void setZWeightStr_(        const std::string & val ){ setZWeight( toDouble( val ) ); }
    void setConstraintTypeStr_( const std::string & val ){ setConstraintType( toInt( val ) ); }
    void setLowerBoundStr_(     const std::string & val ){ setLowerBound( toDouble( val ) ); }
    void setUpperBoundStr_(     const std::string & val ){ setUpperBound( toDouble( val ) ); }
    void setSingleStr_(         const std::string & val ){ markSingle( (bool)toInt( val ) ); }
    void setBackgroundStr_(     const std::string & val ){ markBackground( (bool)toInt( val ) ); }

protected:
    void init_();
    void copy_( const Region & region );

    int marker_;
    RegionManager * parent_;

    std::vector < Cell * > cells_;
    mutable std::vector < Boundary * > bounds_;
    
    bool isBackground_;
    bool isSingle_;

    uint parameterCount_;
    uint startParameter_;
    uint endParameter_;

    uint constraintType_;

    RVector startVector_;
    RVector modelControl_;
    RVector constraintsWeight_;

    double zPower_;
    double zWeight_;
    
    double mcDefault_;
    double startDefault_;

//    double minZWeight_;

    double lowerBound_;
    double upperBound_;

    Trans< RVector >        * tM_;
    bool ownsTrans_;        // smart ptr. would be nice here.
    std::string transString_;
};

class DLLEXPORT RegionManager{
public:
    RegionManager( bool verbose = true );

    RegionManager & operator = ( const RegionManager & rm );

    ~RegionManager();

    void clear();

    void setMesh( const Mesh & mesh, bool holdRegionInfos = false );

    const Mesh & mesh() const { return *mesh_; }

    //inline Mesh * mesh(  ) { return mesh_; }

    Region * createRegion( int marker, const Mesh & mesh );

    Region * createSingleRegion( int marker, const std::vector < Cell * > & cells );

    const std::map < int, Region * > & regions() const { return regionMap_; }

    std::map < int, Region * >  * regions() { return &regionMap_; }

    uint regionCount( ) const { return regionMap_.size(); }

    /*! Returns a ptr to the region with the given marker. If no region exist an exception is thrown. */
    Region * region( int marker );

    inline bool regionExists( int marker ) { return ( regionMap_.count( marker ) > 0 ); }

    void setInterRegionConstraint( int a, int b, double c );

    void loadMap( const std::string & fname );

    void saveMap( const std::string & fname );

    /*! Set the amount of parameter, will be override if regions are defined */
    inline void setParameterCount( uint count ) { parameterCount_ = count; }

    /*! Return global amount of parameter */
    uint parameterCount() const;

    /*! Return global amount of constraints */
    uint constraintCount() const;

    uint interRegionConstraintsCount() const;

    void fillStartVector( RVector & vec );

    RVector createStartVector( );

    /*! Create and fill global model-weight vector */
    RVector createModelControl( );

    /*! Fill global model-weight vector */
    void fillModelControl( RVector & vec );
    
    /*! Create and fill global constraints-weight vector */
    RVector createConstraintsWeight( );

    /*! Fill global constraints-weight vector */
    void fillConstraintsWeight( RVector & vec );
    
    /*! Fill global constraints-matrix
        no regions: fill with 0th-order constraints */
    void fillConstraints( DSparseMapMatrix & C );

    /*! Syntactic sugar: set zweight/constraintType to all regions. */
    void setZWeight( double z );

    void setConstraintType( uint type );

    void fillBoundarySize( RVector & vec );

    inline const Mesh & paraDomain() const { return *paraDomain_; }

    inline Mesh & paraDomain() { return *paraDomain_; }

    std::vector < RVector3 > boundaryNorm() const;

//     RVector createFlatWeight( double zPower, double zWeight ) const;

    void createParaDomain_();
    void recountParaMarker_();

    CumulativeTrans< RVector > * transModel();

    void setLocalTransFlag( bool flag ) { haveLocalTrans_ = flag; }

    bool haveLocalTrans( ) const { return haveLocalTrans_; }

protected:
    void copy_( const RegionManager & rm );

    /*! Fill \ref interRegionInterfaceMap_ */
    void findInterRegionInterfaces_();

    std::vector < int > allRegionMarker_( bool exludeBoundary = false ){
        std::vector < int > tmp;
        for ( std::map < int, Region * > ::const_iterator it  = regionMap_.begin();
                                                          it != regionMap_.end(); it++ ){
            if ( !it->second->isBackground() ) tmp.push_back( it->first );
        }
        return tmp;
    }


    bool verbose_;

    uint parameterCount_;

    Mesh * mesh_;
    Mesh * paraDomain_;

    std::map < int, Region * > regionMap_;
    std::map< std::pair< int, int >, std::list < Boundary * > > interRegionInterfaceMap_;
    std::map< std::pair< int, int >, double > interRegionConstraints_;
    std::map< int, double > interfaceConstraint_;

    double interRegionConstraintsZWeight_;

    CumulativeTrans< RVector >  localTrans_;
    bool haveLocalTrans_;
};


} // namespace GIMLI

#endif // _GIMLI_REGIONMANAGER__H
