/******************************************************************************
 *   Copyright (C) 2008-2017 by the GIMLi development team                    *
 *   Carsten Rücker carsten@resistivity.net                                   *
 *   Thomas Günther thomas@resistivity.net                                    *
 *                                                                            *
 *   Licensed under the Apache License, Version 2.0 (the "License");          *
 *   you may not use this file except in compliance with the License.         *
 *   You may obtain a copy of the License at                                  *
 *                                                                            *
 *       http://www.apache.org/licenses/LICENSE-2.0                           *
 *                                                                            *
 *   Unless required by applicable law or agreed to in writing, software      *
 *   distributed under the License is distributed on an "AS IS" BASIS,        *
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 *   See the License for the specific language governing permissions and      *
 *   limitations under the License.                                           *
 *                                                                            *
 ******************************************************************************/

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
    virtual const double operator()(const RVector3 & pos) const = 0;
};

class DLLEXPORT Region{
public:
    Region(SIndex marker, RegionManager * parent, bool single=false);

    Region(SIndex marker, const Mesh & mesh, RegionManager * parent);

    Region(const Region & region);

    Region & operator = (const Region & region);

    ~Region();

    /*! Set the region marker id */
    inline void setMarker(SIndex marker) { marker_ = marker; }

    /*! Return the region marker id. */
    inline SIndex marker() const { return marker_; }

    /*! Returns read_only acccess to the boundaries of this region */
    const std::vector < Boundary * > & boundaries() const { return bounds_; }

    /*! Set new parameter cells, i.e. update the related mesh and all sizes. */
    void resize(const Mesh & mesh);

    /*! Set new parameter cells, i.e. update the related mesh and all sizes. Only for single region.*/
    void resize(const std::vector < Cell * > & cells);

    /*! Mark this region to be a background region, need RegionManger::recount */
    inline void markBackground(bool background){ isBackground_ = background; }

    /*! Set this region to be a background region. Forces the regionmanager to recount all regions. */
    void setBackground(bool background = true);

    /*! Return true if this region is a parameter region */
    inline bool isBackground() const { return isBackground_ ; }

    /*! Mark this region to be a single parameter region, need RegionManger::recount */
    inline void markSingle(bool issingle){ isSingle_ = issingle; }

    /*! Set this region to be a single parameter region. Forces the regionmanager to recount all regions. */
    void setSingle(bool background = true);

    /*! Return true if this region is a single parameter region */
    inline bool isSingle() const { return isSingle_ ; }

    /*! Return number of parameters for this region, 1 on single region */
    inline Index parameterCount() const { return parameterCount_; }

    /*! Return the first parameter id by means of a global count, defined in countParameter called from RegionManager */
    inline Index startParameter() const { return startParameter_; }

    /*! Return the last parameter id by means of a global count, defined in countParameter called from RegionManager */
    inline Index endParameter() const { return endParameter_; }

    /*! Set the constraint type for this region available(0, 1) */
    void setConstraintType(Index type);

    /*! Return constraint type for this region. */
    inline Index constraintType() const { return constraintType_; }

    /*! Returns number of constraints defined for this region,
        For single region return 1,\n
        for constraintstype == 0 return amount of parameter cells,\n
        else return amount of inner boundaries */
    Index constraintCount() const;

    /*! Fill given constraints matrix with local constraints.
        For single region (startConstraintsID, startParameter_) = 1, \n
        for constraintstype == 0 fill  (startConstraintsID + i, startParameter_ + i) = 1, i=1..constraintCount()\n
        else fill (startConstraintsID + i, Boundary_i_leftNeightbourParameterID) = 1,
                  (startConstraintsID + i, Boundary_i_rightNeightbourParameterID) = -1, i = 1..nBoundaries.*/
    void fillConstraints(RSparseMapMatrix & C, Index startConstraintsID);

    /*! Set region wide constant constraints weight, (default = 1). If this method is called background is forced to false. */
    void setConstraintsWeight(double bc);

    /*! Set region wide variable constraints weight from RVector cw. If this method is called background is forced to false. */
    void setConstraintsWeight(const RVector & cw);

    /*! Return Read-only vector of given constraints weights as RVector (default RVector(constraintCount, 1.0) */
    inline const RVector & constraintsWeight() const { return constraintsWeight_; }

    /*! Fill global constraints weight vector started at constraintStart. */
    void fillConstraintsWeight(RVector & vec, Index constraintStart);

    /*! Set Region-Wide horizontal(z) weighting parameter for anisotropic smoothing \n
        1 - isotrope, 0 -- no vertical smoothing
    */
    inline void setZWeight(double zw){
        zWeight_ = zw;
        this->fillConstraintsWeightWithFlatWeight();
    }
    /*! Return Region-Wide horizontal(z)-weighting parameter*/
    inline double zWeight() const { return zWeight_; }

    /*! Possible DEPRECATED ?? */
    inline void setZPower(double zp){
        __MS("DEPRECATED")
        zPower_ = zp; zWeight_ = 0.01;
        fillConstraintsWeightWithFlatWeight();
    }
    inline double zPower() const { return zPower_; }

    /*! Set fixed value for background regions that will not
     * part of any value prolongation.*/
    void setFixValue(double val);
    inline double fixValue() const { return fixValue_;}

    /*! Helper method that convert cWeight parameter into individual
     * constraintsWeights depending on the associated boundary norm.
     * At the moment only zWeight is considered. */
    void fillConstraintsWeightWithFlatWeight();

//DEPRECATED ??
//      inline RVector * constraintsWeight() { return & boundaryControl_; }

    void fillBoundaryNorm(std::vector< RVector3 > & vnorm, Index boundCount);

    void fillBoundarySize(RVector & vec, Index boundStart);

    void fillStartModel(RVector & vec);

//     void fillStartVector(RVector & vec);

    void fillModelControl(RVector & vec);

    const std::vector < Cell * > & cells() const { return cells_; }

    std::vector < Cell * > & cells() { return cells_; }

    void countParameter(Index start);

    void permuteParameterMarker(const IVector & p);


    /*! Set the values of start for the start model of this region. */
    void setStartModel(const RVector & start);
    /*! Set the value of start into the start model vector for this region. */
    void setStartModel(double start);

//     /*! DEPRECATED use setStartModel */
//     void setStartVector(const RVector & start);
//     /*! DEPRECATED use setStartModel */
    void setStartValue(double start){ DEPRECATED setStartModel(start);}

    void setModelControl(double val);
    void setModelControl(const RVector & vec);
    void setModelControl(PosFunctor * mcF);

    inline const RVector & modelControl() const { return modelControl_; }
    inline RVector * modelControl() { return & modelControl_; }

    void setModelTransformation(const Trans< RVector > & tM);
    void setTransModel(Trans< RVector > & tM);

    /*! Return a cumulative transform function based on transform functions for each region.
        Return NULL if no region is defined. */
    Trans< RVector > * transModel() { return tM_; }

    /*! set lower parameter bound for region */
    void setLowerBound(double lb);
    /*! set lower parameter bound for region */
    void setUpperBound(double ub);

    /*! set start and upper/lower bounds for region */
    void setParameters(double start, double lb, double ub, std::string transString="");

    /*! Create Constraints weight values from constrains type and weighting.*/
    void createConstraintsWeight_();

    void setModelTransStr_(const std::string & val);

    void setModelControlStr_(  const std::string & val){ setModelControl(toDouble(val)); }
    void setStartModelStr_(    const std::string & val){ setStartModel(toDouble(val)); }
    void setZPowerStr_(        const std::string & val){ setZPower(toDouble(val)); }
    void setZWeightStr_(       const std::string & val){ setZWeight(toDouble(val)); }
    void setFixValueStr_(      const std::string & val){ setFixValue(toDouble(val)); }
    void setConstraintTypeStr_(const std::string & val){ setConstraintType(toInt(val)); }
    void setLowerBoundStr_(    const std::string & val){ setLowerBound(toDouble(val)); }
    void setUpperBoundStr_(    const std::string & val){ setUpperBound(toDouble(val)); }
    void setSingleStr_(        const std::string & val){ markSingle(toInt(val) != 0); }
    void setBackgroundStr_(    const std::string & val){ markBackground(toInt(val) != 0); }

protected:
    void init_();
    void copy_(const Region & region);

    SIndex marker_;
    RegionManager * parent_;

    std::vector < Cell * > cells_;
    mutable std::vector < Boundary * > bounds_;

    bool isBackground_;
    bool isSingle_;
    bool isPermuted_;

    Index parameterCount_;
    Index startParameter_;
    Index endParameter_;

    Index constraintType_;

    RVector startModel_;
    RVector modelControl_;
    RVector constraintsWeight_;

    double zPower_;
    double zWeight_;
    double fixValue_;

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
    RegionManager(bool verbose=true);

    RegionManager & operator = (const RegionManager & rm);

    ~RegionManager();

    void setVerbose(bool verbose){ verbose_ = verbose; }
    bool verbose() const { return verbose_; }

    void clear();

    void setMesh(const Mesh & mesh, bool holdRegionInfos=false);

    const Mesh & mesh() const;

    inline Mesh * pMesh() { return mesh_; }

    /*!
     * Add a external region to the RegionManager.
     */
    Region * addRegion(SIndex marker, const Mesh & mesh);
#ifndef PYGIMLI_CAST
    const std::map < SIndex, Region * > & regions() const { return regionMap_; }

    std::map < SIndex, Region * >  * regions() { return &regionMap_; }
#endif
    Index regionCount() const { return regionMap_.size(); }

    /*!Return all region Indieces.*/
    IVector regionIdxs() const;

    /*! Returns a ptr to the region with the given marker. If no region exist an exception is thrown. */
    Region * region(SIndex marker);

    inline bool regionExists(SIndex marker) { return (regionMap_.count(marker) > 0); }

    void setInterRegionConstraint(SIndex a, SIndex b, double c);

    /*! load region parameters from region control file */
    void loadMap(const std::string & fname);

    /*! save region parameters to region control file */
    void saveMap(const std::string & fname);

    /*! Set the amount of parameter, will be override if regions are defined */
    inline void setParameterCount(Index count) { parameterCount_ = count; }

    /*! Return global amount of parameter */
    Index parameterCount() const;

    /*! Return global amount of constraints */
    Index constraintCount() const;

    Index interRegionConstraintsCount() const;

    /*! Create starting model by iterating over all regions.*/
    RVector createStartModel();

    /*! Fill vec with starting model values by iterating over all regions.*/
    void fillStartModel(RVector & vec);

    /*! DEPRECATED use setStartModel */
     void fillStartVector(RVector & vec){DEPRECATED;
         fillStartModel(vec);}

    /*! DEPRECATED use setStartModel */
     RVector createStartVector(){DEPRECATED;
         return createStartModel();}

    /*! Create and fill global model-weight vector */
    RVector createModelControl();

    /*! Fill global model-weight vector */
    void fillModelControl(RVector & vec);

    /*! Create and fill global constraints-weight vector */
    RVector createConstraintsWeight();

    /*! Fill global constraints-weight vector */
    void fillConstraintsWeight(RVector & vec);

    /*! Fill global constraints-matrix
        no regions: fill with 0th-order constraints */
    void fillConstraints(RSparseMapMatrix & C);

    /*! Syntactic sugar: set zweight/constraintType to all regions. */
    void setZWeight(double z);

    void setConstraintType(Index type);

    void fillBoundarySize(RVector & vec);

    inline const Mesh & paraDomain() const { return *paraDomain_; }

    inline Mesh & paraDomain() { return *paraDomain_; }

    std::vector < RVector3 > boundaryNorm() const;

//     RVector createFlatWeight(double zPower, double zWeight) const;

    /*!Permute all parameter marker after successful filled the
     * Regionmanager. */
    void permuteParameterMarker(const IVector & p);

    void createParaDomain_();

    void recountParaMarker_();

    TransCumulative < RVector > * transModel();

    void setLocalTransFlag(bool flag);

    bool haveLocalTrans() const { return haveLocalTrans_; }

protected:
    void copy_(const RegionManager & rm);

    /*! Fill \ref interRegionInterfaceMap_ */
    void findInterRegionInterfaces_();

    /*!
     * Internal method to create a region. The method is called from \ref setMesh()
     */
    Region * createRegion_(SIndex marker, const Mesh & mesh);

    /*!
     * Internal method to create a single parameter region. The method is called from \ref setMesh()
     */
    Region * createSingleRegion_(SIndex marker, const std::vector < Cell * > & cells);

    std::vector < SIndex > allRegionMarker_(bool exludeBoundary=false){
        std::vector < SIndex > tmp;
        for (std::map < SIndex, Region * > ::const_iterator
                it = regionMap_.begin(); it != regionMap_.end(); it++){
            if (!it->second->isBackground()) tmp.push_back(it->first);
        }
        return tmp;
    }


    bool verbose_;

    Index parameterCount_;

    Mesh * mesh_;
    Mesh * paraDomain_;

    std::map < SIndex, Region * > regionMap_;
    std::map< std::pair< SIndex, SIndex >, std::list < Boundary * > > interRegionInterfaceMap_;
    std::map< std::pair< SIndex, SIndex >, double > interRegionConstraints_;
    std::map< SIndex, double > interfaceConstraint_;

    double interRegionConstraintsZWeight_;

    TransCumulative < RVector > localTrans_;
    bool haveLocalTrans_;
    bool localTransHaveChanges_;
};


} // namespace GIMLI

#endif // _GIMLI_REGIONMANAGER__H
