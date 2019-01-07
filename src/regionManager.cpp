/******************************************************************************
 *   Copyright (C) 2008-2019 by the GIMLi development team                    *
 *   Carsten Rücker carsten@resistivity.net                                   *
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

#include "regionManager.h"

#include "mesh.h"
#include "node.h"
#include "shape.h"
#include "sparsematrix.h"
#include "stopwatch.h"
#include "trans.h"
#include "vectortemplates.h"


namespace GIMLI{

Region::Region(SIndex marker, RegionManager * parent, bool single)
    : marker_(marker), parent_(parent),
    isBackground_(false), isSingle_(single), parameterCount_(0)
    , tM_(NULL) {
    init_();
    if (isSingle_) {
        parameterCount_ = 1;
        constraintType_ = 0;
    }
}

Region::Region(SIndex marker, const Mesh & mesh, RegionManager * parent)
    : marker_(marker), parent_(parent)
        , isBackground_(false), isSingle_(false)
        , parameterCount_(0), tM_(NULL) {
    init_();
    this->resize(mesh);
}

Region::Region(const Region & region){
    copy_(region);
}

Region::~Region(){
   if (tM_ && ownsTrans_) delete tM_;
}

Region & Region::operator = (const Region & region){
    if (this != & region){
        copy_(region);
    } return *this;
}

void Region::init_() {
    constraintType_ = 1;
    startParameter_ = 0;
    endParameter_   = 0;
    parameterCount_ = 0;
    zPower_         = 0.0;
    zWeight_        = 1.0;
    fixValue_       = 0.0;
    lowerBound_     = 0.0;
    upperBound_     = 0.0;
    mcDefault_      = 1.0;
    startDefault_   = 1.0;
    ownsTrans_ = true;

    transString_    = "Log";
    tM_ = new Trans < RVector >;
}

void Region::copy_(const Region & region){
    THROW_TO_IMPL
    marker_ = region.marker();
    isBackground_ = region.isBackground();
    isSingle_ = region.isSingle();
}

void Region::setBackground(bool background){
    if (background != isBackground_) {
        markBackground(background);
        parent_->recountParaMarker_();
        parent_->createParaDomain_();
        bounds_.clear();
    }
}

void Region::setFixValue(double val){
    fixValue_ = val;
    markBackground(false);
    setBackground(true);
}

void Region::setSingle(bool single){
    if (single != isSingle_) {
        markSingle(single);
        parent_->recountParaMarker_();
        parent_->createParaDomain_();
        bounds_.clear();
    }
}

void Region::resize(const Mesh & mesh){
    cells_ = mesh.findCellByMarker(marker_);
    bounds_.clear();

    if (!isBackground_ && !isSingle_){

        Index nBounds = mesh.boundaryCount();
        if (nBounds == 0) {
            std::cerr << "WARNING! no boundaries defined! run mesh.createNeighbourInfos()" << std::endl;
        }
        bool leftParaId = false;
        bool rightParaId = false;

        for (Index i = 0; i < nBounds; i ++){
            leftParaId = false;
            rightParaId = false;

            if (mesh.boundary(i).marker() == 0){
                if (mesh.boundary(i).leftCell() != NULL) {
                    leftParaId = (mesh.boundary(i).leftCell()->marker() == marker_);
                }
                if (mesh.boundary(i).rightCell() != NULL){
                    rightParaId = (mesh.boundary(i).rightCell()->marker() == marker_);
                }

                if (leftParaId && rightParaId) bounds_.push_back(&mesh.boundary(i));
            }
        }
    }
//    std::cout <<WHERE_AM_I << " " << marker_ << " " << bounds_.size() << std::endl;
    if (zWeight_ != 1.0){
        fillConstraintsWeightWithFlatWeight();
    } else {
//         std::cout << "Region::resize(const Mesh & mesh)" <<  parameterCount_ << " " << constraintType_ << " " << this->constraintCount() << std::endl;
        constraintsWeight_.resize(constraintCount(), 1.0);
    }
}

void Region::resize(const std::vector < Cell * > & cells){
    cells_ = cells;
    bounds_.clear();
    if (zWeight_ != 1.0){
        fillConstraintsWeightWithFlatWeight();
    } else {
        constraintsWeight_.resize(constraintCount(), 1.0);
    }
}

void Region::countParameter(Index start){
//     __MS(marker()<< " " << fixValue())
    startParameter_ = start;
    if (isBackground_) {
        for (Index i = 0, imax = cells_.size(); i < imax; i ++) {
//             __MS(cells_[i]->marker() << " " << fixValue_)
            if (abs(fixValue_) > TOLERANCE) {
                if (cells_[i]->marker() >= -1) {
                    // set only for positive marker that not already fixed regions
                    cells_[i]->setMarker(MARKER_FIXEDVALUE_REGION - marker());
                }
            } else {
                cells_[i]->setMarker(-1);
            }
        }
        bounds_.clear();
        parameterCount_ = 0;
    } else if (isSingle_) {
        for (Index i = 0, imax = cells_.size(); i < imax; i ++) cells_[i]->setMarker(start);
        bounds_.clear();
        parameterCount_ = 1;
    } else {
        for (Index i = 0, imax = cells_.size(); i < imax; i ++) cells_[i]->setMarker(start + i);
        parameterCount_ = cells_.size();
    }

    // reset attributes to allow for clean mapping
    for (Index i = 0, imax = cells_.size(); i < imax; i ++) cells_[i]->setAttribute(0.0);

    endParameter_ = start + parameterCount_;
    modelControl_.resize(parameterCount_, mcDefault_);
    startModel_.resize(parameterCount_, startDefault_);
    //std::cout << WHERE_AM_I << " " << marker_ << " " << parameterCount_ << " " << startParameter_ << " " << endParameter_ <<  std::endl;
}

//################ Start values
void Region::setStartModel(const RVector & start){
    setBackground(false);
    if (start.size() == parameterCount_){
       startModel_ = start;
    } else {
        throwLengthError(1, WHERE_AM_I + " " + toStr(start.size()) + " != " + toStr(parameterCount_));
    }
}

void Region::setStartModel(double start){
    startDefault_ = start;
    this->setStartModel(RVector(parameterCount_, start));
}

void Region::fillStartModel(RVector & vec){
    if (isBackground_) return;
    if (startModel_.size() != parameterCount_){
        std::cerr << "WARNING! starting value for region " << marker_ << " not set. "<< std::endl;
    } else {
        if (isSingle_){
            vec[startParameter_] = startModel_[0];
        } else {
            for (Index i = 0, imax = cells_.size(); i < imax; i ++) {
                vec[cells_[i]->marker()] = startModel_[i];
            }
        }
    }
}

//################ Model behaviour

void Region::setModelControl(double val){
    if (val < TOLERANCE) val = 1.0;
    mcDefault_ = val;
    setBackground(false);
    modelControl_.resize(parameterCount_);
    modelControl_.fill(val);
}

void Region::setModelControl(const RVector & mc){
    setBackground(false);
    if (mc.size() == parameterCount_){
       modelControl_ = mc;
    } else {
        throwLengthError(1, WHERE_AM_I + " " + toStr(mc.size()) + " != " + toStr(parameterCount_));
    }
}

void Region::setModelControl(PosFunctor * mcF){
    setBackground(false);
    modelControl_.resize(parameterCount_);
    if (isSingle_){
        THROW_TO_IMPL
    }
    for (size_t i = 0; i < parameterCount_; i ++) modelControl_[i] = (*mcF)(cells_[i]->center());
}

void Region::fillModelControl(RVector & vec){
    if (isBackground_) return;
    if (isSingle_){
        vec[startParameter_] = modelControl_[0];
    } else{
        for (Index i = 0, imax = cells_.size(); i < imax; i ++) {
            vec[cells_[i]->marker()] = modelControl_[i];
        }
    }
}

//################ constraints behaviour
Index Region::constraintCount() const {
    if (isSingle_ && constraintType_ == 1) return 0;

    if (constraintType_ == 0 || constraintType_ == 2 || constraintType_ == 20) return parameterCount();
    if (constraintType_ == 10) return bounds_.size() + parameterCount();
    return bounds_.size();
}

void Region::fillConstraints(RSparseMapMatrix & C, Index startConstraintsID){
    if (isBackground_) return;

    if (isSingle_ && constraintType_ == 1) return;

    double cMixRatio = 1.0; // for mixing 1st or 2nd order with 0th order (constraintTypes 10 and 20)
    if (constraintType_ == 10 || constraintType_ == 20) cMixRatio = 1.0; //**retrieve from properties!!!

    if (constraintType_ == 0 || constraintType_ == 20){ //purely 0th or mixed 2nd+0th
        for (size_t i = 0; i < parameterCount_; i++) {
            C[startConstraintsID + i][startParameter_ + i] = cMixRatio;
        }
        if (constraintType_ == 0) return;
    }

    SIndex leftParaId = -1, rightParaId = -1;
    if (constraintType_ == 2 || constraintType_ == 20) { //** 2nd order constraints (opt. mixed with 0th)
        for (std::vector < Boundary * >::iterator it = bounds_.begin(), itmax = bounds_.end();
            it != itmax; it ++){
            leftParaId = -1;
            rightParaId = -1;
            if ((*it)->leftCell() ) leftParaId  = (*it)->leftCell()->marker();
            if ((*it)->rightCell()) rightParaId = (*it)->rightCell()->marker();
            if (leftParaId >= (int)startParameter_ && leftParaId < (int)endParameter_ &&
                 rightParaId >= (int)startParameter_ && rightParaId < (int)endParameter_ &&
                    leftParaId != rightParaId){
                C[leftParaId][rightParaId] = -1;
                C[rightParaId][leftParaId] = -1;
                C[leftParaId][leftParaId] += 1;
                C[rightParaId][rightParaId] += 1;
            }
        }
        return;
    }
    //** 1st order constraints (opt. combined with 0th order)
    Index cID = startConstraintsID;
    for (std::vector < Boundary * >::iterator it = bounds_.begin(), itmax = bounds_.end();
        it != itmax; it ++){

        leftParaId = -1;
        rightParaId = -1;
        if ((*it)->leftCell()) leftParaId  = (*it)->leftCell()->marker();
        if ((*it)->rightCell()) rightParaId = (*it)->rightCell()->marker();

        //if ( leftParaId > -1 && rightParaId > -1 && leftParaId != rightParaId) {
        if (leftParaId >= (int)startParameter_ && leftParaId < (int)endParameter_ &&
            rightParaId >= (int)startParameter_ && rightParaId < (int)endParameter_ &&
                leftParaId != rightParaId){
            C[cID][leftParaId] = 1;
            C[cID][rightParaId] = -1;
        }
        cID ++;
    }
    if (constraintType_ == 10) { //** combination with 0th order
        for (size_t i = 0; i < parameterCount_; i++) {
            C[cID][startParameter_ + i] = cMixRatio;
            cID++;
        }

    }
}

void Region::setConstraintType(Index type) {
    constraintType_ = type;
}

void Region::setConstraintsWeight(double val){
    setBackground(false);
    zWeight_ = 1.0;
    //std::cout << "Region::setConstraintsWeight(double val) " << val << " " <<  this->constraintCount() << std::endl;
    constraintsWeight_.resize(this->constraintCount());
    constraintsWeight_.fill(val);
}

void Region::setConstraintsWeight(const RVector & sw){
    //std::cout << "Region::setConstraintsWeight(const RVector & sw) " << sw.size() << " " <<  this->constraintCount() << std::endl;
    setBackground(false);
    if (sw.size() == this->constraintCount()){
        constraintsWeight_ = sw;
        zWeight_ = 1.0;
    } else {
        throwLengthError(1, WHERE_AM_I + " " + toStr(sw.size()) + " != " + toStr(constraintCount()));
    }
}

void Region::createConstraintsWeight_(){
    if (constraintsWeight_.size() != this->constraintCount()){
        constraintsWeight_.resize(this->constraintCount(), 1);
    }
    if (zWeight_ != 1.0) this->fillConstraintsWeightWithFlatWeight();
}

void Region::fillConstraintsWeight(RVector & vec, Index constraintStart){
    if (isBackground_) return;
    this->createConstraintsWeight_();

    for (Index i = 0, imax = constraintCount(); i < imax; i ++) {
        //std::cout << i << " " << constraintsWeight_[i] << std::endl;
        vec[constraintStart + i] = constraintsWeight_[i];
    }
}

void Region::fillConstraintsWeightWithFlatWeight(){
    if (isBackground_ || isSingle_ || (constraintType() == 0) || (constraintType() == 2)) return;
    constraintsWeight_.resize(constraintCount(), 1.0);

    for (Index i = 0, imax = bounds_.size(); i < imax; i ++) {

        double zDir = std::fabs(bounds_[i]->norm()[parent_->mesh().dim() -1]); //! z-component
//        double horizDir = std::sqrt(1.0 - zDir * zDir); //! horizontal component

        if (zPower_ != 0.0){ //! zPower controls and zweight is minimum zweight
            __MS("DEPRECATED")
            constraintsWeight_[i] = max(zWeight_, std::pow(1.0 - zDir, zPower_));
        } else { //! zweight controls and is power factor, includes
            //** for a temporary back conversion hack: this was the old zWeight function:
            //constraintsWeight_[i] = 1.0 + zWeight_ - zDir; //! rather linear for bigger angles

//            constraintsWeight_[i] = std::pow(zWeight_, std::fabs(zDir / 3.0)); //! rather logarithmically
            //** this is the new and better one:
            constraintsWeight_[i] = 1.0 + zDir * (zWeight_ - 1.0); //! rather linear for bigger angles
        }
    }
}

void Region::fillBoundaryNorm(std::vector< RVector3 > & vnorm, Index boundStart){
    if (isBackground_ || isSingle_ || (constraintType() == 0)) return;

    for (Index i = 0, imax = bounds_.size(); i < imax; i ++) {
        vnorm[i + boundStart] = bounds_[i]->norm();
    }
}

void Region::fillBoundarySize(RVector & vec, Index boundStart){
    if (isBackground_ || isSingle_ || (constraintType() == 0)) return;

    for (Index i = 0, imax = bounds_.size(); i < imax; i ++) {
        vec[i + boundStart] = bounds_[i]->shape().domainSize();
    }
}

void Region::setModelTransformation(const Trans< RVector > & tM){
THROW_TO_IMPL

}

void Region::setTransModel(Trans< RVector > & tM){
    if (tM_ && ownsTrans_) delete tM_;
    tM_ = & tM;
    parent_->setLocalTransFlag(true);
    ownsTrans_ = false;
}

void Region::setModelTransStr_(const std::string & val){
    transString_ = val;
    delete tM_; tM_ = NULL;

    if (val == "lin" || val == "Lin"){
        tM_ = new Trans< RVector >();
    } else if (val == "log" || val == "Log"){
        tM_ = new TransLogLU< RVector >(lowerBound_, upperBound_);
    } else if (val == "cot" || val == "Cot" || val == "tan" || val == "Tan"){
        tM_ = new TransCotLU< RVector >(lowerBound_, upperBound_);
    } else  {
        throwLengthError(1, WHERE_AM_I + " : " + val + ". Available are: lin, log, cot/tan.");
    }
    parent_->setLocalTransFlag(true);

    ownsTrans_ = true;
}

void Region::setLowerBound(double lb){
    lowerBound_ = lb;
    setModelTransStr_(transString_);
}

void Region::setUpperBound(double ub){
    upperBound_ = ub;
    setModelTransStr_(transString_);
}

void Region::setParameters(double start, double lb, double ub, std::string transString){
    if (lb < ub) {
        if ((start <= lb) | (start >= ub)) {
            std::cout << "WARNING! starting model not within bounds! readjusting" << std::endl;
            setStartModel(std::sqrt(lb * ub));
        } else {
            setStartModel(start);
        }
        lowerBound_ = lb;
        upperBound_ = ub;
        if (transString.size() > 0) { // any given
            setModelTransStr_(transString);
        } else { // use preset otherwise
            setModelTransStr_(transString_);
        }
    } else {
        throwError(EXIT_FAILURE, WHERE_AM_I + " bounds not matching: " + toStr(lb) + ">=" + toStr(ub));
    }
}

//******************************************************************************
RegionManager::RegionManager(bool verbose) : verbose_(verbose), mesh_(NULL){
    paraDomain_ = new Mesh();
    parameterCount_ = 0;
    haveLocalTrans_ = false;
    localTransHaveChanges_ = true;
    interRegionConstraintsZWeight_ = 1.0;
}

RegionManager & RegionManager::operator = (const RegionManager & rm){
    if (this != &rm){
        copy_(rm);
    } return *this;
}

void RegionManager::copy_(const RegionManager & rm){
    CERR_TO_IMPL
}

RegionManager::~RegionManager(){
    clear();
    if (paraDomain_) delete paraDomain_;
}

const Mesh & RegionManager::mesh() const {
    if (mesh_== 0){
        throwError(1, "RegionManager knows no mesh.");
    }
    return *mesh_;
}

Region * RegionManager::region(SIndex marker){
    if (regionMap_.count(marker) == 0){
        throwError(EXIT_DEFAULT, WHERE_AM_I + " no region with marker " + toStr(marker));
    } else {
        return regionMap_[marker];
    }
    return NULL;
}

void RegionManager::clear(){
    for (std::map< SIndex, Region* >::const_iterator it = regionMap_.begin(), end = regionMap_.end();
          it != end; it ++){
        delete it->second;
    }
    regionMap_.clear();
    interRegionInterfaceMap_.clear();
    interRegionConstraints_.clear();
    interfaceConstraint_.clear();

    if (paraDomain_) { paraDomain_->clear(); }
    if (mesh_) { delete mesh_; mesh_ = NULL; }
}

void RegionManager::setMesh(const Mesh & mesh, bool holdRegionInfos){

    if (!holdRegionInfos){
        if (verbose_) std::cout << "Reset region parameter" << std::endl;
        this->clear();
    }

    Stopwatch swatch(true);
    if (verbose_) std::cout << "RegionManager copying mesh ...";

    if (mesh_) delete mesh_;
    mesh_ = new Mesh(mesh);

    if (verbose_){
        std::cout << swatch.duration(true) << " s " << std::endl;
        std::cout << "create NeighbourInfos ... ";
    }
    mesh_->createNeighbourInfos();
    if (verbose_) {
        std::cout << swatch.duration(true) << " s " << std::endl;
        std::cout << "analysing mesh ... ";
    }

    //** looking for and create regions
    IVector regions(unique(sort(mesh_->cellMarkers())));

    if (verbose_) std::cout << regions.size() << " regions." << std::endl;

    bool singleOnly = false;
    if (regions.size() > 50){
        singleOnly = true;
        log(Info,"More than 50 regions so we assume singles only regions.");
    }

    std::map < int, std::vector< Cell * > > markerCellVectorMap;
    if (singleOnly){
        for(std::vector< Cell * >::const_iterator it = mesh_->cells().begin();
                                               it != mesh_->cells().end(); it++){
            if (! markerCellVectorMap.count((*it)->marker())){
                markerCellVectorMap[(*it)->marker()] = std::vector< Cell * >();
            }
            markerCellVectorMap[(*it)->marker()].push_back((*it));
        }
    }

    for (Index i = 0; i < regions.size(); i ++){
        if (singleOnly){
            createSingleRegion_(regions[i], markerCellVectorMap[regions[i]]);
        } else {
            createRegion_(regions[i], *mesh_);
        }
    }
    //** looking for and create region interfaces

    //** looking for and create inter-region interfaces
    this->findInterRegionInterfaces_();

    if (singleOnly){
        log(Info, "Applying *:* interregion constraints.");
        for (Index i = 0; i < regions.size(); i ++){
            for (Index j = 0; j < regions.size(); j ++){
                if (i != j){
                    setInterRegionConstraint(i, j, 1.0);
                                // if (verbose_) std::cout << minRegion[i] << " <-> "
                                //                         << maxRegion[j] << " weight:"
                                //                         << toDouble(row[2]) << std::endl;
                }
            }
        }
    }
    this->recountParaMarker_();
    this->createParaDomain_();
}

Region * RegionManager::createSingleRegion_(SIndex marker, const std::vector < Cell * > & cells){
    Stopwatch swatch(true);
    Region * region = NULL;
    if (regionMap_.count(marker) == 0){
        region = new Region(marker, this, true);
        regionMap_.insert(std::make_pair(marker, region));
    } else {
        THROW_TO_IMPL
        region = regionMap_[marker];
    }

    region->resize(cells);
    return region;
}

Region * RegionManager::createRegion_(SIndex marker, const Mesh & mesh){
    Region * region = NULL;

    if (regionMap_.count(marker) == 0){
        region = new Region(marker, mesh, this);
        regionMap_.insert(std::make_pair(marker, region));
    } else {
        region = regionMap_[marker];
        region->resize(mesh);
        //std::cerr << WHERE_AM_I << " Region with marker " << marker << " already exists." << std::endl;
    }
    return region;
}

Region * RegionManager::addRegion(SIndex marker, const Mesh & mesh){
    Region * region = createRegion_(marker, mesh);
    recountParaMarker_(); //** make sure the counter is right
    return region;
}

void RegionManager::createParaDomain_(){
    if (verbose_) std::cout << "creating para domain ... ";
    Stopwatch swatch(true);
    IndexArray cellIdx;
    cellIdx.reserve(mesh_->cellCount());

    for (std::map< SIndex, Region* >::const_iterator
         it = regionMap_.begin(), end = regionMap_.end(); it != end; it ++){

        if (!it->second->isBackground()){
            for (std::vector < Cell * >::const_iterator itc = it->second->cells().begin();
                 itc != it->second->cells().end(); itc++){
                cellIdx.push_back((*itc)->id());
            }
//             std::transform(it->second->cells().begin(), it->second->cells().end(),
//                             std::back_inserter(cellIdx), std::mem_fun(&Cell::id));
        }
    }

    paraDomain_->createMeshByCellIdx(*mesh_, cellIdx);
    if (verbose_) std::cout << swatch.duration(true) << " s" << std::endl;
}

void RegionManager::recountParaMarker_(){
    Index count = 0;

    for (std::map< SIndex, Region* >::const_iterator it = regionMap_.begin(),
        end = regionMap_.end();
          it != end; it ++){
//         __MS(it->second->marker() << "  " << it->second->fixValue())
        it->second->countParameter(count);
        count += it->second->parameterCount();
    }
    if (verbose_) std::cout << "Recounted parameter: " << count << std::endl;
}

void RegionManager::findInterRegionInterfaces_(){
    interRegionInterfaceMap_.clear();
    Boundary * bound;
    //** having fun with stl
    std::map< std::pair< SIndex, SIndex >, std::list < Boundary * > > ::iterator iRMapIter;

    for (std::vector < Boundary * >::const_iterator
            iBound  = mesh_->boundaries().begin();
            iBound != mesh_->boundaries().end(); iBound ++){
        bound = (*iBound);
        if (bound->leftCell() && bound->rightCell()){
            if (bound->leftCell()->marker() != bound->rightCell()->marker()){
                SIndex minMarker = min(bound->leftCell()->marker(), bound->rightCell()->marker());
                SIndex maxMarker = max(bound->leftCell()->marker(), bound->rightCell()->marker());

                iRMapIter = interRegionInterfaceMap_.find(std::pair< SIndex, SIndex >(minMarker, maxMarker));
                if (iRMapIter == interRegionInterfaceMap_.end()){
                    interRegionInterfaceMap_.insert(std::pair< std::pair< SIndex, SIndex >,
                                      std::list < Boundary * > > (std::pair< SIndex, SIndex>(minMarker, maxMarker),
                                                                  std::list< Boundary* >()));
                }
                interRegionInterfaceMap_[std::pair< SIndex, SIndex > (minMarker, maxMarker)].push_back(bound);
            }
        }
    }
    // if (verbose_ && interRegionInterfaceMap_.size()){
    //     for (std::map< std::pair< SIndex, SIndex >, std::list < Boundary * > >::iterator
    //             it  = interRegionInterfaceMap_.begin();
    //             it != interRegionInterfaceMap_.end(); it ++){
    //         std::cout << "(" << it->first.first << "," << it->first.second << ") "
    //                          << it->second.size() << " boundaries. " << std::endl;
    //     }
    // }
}

void RegionManager::fillStartModel(RVector & vec){
    if (vec.size() != parameterCount()) vec.resize(parameterCount());
    for (std::map< SIndex, Region* >::const_iterator it = regionMap_.begin(), end = regionMap_.end();
          it != end; it ++){
        it->second->fillStartModel(vec);
    }
}

RVector RegionManager::createStartModel(){
    RVector vec(parameterCount());
    fillStartModel(vec);
    return vec;
}

void RegionManager::fillModelControl(RVector & vec){
    //!** no regions: fill 0th-order constraints
    if (regionMap_.empty()){
        vec.resize(this->parameterCount(), 1.0);
        return;
    }

    if (vec.size() != parameterCount()) vec.resize(parameterCount(), 1.0);

    for (std::map< SIndex, Region* >::const_iterator it = regionMap_.begin(), end = regionMap_.end();
          it != end; it ++){
        it->second->fillModelControl(vec);
    }
}

RVector RegionManager::createModelControl(){
    RVector vec(parameterCount(), 1.0);
    fillModelControl(vec);
    return vec;
}

RVector RegionManager::createConstraintsWeight(){
    RVector vec(constraintCount(), 1.0);
    fillConstraintsWeight(vec);
    return vec;
}

void RegionManager::fillConstraintsWeight(RVector & vec){

    //!** no regions: fill 0th-order constraints
    if (regionMap_.empty()){
        vec.resize(this->parameterCount(), 1.0);
        return;
    }

    if (vec.size() != constraintCount()) vec.resize(constraintCount(), 1.0);

    //!** fill constraints weights from individual regions
    Index cID = 0;
    for (std::map< SIndex, Region* >::const_iterator
            it  = regionMap_.begin(); it != regionMap_.end(); it ++){
        it->second->fillConstraintsWeight(vec, cID);
        cID += it->second->constraintCount();
    }

    //!** fill constraints weights from inter regions constrains
    if (interRegionConstraints_.size() > 0){
        if (verbose_) std::cout << "Applying inter region constraints weights." << std::endl;

        for (std::map< std::pair< SIndex, SIndex >, double >::iterator
                        it = interRegionConstraints_.begin();
                        it != interRegionConstraints_.end(); it ++){

            std::map< std::pair< SIndex, SIndex >, std::list < Boundary * > >::const_iterator iRMapIter = interRegionInterfaceMap_.find(it->first);

//             if (verbose_) {
//                 std::cout << "regions: " << it->first.first << " " << it->first.second;
//                 if (iRMapIter != interRegionInterfaceMap_.end()){
//                     std::cout << " (" << iRMapIter->second.size() << ") " ;
//                 }
//                 std::cout << " = " << it->second << std::endl;
//             }

            if (iRMapIter != interRegionInterfaceMap_.end()){
                for (std::list < Boundary * >::const_iterator  bIter  = iRMapIter->second.begin();
                                                                bIter != iRMapIter->second.end(); bIter++){

                    //** interfaceConstraint_ overwrite interRegionConstraints_
                    if (interfaceConstraint_.count((*bIter)->marker())) {
                        vec[cID] = interfaceConstraint_[(*bIter)->marker()];
                    } else {
                        if (std::fabs(interRegionConstraintsZWeight_ - 1.0) < TOLERANCE){
                            vec[cID] = it->second;
                        } else {
                            RVector3 meanNorm(0.0, 0.0, 0.0);
//                             for (std::list < Boundary * >::const_iterator
//                                 boundIter  = iRMapIter->second.begin();
//                                 boundIter != iRMapIter->second.end(); boundIter++){
//                                 meanNorm += (*boundIter)->norm();
//                             }

                            RVector3 left(0.0, 0.0, 0.0);
                            for (std::vector < Cell * >::iterator
                                cIt = regionMap_.find(it->first.first)->second->cells().begin();
                                cIt != regionMap_.find(it->first.first)->second->cells().end();
                                    cIt ++){
                                left += (*cIt)->center();
                            }
                            left /= regionMap_.find(it->first.first)->second->cells().size();

                            RVector3 right(0.0, 0.0, 0.0);
                            for (std::vector < Cell * >::iterator
                                cIt = regionMap_.find(it->first.second)->second->cells().begin();
                                cIt != regionMap_.find(it->first.second)->second->cells().end();
                                cIt ++){
                                right += (*cIt)->center();
                            }
                            right /= regionMap_.find(it->first.second)->second->cells().size();
                            meanNorm = left-right;

                            meanNorm.normalise();
                            //std::cout << meanNorm.abs() << std::endl;
                            double zDir = std::fabs(meanNorm[mesh_->dim() -1]); //! z-component

                            //! rather linear for bigger angles
                            vec[cID] = (1.0 + (interRegionConstraintsZWeight_ - 1.0) * zDir)
                                            * it->second;

                        }

                    }
                    cID ++;

                    if (regionMap_.find(it->first.first )->second->isSingle() &&
                         regionMap_.find(it->first.second)->second->isSingle()){
                        break;
                    }
                }
            }
        }
    } // if (interRegionConstraints_.size() > 0)

//     if (interfaceConstraintMap_
//     interfaceConstraintMap_[it->second->boundaries()[i]->marker()] = toDouble(row[1]);
}

void RegionManager::fillBoundarySize(RVector & vec){
    if (vec.size() != constraintCount()) vec.resize(constraintCount(), 1.0);
    Index boundCount = 0;
    for (std::map< SIndex, Region* >::const_iterator it = regionMap_.begin(), end = regionMap_.end();
          it != end; it ++){
        it->second->fillBoundarySize(vec, boundCount);
        boundCount += it->second->constraintCount();
    }

}

Index RegionManager::parameterCount() const {
    if (regionMap_.empty()) {
        if (parameterCount_ == 0){
            return parameterCount_;

            // zero should be possible
            throwLengthError(1, WHERE_AM_I + " neither region defined nor parameterCount set.");
        }
        return parameterCount_;
    }

    Index count = 0;
    for (std::map< SIndex, Region* >::const_iterator
        it = regionMap_.begin(), end = regionMap_.end(); it != end; it ++){
        count += it->second->parameterCount();
    }
    return count;
// template <class Pair> struct select2nd {
// typename Pair::second_type operator () (const Pair &x) const
// { return x.second;}
// };
//std::plus<int>
//     return std::accumulate(regions_.begin(), regions_.end(), 0,
//                             std::mem_fun(&Region::modelCount));
}

Index RegionManager::constraintCount() const {
    if (regionMap_.empty()) {
        return parameterCount_;
    }

    Index count = 0;
    for (std::map< SIndex, Region* >::const_iterator it = regionMap_.begin(), end = regionMap_.end();
          it != end; it ++){
        count += it->second->constraintCount();
//         std::cout << count << std::endl;
    }
//     std::cout << count << " " <<  interRegionConstraintsCount() << std::endl;

    return count + interRegionConstraintsCount();
}

Index RegionManager::interRegionConstraintsCount() const {
    Index count = 0;

    for (std::map< std::pair< SIndex, SIndex >, double >::const_iterator
            it = interRegionConstraints_.begin();
            it != interRegionConstraints_.end(); it ++){

        if (regionMap_.find(it->first.first )->second->isSingle() &&
             regionMap_.find(it->first.second)->second->isSingle()){
             count += 1;
        } else {
            if (interRegionInterfaceMap_.find(it->first) != interRegionInterfaceMap_.end()){
                count += interRegionInterfaceMap_.find(it->first)->second.size();
            }
        }
    }

    return count;
}

void RegionManager::fillConstraints(RSparseMapMatrix & C){
//     __MS(&C)
//     __MS(C.rtti())

    Index nModel  = parameterCount();
    Index nConstr = constraintCount();

    C.clear();

    //!** no regions: fill 0th-order constraints
    if (regionMap_.empty() || nConstr == 0){
        C.setCols(nModel);
        C.setRows(nModel);

        for (size_t i = 0; i < parameterCount(); i++) C[i][i] = 1;
        return;
    }

    C.setRows(nConstr);
    C.setCols(nModel);

    Index consCount = 0;

    for (std::map< SIndex, Region* >::const_iterator it = regionMap_.begin(), end = regionMap_.end();
          it != end; it ++){
        it->second->fillConstraints(C, consCount);
        consCount += it->second->constraintCount();
    }

    if (interRegionConstraints_.size() > 0){
        if (verbose_) std::cout << "Applying inter region constraints." << std::endl;

        for (std::map< std::pair< SIndex, SIndex >, double >::iterator
                it  = interRegionConstraints_.begin();
                it != interRegionConstraints_.end(); it ++){
//             if (verbose_) std::cout << "regions: " << it->first.first << " "
//                                       << it->first.second << std::endl;
            std::map< std::pair< SIndex, SIndex >, std::list < Boundary * > >::const_iterator iRMapIter;
            iRMapIter = interRegionInterfaceMap_.find(it->first);

            Region * lR = region(it->first.first);
            Region * rR = region(it->first.second);
            size_t lStart = lR->startParameter();
            size_t rStart = rR->startParameter();
            RVector *lMC = lR->modelControl();
            RVector *rMC = rR->modelControl();
            if (lMC->size() == 0 || rMC->size() == 0){
                __MS(lR->parameterCount() << "  " << rR->parameterCount())
                throwLengthError(1, WHERE_AM_I + " left | right  MC size == 0 " + toStr(lMC->size())
                + " "+ toStr(rMC->size()));
            }
//             if (verbose_){
//                 std::cout << "left : " << lMC->size() << " min " << min(*lMC) << " max " << max(*lMC) << " start " << lStart << std::endl;
//                 std::cout << "right: " << rMC->size() << " min " << min(*rMC) << " max " << max(*rMC) << " start " << rStart << std::endl;
//             }

            if (iRMapIter != interRegionInterfaceMap_.end()){
                for (std::list < Boundary * >::const_iterator bIter  = iRMapIter->second.begin();
                                                               bIter != iRMapIter->second.end(); bIter++){

                    size_t lMarker = (*bIter)->leftCell()->marker();
                    size_t rMarker = (*bIter)->rightCell()->marker();
                    if (lMarker < lStart || rMarker < rStart){ //** interchange left/right
                        size_t dummy = lMarker;
                        lMarker = rMarker;
                        rMarker = dummy;
                    }
//                    std::cout << lMarker << " " << rMarker << " " << lMarker - lStart << " " << rMarker - rStart << std::endl;

                    C[consCount][lMarker] = +1.0 / (*lMC)[size_t(lMarker - lStart)];
                    C[consCount][rMarker] = -1.0 / (*rMC)[size_t(rMarker - rStart)];
                    consCount ++;

                    if (regionMap_.find(it->first.first )->second->isSingle() &&
                         regionMap_.find(it->first.second)->second->isSingle()){
                        break;
                    }
                }
            }
        }
    }
//     __MS(C.rtti())
}

std::vector < RVector3 > RegionManager::boundaryNorm() const {
    std::vector < RVector3 > vnorm(constraintCount(), RVector3(0.0, 0.0, 0.0));
    Index consCount = 0;

    for (std::map< SIndex, Region* >::const_iterator it = regionMap_.begin(), end = regionMap_.end();
          it != end; it ++){
        it->second->fillBoundaryNorm(vnorm, consCount);
        consCount += it->second->constraintCount();
    }
    return vnorm;
}

// RVector RegionManager::createFlatWeight(double zPower, double zWeight) const {
//     DEPRECATED
//     RVector tmp(constraintCount());
//     RVector normDir;
//     if (mesh_->dim() == 2){
//         normDir = y(boundaryNorm());
//     } else {
//         normDir = z(boundaryNorm());
//     }
//     for (size_t i = 0, imax = normDir.size(); i < imax; i++) {
//         if (zPower != 0) { //! zpower controls, zweight is minimal zweight
//             tmp[i] = std::pow(1.0 - std::fabs(normDir[i]), zPower);
//             if (tmp[i] < zWeight) tmp[i] = zWeight;
//         } else if (zWeight != 1.0) { //! zweight controls and is power factor, includes
//             tmp[i] = std::pow(std::fabs(normDir[i]), zWeight);
//         }
//     }
//     return tmp;
// }

void RegionManager::loadMap(const std::string & fname){

    std::map< std::string, void (Region::*)(const std::string & val) > regionAttributeMap;
    regionAttributeMap[lower("MC")]       = &Region::setModelControlStr_;
    regionAttributeMap[lower("start")]    = &Region::setStartModelStr_;
    regionAttributeMap[lower("zPower")]   = &Region::setZPowerStr_;
    regionAttributeMap[lower("zWeight")]  = &Region::setZWeightStr_;
    regionAttributeMap[lower("fix")]      = &Region::setFixValueStr_;
    regionAttributeMap[lower("Ctype")]    = &Region::setConstraintTypeStr_;
    regionAttributeMap[lower("Trans")]    = &Region::setModelTransStr_;
    regionAttributeMap[lower("lBound")]   = &Region::setLowerBoundStr_;
    regionAttributeMap[lower("uBound")]   = &Region::setUpperBoundStr_;
    regionAttributeMap[lower("single")]   = &Region::setSingleStr_;
    regionAttributeMap[lower("background")] = &Region::setBackgroundStr_;

    if (verbose_) std::cout << "Reading region control file" << std::endl;

    std::fstream file; openInFile(fname, &file);
    std::vector < std::string > token;
    std::vector < std::string > row;

    while (!file.eof()){

        row = getNonEmptyRow(file, '@');
        if (row.empty()){
            continue;
//             file.close();
//             return;
        }

        if (row[0][0] == '#'){

            token.clear();

            if (row[0].size() > 1) {
                // tokenline starts with '#XYZ'
                token.push_back(row[0].substr(1, row[0].size()));
            }

            for (Index i = 1; i < row.size(); i ++) {
                token.push_back(row[i]);
            }

            // read next row line to parse them
            continue;
        }

        if (token.size() == 0) {
             std::cerr << WHERE_AM_I << " not a valid region file. looking for leading '#' in " << fname << std::endl;
             file.close();
             return;
        }

        // interpret row as region informations
        if (lower(token[0]) == "no"){

            if (verbose_){
                if (verbose_) std::cout << "Get region property tokens: " << std::endl;
                for (Index i = 0; i < token.size(); i ++) {
                    if (verbose_) std::cout << token[i] << ", ";
                }
                if (verbose_) std::cout << std::endl;
            }

            if (token.size() >= row.size()){
                std::vector < SIndex > regionMarker;

                if (row[0] != "*"){
                    //std::cout << row[0] << std::endl;
                    if (this->regionExists(toInt(row[0]))) {
                        regionMarker.push_back(toInt(row[0]));
                    } else {
                        std::cerr << "Region number " << toInt(row[0]) << " " << row[0] << " does not exist!" << std::endl;
                    }
                } else {
                    regionMarker = allRegionMarker_();
                }

                for (Index j = 0; j < regionMarker.size(); j ++){
                    for (Index i = 1; i < row.size(); i ++){
                        if (regionAttributeMap.count(lower(token[i]))){
                            if (verbose_) {
                                if (j < 5){
                                    std::cout << regionMarker[j] << " : " << token[i]
                                                                 << " " << row[i] << std::endl;
                                } else if (j == 5 && j < regionMarker.size()-1){
                                    std::cout << "..." << std::endl;
                                } else if (j == regionMarker.size()-1){
                                    std::cout << regionMarker[j] << " : " << token[i]
                                                                 << " " << row[i] << std::endl;
                                }
                                
                            }
                            (region(regionMarker[j])->*regionAttributeMap[lower(token[i])])(row[i]);
                        } else {
                            std::cerr << WHERE_AM_I << " no region attribute associated with key:"
                                    << token[i] << std::endl;
                            std::cerr << "valid keys are: " << std::endl;
                            for (std::map< std::string, void (Region::*)(const std::string &) >::iterator it =
                                regionAttributeMap.begin(); it != regionAttributeMap.end(); it ++){
                                std::cerr << "\t" << it->first << std::endl;
                            }
                        }
                    }
                }
            } else {
                std::cerr << WHERE_AM_I << " too few tokens defined in region control file: " << fname << std::endl;
            }
        } else if (lower(token[0]) == "inter-region" || lower(token[0]) == "interregion") {
            if (verbose_){
                if (verbose_) std::cout << "Get inter-region properties" << std::endl;
                for (Index i = 0; i < token.size(); i ++) {
                    if (verbose_) std::cout << token[i] << " ";
                }
                if (verbose_) std::cout << std::endl;
            }

            if (row.size() == 3){
                std::vector < SIndex > minRegion;
                std::vector < SIndex > maxRegion;
                if (row[0] != "*"){
                    minRegion.push_back(toInt(row[0]));
                } else {
                    minRegion = allRegionMarker_(true);
                }
                if (row[1] != "*"){
                    maxRegion.push_back(toInt(row[1]));
                } else {
                    maxRegion = allRegionMarker_(true);
                }

                for (Index i = 0; i < minRegion.size(); i ++){
                    for (Index j = 0; j < maxRegion.size(); j ++){
                        if (minRegion[i] != maxRegion[j]){
                            setInterRegionConstraint(minRegion[i], maxRegion[j], toDouble(row[2]));
                            // if (verbose_) std::cout << minRegion[i] << " <-> "
                            //                         << maxRegion[j] << " weight:"
                            //                         << toDouble(row[2]) << std::endl;
                        }
                    }
                }
            } else {
                std::cerr << WHERE_AM_I << " too few tokens defined for inter-region constraints 3 != " << row.size() << std::endl;
            }

        } else if (lower(token[0]) == "interface"){
            if (verbose_){
                std::cout << "Applying interface properties." << std::endl;
                std::cout << "WARNING! no inner interfaces yet." << std::endl;
            }

            if (row.size() == 2){ // interface constraint
                if (row[0] == "*"){
                    THROW_TO_IMPL
                    for (std::map < SIndex, Region * > ::const_iterator it  = regionMap_.begin();
                                                          it != regionMap_.end(); it++){
                        for (Index i = 0; i < it->second->boundaries().size(); i ++){
                            if (verbose_) std::cout << it->second->boundaries()[i]->marker() << std::endl;
                            if (it->second->boundaries()[i]->marker() != 0){
                                interfaceConstraint_[it->second->boundaries()[i]->marker()] = toDouble(row[1]);
                            }
                        }
                    }
                } else {
                    interfaceConstraint_[toInt(row[0])] = toDouble(row[1]);
                }
            } else {
                for (Index i = 0; i < row.size(); i ++){
                    if (verbose_) std::cout << row[i] << " ";
                }
                if (verbose_) std::cout << std::endl;
                std::cerr << "Format unknown: (interfaceNo constraint) " << std::endl;
            }

        } else {
            std::cerr << WHERE_AM_I << " cannot interpret 1.st token: " << token[0] << std::endl;
            std::cerr << "Available are: " << std::endl
                      << "#no" << std::endl
                      << "#interface" << std::endl
                      << "#inter-region" << std::endl;
            row = getRow(file); if (row.empty()) continue;
        }
    }
    file.close();

    this->recountParaMarker_();
    this->createParaDomain_();
}

void RegionManager::setInterRegionConstraint(SIndex aIn, SIndex bIn, double c){
    SIndex a = min(aIn, bIn);
    SIndex b = max(aIn, bIn);

    if (a == b){
        std::cerr << WHERE_AM_I << " ignoring inter-region constraints "
                << a << " == " << b << std::endl;
    } else {
        if (interRegionInterfaceMap_.find(std::pair< SIndex, SIndex >(a, b))
             != interRegionInterfaceMap_.end()){
            interRegionConstraints_[std::pair< SIndex, SIndex >(a, b)] = c;
            if (debug()){
               std::cout << "regions: " << a << "<->" << b << " " << c << std::endl;
           }
        }
    }
}

void RegionManager::saveMap(const std::string & fname){
    THROW_TO_IMPL
}

void RegionManager::setLocalTransFlag(bool flag) {
    haveLocalTrans_ = flag;
    localTransHaveChanges_ = true;
}

TransCumulative < RVector > * RegionManager::transModel(){
    if (regionMap_.empty()){
        return NULL;
    }

    if (localTransHaveChanges_) localTrans_.clear();

    if (localTrans_.size() != allRegionMarker_(true).size()){

        for (std::map< SIndex, Region* >::const_iterator
             it = regionMap_.begin(); it != regionMap_.end(); it ++){

            if (!it->second->isBackground()){

                localTrans_.add(*it->second->transModel(),
                                it->second->startParameter(),
                                it->second->endParameter());
            }
        }
    }
    return &localTrans_;
}


void RegionManager::setZWeight(double z){
    for (std::map< SIndex, Region * >::const_iterator it  = regionMap_.begin();
         it != regionMap_.end(); it ++){
        it->second->setZWeight(z);
    }
    interRegionConstraintsZWeight_ = z;
}

void RegionManager::setConstraintType(Index type){
    for (std::map< SIndex, Region * >::const_iterator it  = regionMap_.begin();
                                                   it != regionMap_.end(); it ++){
        it->second->setConstraintType(type);
    }
}

} // namespace GIMLI
