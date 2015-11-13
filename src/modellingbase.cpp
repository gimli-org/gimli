/***************************************************************************
 *   Copyright (C) 2005-2015 by the resistivity.net development team       *
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
#include "sparsematrix.h"

namespace GIMLI{

ModellingBase::ModellingBase(bool verbose)
    : dataContainer_(NULL), verbose_(verbose){
    init_();
}

ModellingBase::ModellingBase(DataContainer & data, bool verbose)
    : dataContainer_(NULL), verbose_(verbose){
    init_();
    setData(data);
}

ModellingBase::ModellingBase(const Mesh & mesh, bool verbose)
    : dataContainer_(NULL), verbose_(verbose){
    init_();
    setMesh(mesh);
}

ModellingBase::ModellingBase(const Mesh & mesh, DataContainer & data, bool verbose)
    : dataContainer_(NULL), verbose_(verbose){
    init_();
    setData(data);
    setMesh(mesh);
}

ModellingBase::~ModellingBase() {
    delete regionManager_;
    if (mesh_) delete mesh_;
    if (jacobian_ && ownJacobian_) delete jacobian_;
    if (constraints_ && ownConstraints_) delete constraints_;
}

void ModellingBase::init_() {
    regionManager_      = new RegionManager(verbose_);
    regionManagerInUse_ = false;
    
    mesh_               = 0;
    jacobian_           = 0;
    constraints_        = 0;
    dataContainer_      = 0;
    
    nThreads_           = numberOfCPU();
    
    ownJacobian_        = false;
    ownConstraints_     = false;
    
    initJacobian();
    initConstraints();
}
  
void ModellingBase::setThreadCount(Index nThreads) { 
    nThreads_=max(1, (int)nThreads); 
    GIMLI::setThreadCount(nThreads);
}
  
void ModellingBase::setData(DataContainer & data){
    //if (dataContainer_) {
    dataContainer_ = &data;
//     } else {
//         dataContainer_ = new DataContainer(data);
//     }
    updateDataDependency_();
}

DataContainer & ModellingBase::data() const{
    if (dataContainer_ == 0){
        throwError(1, WHERE_AM_I + " no data defined");
    }
    return * dataContainer_;
}
    
RVector ModellingBase::startModel() {
    //*! Create startmodel from default builder (abstract may be overloaded)
    //*! Create startmodel from regionManger
    if (startModel_.size() == 0 && regionManager_){
        setStartModel(regionManager_->createStartVector());
    }
    
    if (startModel_.size() == 0){
        setStartModel(createDefaultStartModel());
    }

    if (startModel_.size() == 0){
        std::cout << "Warning! there is no startmodel defined." << std::endl;
    }
    return startModel_;
}

void ModellingBase::setStartModel(const RVector & startModel){
    startModel_=startModel; 
}
        
void ModellingBase::createRefinedForwardMesh(bool refine, bool pRefine){
	this->initRegionManager();

    if (!mesh_) {
        mesh_ = new Mesh();
    } else {
        mesh_->clear();
    }
    
    if (refine){
        if (pRefine){
            *mesh_ = regionManager_->mesh().createP2();
        } else {
            *mesh_ = regionManager_->mesh().createH2();
        }
        updateMeshDependency_();
    } else {
        setMesh_(regionManager_->mesh());
    }
}
 
void ModellingBase::setRefinedMesh(const Mesh & mesh){
    if (verbose_) {
        std::cout << "set extenal secondary mesh:" << std::endl;
    }
    setMesh_(mesh);
    if (verbose_) {
        std::cout << "nModel = " << regionManager_->parameterCount() << std::endl;
        IVector m(unique(sort(mesh_->cellMarker())));
        std::cout << "secMesh marker = [" << m[0] <<", " << m[1] << ", " << m[2]  
         << ", ... ,  " << m[-1] << "]" << std::endl;
    }
}
 
void ModellingBase::setMesh(const Mesh & mesh, bool holdRegionInfos) {
    deleteMeshDependency_();

    Stopwatch swatch(true);
    if (regionManagerInUse_){
        regionManager_->setMesh(mesh, holdRegionInfos);
        if (verbose_) std::cout << "ModellingBase::setMesh() switch to regionmanager mesh" << std::endl;
        this->setMesh_(regionManager_->mesh());
    } else {
        if (verbose_) std::cout << "ModellingBase::setMesh() copying new mesh ... ";
        this->setMesh_(mesh);
        if (verbose_) std::cout << swatch.duration(true) << " s" << std::endl;
    }

    if (verbose_) std::cout << "FOP updating mesh dependencies ... ";
    startModel_.clear();
    updateMeshDependency_();
    
    if (verbose_) std::cout << swatch.duration(true) << " s" << std::endl;
}

void ModellingBase::setMesh_(const Mesh & mesh){
    if (!mesh_) mesh_ = new Mesh();
    (*mesh_) = mesh;
}

void ModellingBase::initJacobian(){
    if (!jacobian_){
        jacobian_ = new RMatrix();
        ownJacobian_ = true;
    }
}

void ModellingBase::setJacobian(MatrixBase * J){
    if (!jacobian_ && ownJacobian_){
        delete jacobian_;
    }
    jacobian_ = J;
    ownJacobian_ = false;
}

void ModellingBase::createJacobian(const RVector & model){
    if (verbose_) std::cout << "Create Jacobian matrix (brute force) ...";

    Stopwatch swatch(true);
    RVector resp(response(model));

    if (!jacobian_){
        this->initJacobian();
    }
 
    RMatrix *J = dynamic_cast< RMatrix * >(jacobian_);
    
    if (J->rows() != resp.size()){ J->resize(resp.size(), model.size()); }

    double fak = 1.05;
    for (size_t i = 0; i < model.size(); i++) {
        RVector modelChange(model);
        modelChange[i] *= fak;
        
        RVector respChange(response(modelChange));

        for (size_t j = 0; j < resp.size(); j++){
            if (::fabs(modelChange[i] - model[i]) > TOLERANCE){
                (*J)[j][i] = (respChange[j] - resp[j]) / (modelChange[i] - model[i]);
            } else {
                (*J)[j][i] = 0.0;
            }
        }
    }

    swatch.stop();
    if (verbose_) std::cout << " ... " << swatch.duration() << " s." << std::endl;
}
    
void ModellingBase::initConstraints(){
    if (constraints_ == 0){
        constraints_ = new RSparseMapMatrix(0, 0, 0);
        ownConstraints_ = true;
    } 
}

void ModellingBase::setConstraints(MatrixBase * C){
    if (!constraints_ && ownConstraints_){
        delete constraints_;
    }
    constraints_ = C;
    ownConstraints_ = false;
}

void ModellingBase::createConstraints(){
//     __MS(constraints_->rtti())
    this->regionManager().fillConstraints(constraintsRef());
}
    
void ModellingBase::clearConstraints(){
    if (constraints_) constraints_->clear();
}
    
MatrixBase * ModellingBase::constraints() { 
    return constraints_;
}
    
MatrixBase * ModellingBase::constraints() const { 
    return constraints_;
}
    
RSparseMapMatrix & ModellingBase::constraintsRef() const { 
    if (!constraints_) throwError(1, WHERE_AM_I + " constraints matrix is not initialized.");
    return *dynamic_cast < RSparseMapMatrix *>(constraints_); 
}
    
RSparseMapMatrix & ModellingBase::constraintsRef() { 
    if (!constraints_) throwError(1, WHERE_AM_I + " constraints matrix is not initialized.");
    return *dynamic_cast < RSparseMapMatrix *>(constraints_); 
}
        
    
void ModellingBase::mapModel(const RVector & model, double background){
    int marker = -1;
    std::vector< Cell * > emptyList;
    mesh_->createNeighbourInfos();

    for (Index i = 0, imax = mesh_->cellCount(); i < imax; i ++){
        marker = mesh_->cell(i).marker();
        if (marker >= 0) {
            if ((size_t)marker >= model.size()){
                mesh_->exportVTK("mapModelfail");
                std::cerr << WHERE_AM_I << std::endl
                          << "Wrong mesh here .. see mapModelfail.vtk" << std::endl
                          << *mesh_ << std::endl
                          << "mesh contains " << unique(sort(mesh_->cellMarker())).size() << " unique marker. " << std::endl;
                throwLengthError(1, WHERE_AM_I + " marker greater = then model.size() " + toStr(marker)
                       + " >= " + toStr(model.size()));
            }
            if (model[marker] < TOLERANCE){
                emptyList.push_back(&mesh_->cell(i));
            }
            mesh_->cell(i).setAttribute(model[marker]);

        } else {
            // general background without fixed values, fixed values will be set at the end
            if (marker == -1) { 
                mesh_->cell(i).setAttribute(0.0);
                emptyList.push_back(&mesh_->cell(i));
            }
        }
    }

    // if background == 0.0 .. empty cells are allowed
    if (emptyList.size() == mesh_->cellCount() && background != 0.0){
        throwLengthError(1, WHERE_AM_I + " too many empty cells" + toStr(emptyList.size())
                       + " == " + toStr(mesh_->cellCount()));
    }

    if (background != 0.0){
        mesh_->fillEmptyCells(emptyList, background);
    }
    
    // setting fixed values
    if (regionManagerInUse_){
        for (Index i = 0, imax = mesh_->cellCount(); i < imax; i ++){
            if (abs(mesh_->cell(i).attribute()) < TOLERANCE){
                if (mesh_->cell(i).marker() <= MARKER_FIXEDVALUE_REGION){
                    SIndex regionMarker = -(mesh_->cell(i).marker() - MARKER_FIXEDVALUE_REGION);
                    double val = regionManager_->region(regionMarker)->fixValue();
                    __DS("fixing region: " << regionMarker << " to: " << val)
                    mesh_->cell(i).setAttribute(val);
                }
            }
        }
    }
}

   
void ModellingBase::initRegionManager() {
    if (!regionManagerInUse_){
        if (mesh_){
            regionManager_->setMesh(*mesh_);
            this->setMesh_(regionManager_->mesh());
        }
        regionManagerInUse_   = true;
    }
}

const RegionManager & ModellingBase::regionManager() const {
    if (regionManager_ == 0) throwError(1, "No RegionManager initialized");
    return *regionManager_;
}

RegionManager & ModellingBase::regionManager(){
    initRegionManager();
    return *regionManager_;
}

Region * ModellingBase::region(int marker){
    return regionManager().region(marker);
}

RVector ModellingBase::createStartVector() {
    return regionManager().createStartVector();
}


LinearModelling::LinearModelling(MatrixBase & A, bool verbose)
    : ModellingBase(verbose){//, A_(& A) {
        setJacobian(&A);
        this->regionManager().setParameterCount(A.cols());
}

RVector LinearModelling::response(const RVector & model) {
    if (jacobian_->cols() != model.size()){
        throwLengthError(1, WHERE_AM_I + " Jacobian col size != model.size()"
                                + toStr(jacobian_->cols()) + " != " + toStr(model.size()));
    }
    return *jacobian_ * model;
}

RVector LinearModelling::createDefaultStartModel() {
    return RVector(jacobian_->cols(), 1.0);
}

} // namespace GIMLI

