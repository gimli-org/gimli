/***************************************************************************
 *   Copyright (C) 2006-2015 by the resistivity.net development team       *
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

#ifndef _GIMLI_INVERSION__H
#define _GIMLI_INVERSION__H

#include "vector.h"
#include "inversionBase.h"
#include "mesh.h"
#include "modellingbase.h"
#include "numericbase.h"
#include "regionManager.h"
#include "solver.h"
#include "stopwatch.h"
#include "sparsematrix.h"
#include "trans.h"
#include "vector.h"

#include "ipcClient.h"

namespace GIMLI{

#define DOSAVE if (dosave_)

//#define PLUS_TMP_VECSUFFIX + ".vec"
#define PLUS_TMP_VECSUFFIX

template < class Vec > Vec getIRLSWeights(const Vec & a, double locut = 0.0, double hicut = 0.0) {
    double suabs = sum(abs(a));
    double suabsq = dot(a, a);

    Vec tmp(suabsq / suabs / (abs(a) + TOLERANCE));
    for(uint i = 0; i < a.size(); i++) {
        if((locut > 0.0) && (tmp[ i ] < locut)) tmp[ i ] = locut;
        if((hicut > 0.0) && (tmp[ i ] > hicut)) tmp[ i ] = hicut;
    }
    return tmp;
}

template < class Vec > Vec getIRLSWeightsP(const Vec & a, int p, double locut = 0.0, double hicut = 0.0) {
    Vec ap = pow(abs(a), p);
    Vec aq = pow(abs(a), 2);
    double suabsp = sum(ap);
    double suabsq = sum(aq);

    Vec tmp((ap + 1e-12) / (aq + 1e-12) * suabsp / suabsq);
    for(uint i = 0; i < a.size(); i++) {
        if((locut > 0.0) && (tmp[ i ] < locut)) tmp[ i ] = locut;
        if((hicut > 0.0) && (tmp[ i ] > hicut)) tmp[ i ] = hicut;
    }
    return tmp;
}

/*! Inversion template using a Generalized Minimization Approach, atm fixed to Gauss-Newton solver
    Inversion(bool verbose, bool dosave
    Inversion(RVector data, FOP f, bool verbose, bool dosave
    Inversion(RVector data, FOP f, transData, transModel, bool verbose, bool dosave */
template < class ModelValType > class Inversion : public InversionBase< ModelValType > {
public:
    typedef Vector < ModelValType > Vec;
    typedef Vector < ModelValType > ModelVector;

    /*! Minimal constructor. verbose -- gives some output and status information about the inversion progress; dosave -- advanced debuging, saves alot of stuff and temporary data */
    Inversion(bool verbose=false, bool dosave=false)
    : InversionBase< ModelValType >(), verbose_(verbose), 
        dosave_(dosave), saveModelHistory_(dosave){
        this->init_();
    }

    /*! Usual constructor. data -- vector of the given data; forward -- the associated forward operator */
    Inversion(const Vec & data, ModellingBase & forward,
              bool verbose=false, bool dosave=false)
    : InversionBase< ModelValType >(), verbose_(verbose),
        dosave_(dosave), saveModelHistory_(dosave) {
        this->init_();

        //** init: data_
        this->setData(data);

        //** init: model_, modelRef, response_
        this->setForwardOperator(forward);
    }

    /*! Intermediate constructor including transforms and data transformation function */
    Inversion(const Vec & data, ModellingBase & forward,
              Trans< Vec > & transData, bool verbose=true, bool dosave=false)
    : InversionBase< ModelValType >(), verbose_(verbose),
        dosave_(dosave), saveModelHistory_(dosave) {
        //** set: default values
        this->init_();
        //** set: paraDomain init: modelWeight_, constraintsWeight_, ConstraintMatrix
        //** init: data_
        this->setData(data);
        //** init: model_, modelRef, response_
        this->setForwardOperator(forward);

        this->setTransData(transData);
    }

    /*! Full constructor including transforms. transData -- data transformation function; transModel -- model transformation function */
    Inversion(const Vec & data, ModellingBase & forward,
              Trans< Vec > & transData, Trans< Vec > & transModel,
              bool verbose = true, bool dosave = false)
    : InversionBase< ModelValType >(), verbose_(verbose), dosave_(dosave), saveModelHistory_(dosave) {
        //** set: default values
        this->init_();
        //** set: paraDomain init: modelWeight_, constraintsWeight_, ConstraintMatrix
        //** init: data_
        this->setData(data);
        //** init: model_, modelRef, response_
        this->setForwardOperator(forward);

        this->setTransData(transData);
        this->setTransModel(transModel);
    }

    /*! Destructor. Frees allocated memory */
    virtual ~Inversion(){
        delete transDataDefault_;
        delete transModelDefault_;
    }

private:
    /*! Copyconstructor */
    Inversion(const Inversion< ModelValType > & inv){
        THROW_TO_IMPL
    }


protected:
    /*! Internal initialization function, which is called from constructor. Set default paramater and allocate required memory */
    void init_(){
        transDataDefault_   = new Trans< Vec >;
        transModelDefault_  = new Trans< Vec >;
        tD_                 = transDataDefault_;
        tM_                 = transModelDefault_;
        isRobust_           = false;
        isBlocky_           = false;
        useLinesearch_      = true;
        optimizeLambda_     = false;
        recalcJacobian_     = true;
        doBroydenUpdate_    = false;
        localRegularization_= false;
        abort_              = false;
        stopAtChi1_         = true;
        haveReferenceModel_ = false;
        isRunning_          = false;
        fixError_           = true;
        activateFillConstraintsWeight_ = true; //jointinv hack!!!

        iter_               = 0;
        maxiter_            = 20;
        maxCGLSIter_        = 200;
        lambda_             = 50.0;
        lambdaFactor_       = 1.0;
        dPhiAbortPercent_   = 2.0;
        phiD_               = 0.0;

        CGLStol_            = -1.0; //** -1 means automatic scaled
    }

public:

    /*! Set data vector */
    inline void setData(const Vec & data) {
        data_ = data;
        //** maybe a data validation here
    }

    /*! Get data vector */
    inline const Vec & data() const { return data_; }

    /*! Set the relative data error to the vector err. */
    void setRelativeError(const Vec & e){
        error_ = e;
        checkError();
    }

    /*! Set relative data error to scalar error value relerr.
     * If you force relerr == 0 here. Data will not corrected or weighted. */
    inline void setRelativeError(double relerr) {
        if (relerr == 0.0) fixError_ = false;
        setRelativeError(RVector(data_.size(), relerr));
    }

    /*! Set the absolute data error to the vector abserr */
    inline void setAbsoluteError(const Vec & abserr) {
        setRelativeError(abs(abserr) / abs(fixZero(data_, TOLERANCE)));
    }

    /*! Set absolute data error to the scalar value abserr */
    inline void setAbsoluteError(double abserr) {
        setRelativeError(std::fabs(abserr) / abs(fixZero(data_, TOLERANCE)));
    }

    /*! Old vector error setter still fixed to relative error, should not be used due to ambiguity */
    inline void setError(const Vec & err) { setRelativeError(err); }

    /*! Old scalar error setter still fixed to relative error, should not be used due to ambiguity */
    inline void setError(double err) { setRelativeError(err); }

    /*! Return the used data error */
    inline const Vec & error() const { return error_; }

    /*! Validate the data error. This check will be called at every new run or change of the error;
        Check size and content unequal to zero. If problems found, warn and set the data error to errorDefault_()  */
    void checkError(){
        if (error_.size() != data_.size()) {
            std::cerr << WHERE_AM_I << " Warning error has the wrong size, reset to default. "
                      << error_.size() << " != " << data_.size() << std::endl;
            error_ = errorDefault_();
        } else if (min(abs(error_)) < TOLERANCE && fixError_) {
            std::cerr << WHERE_AM_I << " Warning error contains zero values, reset to default. " << std::endl;
            error_ = errorDefault_();
        }


        dataWeight_ = 1.0 / tD_->error(fixZero(data_, TOLERANCE), error_);

        if (verbose_) std::cout << "min/max(dweight) = " << min(dataWeight_) << "/"
                                  << max(dataWeight_) << std::endl;

        if (haveInfNaN(dataWeight_)){
            DOSAVE save(dataWeight_,          "Nan_dataWeight_dweight");
            DOSAVE save(error_,               "Nan_dataWeight_error");
            DOSAVE save(data_,                 "Nan_dataWeight_data");

            throwError(1, WHERE_AM_I + " dataWeight_ contains inf or nan" );
        }

    }

    /*! Return a copy of the default data error array, [ Relative data error of 1% ] */
    inline Vec errorDefault_(){ return Vec(data_.size(), 0.01); }

    /*! Set the forward operator and the model-transform-function derived from fop-regionManager if not set before. */
    inline void setForwardOperator(ModellingBase & forward) {
        forward_   = & forward;
        forward_->clearConstraints();

        //! Always use a region manager
        forward_->initRegionManager();
//         model_ = forward_->startModel();
    }

    /*! Return forward operator.*/
    inline ModellingBase * forwardOperator() { return forward_; }

    /*TODO Change interface to Trans < vec > *tD, or copy the trans ,giving non const reference is dangerous.*/
    /*! Set and get data transform. 
     * WARNING! we just store a reference to the trans .. U have to ensure that the trans is and stayes valid. 
     * TODO change this!!!
     */
    inline void setTransData(Trans< Vec > & tD) { tD_ = & tD; }
    inline Trans< Vec > & transData() { return * tD_; }

    /*TODO Change interface to Trans < vec > *tD, or copy the trans ,giving non const reference is dangerous. */
    /*! Set and get Model transform 
     * WARNING! we just store a reference to the trans .. U have to ensure that the trans is and stayes valid. 
     * TODO change this!!!
     */
    inline void setTransModel(Trans< Vec > & tM) { tM_ = & tM; }
    inline Trans< Vec > & transModel() { return * tM_; }

    /*! Return number of contraints and data */
    inline uint constraintsCount() const { return constraintsWeight_.size(); }
    inline uint dataCount()        const { return data_.size(); }

    /*! Set and get verbose behaviour */
    inline void setVerbose(bool verbose){ verbose_ = verbose; }
    inline bool verbose() const { return verbose_; }

    /*! Set and get verbose behaviour */
    inline void saveModelHistory(bool doSaveModelHistory){ saveModelHistory_ = doSaveModelHistory; }

    /*! Set and get line search */
    inline void setLineSearch(bool linesearch) { useLinesearch_ = linesearch; }
    inline bool lineSearch() const { return useLinesearch_; }

    /*! Set and get blocky model behaviour (by L1 reweighting of constraints) */
    inline void setBlockyModel(bool isBlocky) { isBlocky_ = isBlocky; }
    inline bool blockyModel() const { return isBlocky_; }

    /*! Set and get robust data weighting by L1 reweighting scheme */
    inline void setRobustData(bool isRobust) { isRobust_ = isRobust; }
    inline bool robustData() const { return isRobust_; }

    /*! Set and get regularization parameter and its change over iteration */
    inline void setLambda(double lambda) { lambda_ = lambda; }
    inline double lambda() const { return lambda_; }
    inline void setLambdaFactor(double lambdaFactor) { lambdaFactor_ = lambdaFactor; }
    inline double lambdaFactor() const { return lambdaFactor_; }

    /*! Set whether regularization is global or local (e.g. Marquardt method) */
    void setLocalRegularization(bool localReg){ localRegularization_ = localReg; }
    /*! Return whether regularization is global or local (e.g. Marquardt method) */
    bool localRegularization() const { return localRegularization_; }

    /*! Set and get whether lambda is being optimized by Lcurve */
    inline void setOptimizeLambda(bool opt) { optimizeLambda_ = opt; }
    inline bool optimizeLambda() const { return optimizeLambda_; }

    /*! Set and get maximum iteration number */
    inline void setMaxIter(int maxiter) { maxiter_ = maxiter; }
    inline int maxIter() const { return maxiter_; }

    /*! Set and get maximum iteration number for inverse sub step*/
    inline void setMaxCGLSIter(int iter){ maxCGLSIter_ = iter; }
    inline int maxCGLSIter() const { return maxCGLSIter_; }

    /*! Set and get abort tolerance for cgls solver. -1 means default scaling [default] */
    inline void setCGLSTolerance(double tol){ CGLStol_ = tol; }
    inline double maxCGLSTolerance() const { return CGLStol_; }

    /*! Return curent iteration number */
    inline uint iter() const { return iter_; }

    /*! Return true if inversion is running */
    inline bool isRunning() const { return isRunning_ ; }

    /*! Force abort at next iteration */
    inline void abort() { abort_ = true; }

    /*! Define whether minimization stops at chi^2=1 */
    inline void stopAtChi1(bool stopAtChi1) { stopAtChi1_ = stopAtChi1; }

    /*! Set and get relative objective function decrease */
    inline void setDeltaPhiAbortPercent(double dPhi){ dPhiAbortPercent_ = dPhi; }
    inline double deltaPhiAbortPercent() const { return dPhiAbortPercent_; }

    /*! Marquardt scheme (local damping with decreasing regularization strength */
    void setMarquardtScheme(double lambdaFactor=0.8){
        setLocalRegularization(true);  //! no contribution of regularization to objective function
        setLambdaFactor(lambdaFactor); //! lambda is decreased towards zero
        stopAtChi1(false);             //! let the solution fully converge (important for statistics)
        forward_->regionManager().setConstraintType(0);
    }

    /*! Check if transfunctions are valid, set default if no transfunctions are given or override given functions by regionManager.transfunctions if available*/
    void checkTransFunctions(){
        if (forward_->regionManager().haveLocalTrans()) {
            if (verbose_) std::cout << "use cumulative model trans " << std::endl;
            tM_ = forward_->regionManager().transModel();
        }
    }

    /*! Create constraints, check and compare size of constraint matrix with model/boundary control */
    void checkConstraints() {
//         __MS(forward_->constraints()->rtti())
        
        if (forward_->constraints()->cols() == 0 ||
            forward_->constraints()->rows() == 0){
            if (verbose_) std::cout << "Building constraints matrix" << std::endl;
            //forward_->regionManager().fillConstraints(forward_->constraints());
            forward_->createConstraints();
        } else {
            if (verbose_) std::cout << " found valid constraints matrix. omit rebuild" << std::endl;
        }
        Index nModelC = forward_->constraints()->cols();
        Index nCWeightC = forward_->constraints()->rows();

        if (verbose_){
            std::cout << "constraint matrix of size(nBounds x nModel) " 
                << nCWeightC << " x " << nModelC << std::endl;
        }
        
        if (model_.size() != nModelC){
            //** vielleicht irgendwo ne schöne resize funktion
            std::cout << WHERE_AM_I << " resize model " << model_.size()
                      << " to fit constrain size: " 
                      << nModelC << std::endl;
            model_.resize(nModelC);
        }

        if (modelWeight_.size() != model_.size()){
            modelWeight_.resize(model_.size(), 1.0);
        }

        if (constraintsWeight_.size() != nCWeightC) {
            constraintsWeight_.resize(nCWeightC, 1.0);
        }
    }

    /*! Check size of Jacobian matrix against data and model number */
    void checkJacobian() {
        if (forward_->jacobian()->rows() == data_.size() &&
             forward_->jacobian()->cols() == model_.size()) return;

        if (verbose_ && forward_->jacobian()->rows() + forward_->jacobian()->cols() > 0){
            std::cout << "check Jacobian: wrong dimensions: "
                        << "(" << forward_->jacobian()->rows()  << "x" << forward_->jacobian()->cols() << ") == "
                        << "(" << data_.size() << "x" << model_.size()  << ") "  << std::endl;
            std::cout << "jacobian size invalid, forced recalc" << std::endl;
        }
        Stopwatch swatch(true);
        if (verbose_) std::cout << "calculating jacobian matrix ...";
        forward_->createJacobian(model_);
        if (verbose_) std::cout << "... " << swatch.duration(true) << " s" << std::endl;
    }

    /*! Define and find out whether Jacobian is recalculated in each iteration */
    void setRecalcJacobian(bool recalc){
        recalcJacobian_ = recalc;
        if (recalc) doBroydenUpdate_ = false;
    }
    bool recalcJacobian() const { return recalcJacobian_; }

    /*! Enable/disable broyden update (which enforces scaling by transform function derivatives */
    void setBroydenUpdate(bool broydenUpdate){
        doBroydenUpdate_ = broydenUpdate;
        if (doBroydenUpdate_) recalcJacobian_ = false;
        //** puts transformation functions into jacobian since it is not changed anymore;
        if (broydenUpdate) {
            THROW_TO_IMPL
            //scaleMatrix(*J_, tD_->deriv(response_), Vec(1.0 / tM_->deriv(model_)));
        }
    }
    inline bool broydenUpdate() const { return doBroydenUpdate_; }

    /*! Set model vector. */
    void setModel(const Vec & model){ model_ = model; }

    /*! Return a const reference to the current model vector */
    inline const ModelVector & model() const { return model_; }

    /*! Set reference model vector */
    inline void setReferenceModel(const Vec & model){
        haveReferenceModel_ = true; modelRef_ = model;
    }

    /*! Set current model response (e.g. in order to avoid time-consuming forward calculation */
    inline void setResponse(const Vec & response){ response_ = response; }

    /*! Return a reference to the current response vector */
    inline const Vec & response() const { return response_; }

    /*! Return IRLS function of roughness vector */
    Vec getIRLS() const {
        return getIRLSWeights(Vec(*forward_->constraints() * tM_->trans(model_) * constraintsWeight_), 0.0, 1.0);
    }

    /*! Set the constraint weight (boundary control) vector */
    void setCWeight(const Vec & cweight){
        constraintsWeight_ = cweight;
        activateFillConstraintsWeight_ = false; //jointinv hack!!!
        if (verbose_) std::cout << "min/max(cWeight) = " << min(constraintsWeight_) << "/" << max(constraintsWeight_) << std::endl;
    }

    /*! Return reference to the current constraints weight vector */
    inline const Vec & cWeight() const { return constraintsWeight_; }

    /*! Set the model weight vector */
    void setMWeight(const Vec & mweight){
        modelWeight_ = mweight;
        if (verbose_) std::cout << "min/max(mWeight) = " << min(modelWeight_) << "/" << max(modelWeight_) << std::endl;
    }

    /*! Return reference to the current model weight vector */
    inline const Vec & mWeight() const { return modelWeight_; }

    /*! Call robust date reweighting (see also setRobustData) */
    void robustWeighting() {
        if (verbose_) std::cout << "Robust reweighting " << std::endl;
        Vec deltaData((tD_->trans(data_) - tD_->trans(response_)) * dataWeight_);
        //** reweight errors according to IRLS scheme
        // cr: error wird nur noch von getPhiD genutzt und dataweight wird im Laufe der Inversion nicht mehr aktualisiert, soll das so sein? DRY?
        // tg: checkerror ist doch genau das richtige oder?
        error_ /= (getIRLSWeights(deltaData, 0.0, 1.0) + TOLERANCE);
        checkError(); //** macht implizit:
        // dataWeight_ = (1.0 / tD_->error(data_, error_));
        // tg: maybe better just reweight dweight and hold error constant
        // tg: geht nicht weil phiD auf Fehlern beruht (und das auch muss wegen Trafowechsel)
    }

    /*! Apply blocky model constraints (see also setBlockyModel) */
    void constrainBlocky() {
        if (verbose_) std::cout << "Blocky model constraints " << std::endl;
        setCWeight(getIRLSWeights(Vec((*forward_->constraints() * tM_->trans(model_)) * constraintsWeight_), 0.0, 1.0));
    }

    /*! Shortcut for \ref echoStatus(response_, model_). */
    inline void echoStatus() const { echoStatus(response_, model_); }

    /*! Echo PhiD/PhiM/Phi and Chi^2 for given model and response */
    void echoStatus(const Vec & response, const Vec & model, const std::string & xtra = "") const {
        double chi2 = getPhiD(response) / data_.size();
        std::cout << iter_ << ": " << xtra
                  << "Model: min = " << min(model) << "; max = " << max(model) << std::endl;
        std::cout << iter_ << ": " << xtra
                  << "Response: min = " << min(response) << "; max = " << max(response) << std::endl;
        std::cout << iter_ << ": rms/rrms(data, " << xtra << "Response) = "
                        << rms(data_, response) << "/"
                        << rrms(data_, response) * 100.0 << "%" << std::endl;
        std::cout << iter_ << ": chi^2(data, " << xtra
                  << "Response, error, log) = " << chi2 <<        std::endl;
        std::cout << iter_ << ": Phi = " << getPhiD(response) << "+" << getPhiM(model)
                            << "*" << lambda_ << "=" << getPhi(model, response) << std::endl;
    }

    /*! Compute model cell resolution (specific column of resolution matrix) by an LSCG solver (Guenther, 2004) */
    Vec modelCellResolution(int iModel) {
        Vec resolution(model_.size(), 0.0);
        resolution[ iModel ] = 1.0; //** both for retrieving solumn and as starting model for resolution
        //** retrieve one colomn from the jacobian matrix and scale according to data/model transformation
        Vec sensCol = (*forward_->jacobian()) * resolution * tD_->deriv(response_) / tM_->deriv(model_)[ iModel ];
        //** call inverse substep with sensCol on right hand side (see Gnther, 2004)
        Vec deltaModel0(model_.size());// !!! h   zvariante
        
        solveCGLSCDWWtrans(*forward_->jacobian(), *forward_->constraints(),
                           dataWeight_, sensCol, resolution,
                           constraintsWeight_, modelWeight_,
                           tM_->deriv(model_), tD_->deriv(response_),
                           lambda_, deltaModel0, maxCGLSIter_, false);

        return resolution;
    }

    /*! Compute cell resolution for closest cell to given position */
    Vec modelCellResolution(const RVector3 & pos){
        Index ipos = 0;
        double mindist = 9e99;
        R3Vector cellCenters = forward_->mesh()->cellCenter();
        
        for (Index i = 0 ; i < model_.size() ; i++){
            double dist = cellCenters[i].distance(pos);
            if (dist < mindist) {
                mindist = dist;
                ipos = i;
            }
        }
        // Achtung Parametermapping und Backgroundzellen!
        return modelCellResolution(ipos);
    }

    /*! Compute whole resolution matrix by single cell resolutions */
    RMatrix modelResolutionMatrix() {
        RMatrix resMat;
        for (size_t i = 0 ; i < model_.size(); i++) resMat.push_back(modelCellResolution(i));
        return resMat;
    }

    /*! Return (C * m * m_w) * c_w */
    RVector roughness(const RVector & model) const {
       return *forward_->constraints() * Vec(tM_->trans(model) * modelWeight_) * constraintsWeight_;
    }

    /*! Shortcut for roughness for the current model vector */
    RVector roughness() const {
       return roughness(model_);
    }

    /*! Return data objective function (sum of squared data-weighted misfit) */
    double getPhiD(const Vec & response) const {
        Vec deltaData((tD_->trans(data_) - tD_->trans(response)) /
                       tD_->error(fixZero(data_, TOLERANCE), error_)) ;

        double ret = dot(deltaData, deltaData);
        if (isnan(ret) || isinf(ret)){
            save(tD_->trans(data_),          "Nan_PhiD_tD_data");
            save(response,                     "Nan_PhiD_response");
            save(tD_->trans(response),       "Nan_PhiD_tD_response");
            save(tD_->error(data_, error_), "Nan_PhiD_tD_error");

            throwError(1, WHERE_AM_I + " getPhiD == " + str(ret));
        }
        return ret;
    }

    /*! Return model objective function (squared model roughness) */
    double getPhiM(const Vec & model) const {
//        Vec dModel(tM_->trans(model));
//        if (haveReferenceModel_) dModel = dModel - tM_->trans(modelRef_);
//        Vec roughness(Vec(forward_->constraints() * dModel) * constraintsWeight_);
        Vec rough(this->roughness(model));
        if (haveReferenceModel_) rough = rough - constraintsH_;

        double ret = dot(rough, rough);
        if (isnan(ret) || isinf(ret)){
            std::cerr << "haveReferenceModel_: " << haveReferenceModel_<< std::endl;
            DOSAVE save(model,       "Nan_PhiM_model");
            DOSAVE save(modelRef_,  "Nan_PhiM_modelref");
//            DOSAVE save(dModel,      "Nan_PhiM_tM_dmodel");
            DOSAVE save(rough,   "Nan_PhiM_roughness");
            DOSAVE save(constraintsWeight_,   "Nan_PhiM_cweight");

            throwError(1, WHERE_AM_I + " getPhiM == " + str(ret));
        }
        return ret;
    }

    /*! Return total objective function (data OF plus lambda times model OF) */
    double getPhi(const Vec & model, const Vec & response) const {
        return getPhiD(response) + getPhiM(model) * lambda_ * (1.0 - double(localRegularization_));
    }

    /*! Shortcut for \ref getPhiD(response_) necessary ? */
    inline double getPhiD() const { return getPhiD(response_); }

    /*! Shortcut for \ref getPhiM(model_) necessary ? */
    inline double getPhiM() const { return getPhiM(model_); }



    /*! Shortcut for \ref getPhi(model_, response_) necessary ? */
    inline double getPhi() const { return getPhiD() + getPhiM() * lambda_ * (1.0 - double(localRegularization_)); }

    /*! Return the chi-squared data misfit */
    inline double getChi2() const { return getPhiD() / data_.size(); }

    /*! Return the chi-squared data misfit for given forward response */
    inline double getChi2(const Vec & response) const { return getPhiD(response) / data_.size(); }

    /*! Return last chi-squared data misfit */
    inline double chi2() const { return phiD_ / data_.size(); }

    /*! Return last absolute RMS misfit */
    inline double absrms() const { return rms(data_, response_); }

    /*! Return last relative RMS misfit */
    inline double relrms() const { return rrms(data_, response_) * 100.; }

    /*! Start with linear interpolation, followed by quadratic fit if linesearch parameter tau is lower than 0.03. Tries to return values between 0.03 and 1 */
    double linesearch(const Vec & modelNew, const Vec & responseNew) const {
        Vec phiVector(101, getPhi());
        Vec phiDVector(101 , getPhiD());

        Vec dModel(tM_->trans(modelNew)    - tM_->trans(model_));
        Vec dData( tD_->trans(responseNew) - tD_->trans(response_));

        double tau = 0.0, minTau = 0.0;
        double minPhi = phiVector[ 0 ];

        if (localRegularization_) minPhi = phiDVector[ 0 ];

        double thisPhi = minPhi;
        for (int i = 1; i < 101; i++) {
            tau = 0.01 * (double) i;
            Vec appModel (tM_->update(model_, dModel * tau));
            Vec appResponse (tD_->update(response_, dData * tau));
            phiVector[ i ]  = getPhi(appModel, appResponse);
            phiDVector[ i ] = getPhiD(appResponse);
            thisPhi = phiVector[ i ]; //! rigorous minimization of the total objective function
            //** this could also be controlled by another switch (which is enforced by local reg.)
            if (localRegularization_) thisPhi = phiDVector[ i ];
            if (thisPhi < minPhi){
                minPhi = thisPhi;
                minTau = tau;
            }
        }
        DOSAVE save(phiVector,  "linesearchPhi");
        DOSAVE save(phiDVector, "linesearchPhiD");

        tau = minTau;

        /*! parabolic line search using step 0, 1 and tauquad to fit parabola */
        if (tau < 0.03) {
            double tauquad = 0.3;
            if (verbose_) std::cout << "tau = " << tau
                            << ". Trying parabolic line search with step length " << tauquad;
            RVector modelQuad(tM_->update(model_, dModel * tauquad));
            RVector responseQuad (forward_->response(modelQuad));
            tau = linesearchQuad(modelNew, responseNew, modelQuad, responseQuad, tauquad);
            if (verbose_) std::cout << " ==> tau = " << tau;
            if (tau > 1.0) { //! too large
                tau = 1.0;
                if (verbose_) std::cout << " resetting to " << tau;
            }
            if (verbose_) std::cout << std::endl;

            if (tau < 0.03) { //! negative or nearly zero (stucked) -> use small step instead
                tau = 0.03;
                if (verbose_) std::cout << " tau < 0.03 ==> tau = " << tau << std::endl;
            }
        } // else tau > 0.03

        if (tau < 0.95) { //! only if not nearly one (which saves a forward run
            if (verbose_) echoStatus(responseNew, modelNew, "LS new");
        }

        if (verbose_) std::cout << "Linesearch tau = " << tau << std::endl;

        return tau;
    }

    /*! Compute objective function for old (tau=0), new (tau=1) and another model */
    double linesearchQuad(const Vec & modelNew, const Vec & responseNew,
                           const Vec & modelQuad, const Vec & responseQuad,
                           double tauquad) const {
        double phi0  = getPhi();
        double phi10 = getPhi(modelNew, responseNew)   - phi0;
        double phit0 = getPhi(modelQuad, responseQuad) - phi0;
        double dphit = phit0 - phi10 * tauquad;
        if (abs(dphit) < TOLERANCE) return 0.0;
        double tauopt = (phit0 - phi10 * tauquad * tauquad) / dphit / 2.0;

        DOSAVE std::cout << "LineSearchQuad: Phi = " << phi0 << " - " << phit0 + phi0
                         << " - " << phi10 + phi0 << " -> tau= " << tauopt << std::endl;
        return tauopt;
    }

    /*! Return the single models for each iteration. For debugging.*/
    inline const std::vector < RVector > & modelHistory() const { return modelHist_; }

    /*! Compute model update by solving one inverse sub-step
       \param rhs The right-hand-side vector of the system of equation
    */
    Vec invSubStep(const Vec & rhs) {
        Vec deltaModel0(model_.size());//!!! h-variante
        Vec solution(model_.size());
        solveCGLSCDWWtrans(*forward_->jacobian(), *forward_->constraints(), 
                           dataWeight_, rhs, solution, constraintsWeight_,
                           modelWeight_,
                           tM_->deriv(model_), tD_->deriv(response_),
                           lambda_, deltaModel0, maxCGLSIter_, verbose_);
        return solution;
    }

    /*! Optimization of regularization parameter by L-curve */
    Vec optLambda(const Vec & deltaData, const Vec & deltaModel0); //!!! h-variante

    /*! One iteration step. Return true if the step can be calculated successfully else false is returned. */
    bool oneStep();

    /*! Start the inversion procedure from staring model.*/
    const Vec & start();
    
    /*! Start the inversion procedure from the last model and return the final model vector.*/
    const Vec & run();

    /*! Specialized run function that tries to reach a datafit chi^2=1 by varying the regularization paramater lambda */
    Vec runChi1(double acc = 0.01, int maxiter = 50){
        stopAtChi1(false);
        Vec model = run();

        double lambda = lambda_;
        double chi2 = getChi2();
        double fak(2.0);
        double dir = 0, olddir = 0;
        bool verbose = verbose_;
//        setVerbose(false);
        forward_->setVerbose(false);
        if (verbose) std::cout << "Optimizing lambda subject to chi^2=1." << std::endl;
        if (verbose) std::cout << "chi^2 = " << chi2 << " lambda = " << lambda << " dir = " << dir << std::endl;
        int iter = 0;
        while (std::fabs(chi2 - 1.0) > acc and iter < maxiter){
            dir = - sign(chi2 - 1.0);                           //** direction: up (1) or down (-1)
            if (dir * olddir == -1) fak = std::pow(fak, 0.6); //** change direction: decrease step
            lambda *= std::pow(fak, dir);                       //** increase or decrease lambda
            setLambda(lambda);
            model = run();
            chi2 = getChi2();
            if(verbose) std::cout << "chi^2 = " << chi2 << " lambda = " << lambda << " dir = " << dir << std::endl;
            olddir = dir;                                         //** save old direction for step length
            iter++;
        }
        return model;
    }

    const RVector & modelWeight() const { return modelWeight_; }
    const RVector & modelRef() const { return modelRef_; }
    const RVector & dataWeight() const { return dataWeight_; }

    const RVector & deltaDataIter() const { return deltaDataIter_; }
    const RVector & deltaModelIter() const { return deltaModelIter_; }

    IPCClientSHM & ipc() { return ipc_; }
    
protected:

    Vec                   data_;
    ModellingBase       * forward_;

    Trans< Vec >        * tD_;
    Trans< Vec >        * tM_;
    Trans< Vec >        * transModelDefault_;
    Trans< Vec >        * transDataDefault_;

    bool verbose_;
    bool dosave_;
    bool saveModelHistory_;

    Vec                 error_;
    Vec                 response_;

    Vec                 model_;
    Vec                 modelRef_;

    Vec                 constraintsH_;
    Vec                 constraintsWeight_;
    Vec                 modelWeight_;
    Vec                 dataWeight_;

    Vec  deltaDataIter_;
    Vec  deltaModelIter_;

    int maxiter_;
    int iter_;
    int maxCGLSIter_;

    double lambda_;
    double lambdaFactor_;
    double dPhiAbortPercent_;
    double phiD_;
    double CGLStol_;

    bool isBlocky_;
    bool isRobust_;
    bool isRunning_;
    bool useLinesearch_;
    bool optimizeLambda_;
    bool abort_;
    bool stopAtChi1_;
    bool doBroydenUpdate_;
    bool localRegularization_;
    bool haveReferenceModel_;
    bool recalcJacobian_;
    bool activateFillConstraintsWeight_; //jointinv hack!!!

    /*! Set this to zero if u want to use absolute errors == zero*/
    bool fixError_;

    /*! Hold old models, for debuging */
    std::vector < RVector > modelHist_;

    IPCClientSHM ipc_;
};

/*! Start inversion from starting model. */
template < class ModelValType >
const Vector < ModelValType > & Inversion< ModelValType >::start(){ ALLOW_PYTHON_THREADS
    setModel(forward_->startModel());
    return run();
}
    
/*! Run inversion with current model. */
template < class ModelValType >
const Vector < ModelValType > & Inversion< ModelValType >::run(){ ALLOW_PYTHON_THREADS
    
    if (model_.size() == 0 ) setModel(forward_->startModel());
    
    if (data_.size() == 0 ) {
        throwError(1, WHERE_AM_I + " no data given");
    }
    
    abort_ = false;
    ipc_.setBool("abort", false);

    //** check if transfunctions are valid
    this->checkTransFunctions();

    //! calculation of initial modelresponse
    response_ = forward_->response(model_);
    //response_ = forward_->response(forward_->startModel());

    //! () clear the model history
    modelHist_.clear();

    //! validate and rebuild the data error if necessary
    this->checkError();

    //! validate and rebuild the constraint matrix if necessary
    this->checkConstraints();

    forward_->regionManager().fillModelControl(modelWeight_);

    if (activateFillConstraintsWeight_) {
        forward_->regionManager().fillConstraintsWeight(constraintsWeight_);
    }
    
    if (constraintsWeight_.size() != forward_->constraints()->rows()){
        std::cout << WHERE_AM_I << " Fixing cweight.size()" << std::endl;
        std::cout << constraintsWeight_.size() << " " 
                  << forward_->constraints()->rows() << std::endl;
        constraintsWeight_.resize(forward_->constraints()->rows(), 1.0);
    }
    if (modelWeight_.size() != forward_->constraints()->cols()){
        std::cout << WHERE_AM_I << " Fixing mweight.size()" << std::endl;
        std::cout << modelWeight_.size() << " " 
                  << forward_->constraints()->cols() << std::endl;
        modelWeight_.resize(forward_->constraints()->cols(), 1.0);
    }
    
    //! compute roughness constraint and correct it for inter-region constraints
    size_t cc = forward_->regionManager().constraintCount();
    
    if (constraintsH_.size() != cc) constraintsH_.resize(cc);
    
    if (haveReferenceModel_) {
        constraintsH_ = (*forward_->constraints() * Vec(tM_->trans(modelRef_) * modelWeight_)) * constraintsWeight_; //!!!template
        size_t ircc = forward_->regionManager().interRegionConstraintsCount();
        if (ircc > 0) constraintsH_.setVal(0.0, cc - ircc, long(cc));
    }

    //! validate and rebuild the jacobian if necessary
    this->checkJacobian();

    //** End preparation

    if (saveModelHistory_) { save(model_    , "model_0"   ); }
    DOSAVE save(response_ , "response_0");
    DOSAVE save(modelRef_ , "modelRef_0");
    DOSAVE save(RVector(response_ / data_ -1.0), "deltaData_0");
    DOSAVE forward_->constraints()->save("constraint.matrix");
    DOSAVE save(constraintsWeight_, "cweight_0");
    DOSAVE save(modelWeight_, "mweight_0");
    DOSAVE save(*forward_->jacobian(), "sens.bmat");

    DOSAVE std::cout << "C size: " << forward_->constraints()->cols() 
                     << " x " << forward_->constraints()->rows() << std::endl;

    phiD_ = getPhiD();

    if (verbose_) {
        echoMinMax(data_,  "data");
        echoMinMax(error_,  "error");
        echoMinMax(response_,  "response");
        if (haveReferenceModel_) {
            echoMinMax(modelRef_,  "reference model");
        } else {
            std::cout << "calc without reference model" << std::endl;
        }

        std::cout << 0 << ": rms/rrms(data, response) = " << rms(data_, response_)
                  << "/" << rrms(data_, response_) * 100.0 << "%" << std::endl;
        std::cout << 0 << ": chi^2(data, response, error, log) = " << phiD_ / data_.size() << std::endl;
        std::cout << 0 << ": Phi = " << phiD_ << " + " << getPhiM() << " * "
                                    << lambda_ << " = " << getPhi() << std::endl;
    }

    //** Start iteration
    iter_ = 0;
    double oldphi = phiD_;

    //** store initial model
    modelHist_.push_back(model_);

    isRunning_ = true;
    ipc_.setBool("running", true);

    while (iter_ < maxiter_ && !abort_){
        if (ipc_.getBool("abort")) break;

        if (!oneStep()) break;

        DOSAVE save(*forward_->jacobian() * model_, "dataJac_"  + toStr(iter_) PLUS_TMP_VECSUFFIX);
        DOSAVE save(response_,    "response_" + toStr(iter_) PLUS_TMP_VECSUFFIX);

        modelHist_.push_back(model_);

        phiD_ = getPhiD();

        if (stopAtChi1_ && (phiD_ < data_.size())) {
            if (verbose_) std::cout << "Reached data fit criteria (chi^2 <= 1). Stop." << std::endl;
            break;
        }

        if (stopAtChi1_ && (phiD_ < data_.size())) {
            if (verbose_) std::cout << "Reached data fit criteria (chi^2 <= 1). Stop." << std::endl;
            break;
        }

        double phi = getPhi();
        if (phi / oldphi > (1.0 - dPhiAbortPercent_ / 100.0) && iter_ > 2) {
            if (verbose_) std::cout << "Reached data fit criteria (delta phi < " << dPhiAbortPercent_
                        << "%). Stop." << std::endl;
            break;
        }

        oldphi = phi;

        if (isRobust_) robustWeighting();
        if (isBlocky_) constrainBlocky();
        if (lambdaFactor_ > 0.0) lambda_ *= lambdaFactor_;

    } //** while iteration;
    isRunning_ = false;
    ipc_.setBool("running", false);
    return model_;
} //** run

template < class Vec > bool Inversion< Vec>::oneStep() {
    iter_++;
    ipc_.setInt("Iter", iter_);

    if (verbose_) std::cout << "Iter: " << iter_ << std::endl;

    deltaModelIter_.resize(model_.size());
    deltaModelIter_ *= 0.0;
    deltaDataIter_ = (tD_->trans(data_) - tD_->trans(response_));

    if (sum(abs(deltaDataIter_)) < TOLERANCE) {
        if (verbose_) std::cout << "sum(abs(deltaDataIter_)) == Zero" << std::endl;
        return false;
    }
//    Vec deltaModel0(model_.size());
    Vec modelNew(   model_.size());
    Vec responseNew( data_.size());
    Vec roughness(constraintsH_.size(), 0.0);

    if (recalcJacobian_ && iter_ > 1) {
        Stopwatch swatch(true);
        if (verbose_) std::cout << "recalculating jacobian matrix ...";
        forward_->createJacobian(model_);
        if (verbose_) std::cout << swatch.duration(true) << " s" << std::endl;
    }

    if (!localRegularization_) {
        DOSAVE echoMinMax(model_, "model: ");

        roughness = this->roughness();

        if (haveReferenceModel_) {
            DOSAVE echoMinMax(modelRef_,  "reference model");
            roughness = roughness - constraintsH_;
        }
    } else {
        if (verbose_) std::cout << "use local regularization" << std::endl;
    }

    if (iter_ == 1 && optimizeLambda_) { //optimize regularization strength using L-curve
//        deltaModelIter_ = optLambda(deltaDataIter_, deltaModel0); //!!! h-variante

        /////////////*********************
        // fix this!!!!!!!!!!!!!1 constraintsH != deltaModel0
        /////////////*********************
        deltaModelIter_ = optLambda(deltaDataIter_, constraintsH_);
    } else {
        DOSAVE save(deltaDataIter_, "dd_" + toStr(iter_) PLUS_TMP_VECSUFFIX);
        DOSAVE echoMinMax(data_,      "data");
        DOSAVE echoMinMax(dataWeight_,  "dW");
        DOSAVE echoMinMax(deltaDataIter_,  "dd");
        DOSAVE echoMinMax(deltaModelIter_, "dm");
//         save(constraintsWeight_, "cw.tmp");
        DOSAVE echoMinMax(constraintsWeight_,  "cW");
        DOSAVE echoMinMax(modelWeight_,  "mW");
        DOSAVE echoMinMax(model_,    "mod");
        DOSAVE echoMinMax(response_, "resp");
//        DOSAVE echoMinMax(deltaModel0, "dM0");
        DOSAVE echoMinMax(constraintsH_, "constraintsH");
        DOSAVE save(constraintsH_, "constraintsH");
        DOSAVE save(tM_->deriv(model_), "modelTrans");
        DOSAVE save(tD_->deriv(response_), "responseTrans");

        if (doBroydenUpdate_) { //!!! h-variante
           if (verbose_) std::cout << "solve CGLSCDWW with lambda = " << lambda_ << std::endl;
                THROW_TO_IMPL
//                solveCGLSCDWW(*J_, forward_->constraints(), dataWeight_, deltaDataIter_, deltaModelIter_, constraintsWeight_,
//                                modelWeight_, lambda_, deltaModel0, maxCGLSIter_, verbose_);
        } else {
            if (verbose_) std::cout << "solve CGLSCDWWtrans with lambda = " << lambda_ << std::endl;
//             solveCGLSCDWWtrans(*J_, forward_->constraints(), dataWeight_, deltaDataIter_, deltaModelIter_, constraintsWeight_,
//                                  modelWeight_, tM_->deriv(model_), tD_->deriv(response_),
//                                lambda_, deltaModel0, maxCGLSIter_, verbose_);

            //save(forward_->jacobian(), "S"+ toStr(iter_) + ".mat", Ascii);

            // wannebee
//             DoubleWeightedMatrix scaledJacobian (forward_->jacobian(), tM_->deriv(model_), tD_->deriv(response_));
//             DoubleWeightedMatrix weightedConstraints(forward_->constraints(), constraintsWeight_, modelWeight_);
//             solveCGLSCDWWhtransWB(scaledJacobian, weightedConstraints, dataWeight_, deltaDataIter_, deltaModelIter_,
//                                    lambda_, roughness, maxCGLSIter_, verbose_);

            solveCGLSCDWWhtrans(*forward_->jacobian(), *forward_->constraints(),
                                dataWeight_, deltaDataIter_, deltaModelIter_,
                                constraintsWeight_, modelWeight_,
                                tM_->deriv(model_), tD_->deriv(response_),
                                lambda_, roughness, maxCGLSIter_, CGLStol_,
                                verbose_);
        } // else no broyden
    } // else no optimization

    DOSAVE echoMinMax(deltaModelIter_, "dm");

    modelNew = tM_->update(model_, deltaModelIter_);

    DOSAVE save(model_, "oldmodel");
    DOSAVE save(deltaModelIter_, "deltaModel");

    if (dosave_) {
        save(modelNew, "model_" + toStr(iter_) PLUS_TMP_VECSUFFIX);
    } else {
        if (saveModelHistory_) save(modelNew, "modelLS");
    }

    Vec responseLast(response_);
    responseNew = forward_->response(modelNew);

    double tau = 1.0;
    if (useLinesearch_){
        tau = linesearch(modelNew, responseNew);
    }

    if (tau >= 0.95){ //! full step possible;
        response_ = responseNew;
    } else { //! normal line search parameter between 0.03 and 0.94
        modelNew = tM_->update(model_, deltaModelIter_ * tau);
        response_ = forward_->response(modelNew);
    }

    model_ = modelNew;
    if (saveModelHistory_) save(model_, "model_" + toStr(iter_) PLUS_TMP_VECSUFFIX);

    if (verbose_) echoStatus();

    ipc_.setDouble("Chi2", getPhiD() / data_.size());

    if (doBroydenUpdate_) { //** perform Broyden update;
    // did not yet reflect moving jacobian into modellingbase
    THROW_TO_IMPL
//         if (verbose_) std::cout << "perform Broyden update" << std::endl;
//         Vec u(tD_->trans(response_) - tD_->trans(responseLast) - (*J_ * deltaModelIter_));
//         Vec v(deltaModelIter_ / dot(deltaModelIter_, deltaModelIter_));
//         rank1Update(*J_, u, v);
    }

    //!** temporary stuff
    if (forward_->mesh()){
// this forces pygimli/generatecode.py to create a ugly log10 declaration, which overwrites the valid log10 declarion
//        DOSAVE forward_->mesh()->addExportData("F-op-model(log10)", log10(forward_->mesh()->cellAttributes()));
        DOSAVE forward_->mesh()->addExportData("F-op-model", forward_->mesh()->cellAttributes());
        DOSAVE forward_->mesh()->exportVTK("fop-model" + toStr(iter_));
        DOSAVE forward_->mesh()->clearExportData();
    }

    return true;
}

template < class ModelValType >
Vector < ModelValType > Inversion< ModelValType >
    ::optLambda(const Vector < ModelValType > & deltaData,
                const Vector < ModelValType > & deltaModel0) {
        
ALLOW_PYTHON_THREADS

    std::vector< double > phiM, phiD;
    double ys = 0.0, yss = 0.0, curv = 0.0, oldcurv = -1e10;
    Vec oldDModel(model_.size());
    Vec uroldDModel(oldDModel);
    Vec deltaModel(model_.size());
    Vec tModel(tM_->trans(model_));
    Vec tResponse(tM_->trans(response_));
    Vec roughness(constraintsH_.size(), 0.0);

    if (!localRegularization_) {
        DOSAVE echoMinMax(model_, "model: ");
        roughness = this->roughness();

        if (haveReferenceModel_) {
            if (verbose_) echoMinMax(modelRef_,  "reference model");
            roughness = roughness - constraintsH_;
        }
    } else {
        if (verbose_) std::cout << "use local regularization" << std::endl;
    }

    DOSAVE echoMinMax(modelWeight_,  "mW");
    DOSAVE echoMinMax(constraintsH_, "constraintsH");
    DOSAVE save(constraintsH_, "constraintsH");

//    solveCGLSCDWWtrans(*J_, forward_->constraints(), dataWeight_, deltaData, deltaModel, constraintsWeight_,
//                        modelWeight_, tM_->deriv(model_), tD_->deriv(response_),
//                        lambda_, deltaModel0, maxCGLSIter_, verbose_);
    solveCGLSCDWWhtrans(*forward_->jacobian(), *forward_->constraints(),
                        dataWeight_, deltaDataIter_, deltaModel,
                        constraintsWeight_, modelWeight_,
                        tM_->deriv(model_), tD_->deriv(response_),
                        lambda_, roughness, maxCGLSIter_, verbose_);

    Vec appModelStart(tM_->invTrans(tModel + deltaModel));
    DOSAVE save(appModelStart, "appModel");
    double phiMNorm = (getPhiM(appModelStart));
    double phiDNorm = (getPhiD());

    //* normalization
//    phiM.push_back(getPhiM() / phiMNorm);
//    phiD.push_back(1.0);
    //* override normalization
    phiM.push_back(std::log(getPhiM(appModelStart))); phiMNorm = 1.0;
    phiD.push_back(std::log(getPhiD())); phiDNorm = 1.0;
    if(verbose_) std::cout << "lambda(0) = inf" << " PhiD = " << phiD.back() << " PhiM = " << phiM.back()  << std::endl;

    int lambdaIter = 0;
    while (lambdaIter < 30) {
        lambdaIter++;
        if(verbose_) std::cout << lambdaIter << "lambda = " << lambda_ << std::endl;
//        solveCGLSCDWWtrans(*J_, forward_->constraints(), dataWeight_, deltaData, deltaModel, constraintsWeight_,
//                          modelWeight_, tM_->deriv(model_), tD_->deriv(response_),
//                          lambda_, deltaModel0, maxCGLSIter_, verbose_);
        solveCGLSCDWWhtrans(*forward_->jacobian(), *forward_->constraints(),
                            dataWeight_, deltaDataIter_, deltaModel, 
                            constraintsWeight_, modelWeight_,
                            tM_->deriv(model_), tD_->deriv(response_),
                            lambda_, roughness, maxCGLSIter_, verbose_);

        Vec appModel(tM_->invTrans(tModel + deltaModel));
        Vec appResponse(tD_->invTrans(tResponse + *forward_->jacobian() * deltaModel));

        if (lambdaIter == 1) { //* normalize on 1st iteration
//      	 phiMNorm = getPhiM(appModel);
//      	 phiDNorm = getPhiD(appResponse);
        }

        phiM.push_back(std::log(getPhiM(appModel)) / phiMNorm);
        phiD.push_back(std::log(getPhiD(appResponse)) / phiDNorm);

        if (lambdaIter > 1) {
            ys = (phiD[ lambdaIter ] - phiD[ lambdaIter - 2 ]) / (phiM[ lambdaIter ] - phiM[ lambdaIter - 2 ]);
            yss = ((phiD[ lambdaIter ] - phiD[ lambdaIter - 1 ]) / (phiM[ lambdaIter ] - phiM[ lambdaIter - 1]) -
              (phiD[ lambdaIter - 1 ] - phiD[ lambdaIter - 2 ]) / (phiM[ lambdaIter - 1 ] - phiM[ lambdaIter - 2 ])) /
              (phiM[ lambdaIter ] - phiM[ lambdaIter - 2 ]) * 2.0;
            curv = yss / std::pow(1 + ys * ys, 1.5);

            if(verbose_) std::cout << " lambda(" << lambdaIter << ") = " << lambda_ << " PhiD = "
                                     << phiD.back() << " PhiM = " << phiM.back() << " curv = " << curv << std::endl;
            if ((curv < oldcurv) && (lambdaIter > 4)) {
                deltaModel = uroldDModel;
                lambda_ /= (0.8 * 0.8);
                if (verbose_) std::cout << lambdaIter << ": lambdaIter -- " << "Curvature decreasing, choosing lambda = " << lambda_ << std::endl;
                break;
            }
        oldcurv = curv;
        } else { //** lambdaIter > 1
            if(verbose_) std::cout << "lambda(" << lambdaIter << ") = " << lambda_
                                     << " PhiD = " << phiD.back() << " PhiM = " << phiM.back()  << std::endl;
        }
        uroldDModel = oldDModel;
        oldDModel = deltaModel;
        lambda_ *= 0.8;
    }  //** while loop;
    DOSAVE save(phiM, "phiM");
    DOSAVE save(phiD, "phiD");

    return deltaModel;
}

} // namespace GIMLI

#endif
