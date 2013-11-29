/***************************************************************************
 *   Copyright (C) 2006-2013 by the resistivity.net development team       *
 *   Thomas Günther thomas@resistivity.net                                 *
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

#include <gimli.h>
#include <optionmap.h>
#include <mesh.h>
#include <datacontainer.h>
#include <node.h>
#include <inversion.h>

#include <ttdijkstramodelling.h>

using namespace GIMLI;

#include <string>

#define vcout if (verbose) std::cout
#define dcout if (debug) std::cout
#define DEBUG if (debug)

//** MAIN
int main(int argc, char *argv []) {

    //!** Initialisation of mesh and inversion options;
// //     int nSegments = 6;
//     double relativeInnerMaxEdgeLength = 0.01;
    bool lambdaOpt = false, isBlocky = false, isRobust = false, createGradientModel = false;
    bool useAppPar = false, isBertMesh = false, isVelocity = false, refineFwdMesh = false;
    double errTime = 0.001, errPerc = 0.1, lbound = 0.0, ubound = 0.0;
    double lambda = 100.0, zWeight = 1.0;
    int maxIter = 20, verboseCount = 0;
    std::string paraMeshFilename = NOT_DEFINED, startModelFilename = NOT_DEFINED, dataFileName;

    //! Parse command line
    OptionMap oMap;
    oMap.setDescription("Description. TTInv - Travel time inversion using Dijkstra algorithm.\n");
    oMap.addLastArg(dataFileName, "Data file");
    oMap.add(verboseCount,       "v" , "verbose", "Verbose mode (2 times for debug mode).");
    oMap.add(isBertMesh,         "1" , "isBertMesh", "The mesh is a BERT mesh with 1 background");
    oMap.add(useAppPar,          "A" , "UseAppPar", "Use Apparent Parameters.");
    oMap.add(isBlocky,           "B" , "BlockyModel", "Blocky (L1) model constraints.");
    oMap.add(createGradientModel,"G" , "gradientModel", "Gradient starting model.");
    oMap.add(isRobust,           "R" , "RobustData", "Robust (L1) data weighting.");
    oMap.add(lambdaOpt,          "O" , "OptimizeLambda", "Optimize model smoothness using L-curve.");
    oMap.add(refineFwdMesh,      "S" , "refineFwdMesh", "Refine mesh for forward calculation.");
    oMap.add(isVelocity,         "V" , "isVelocity", "boundaries are velocities (not slowness).");
    oMap.add(lambda,             "l:", "lambda", "Regularization parameter lambda.");
    oMap.add(maxIter,            "i:", "iterations", "Maximum iteration number.");
    oMap.add(errPerc,            "e:", "errorPerc", "Percentage error part.");
    oMap.add(errTime,            "t:", "errorTime", "Absolute error for first arrival times.");
    oMap.add(startModelFilename, "s:", "startModel", "Starting model.");
    oMap.add(zWeight,            "z:", "zWeight", "z weight (vertical smoothness)[1=isotropic]");
    oMap.add(paraMeshFilename,   "p:", "paraMeshFile", "Parameter mesh file.");
    oMap.add(lbound,             "b:", "lowerBound", "Lower velocity/slowness bound.");
    oMap.add(ubound,             "u:", "upperBound", "Upper velocity/slowness bound(0=inactive).");
    oMap.parse(argc, argv);

    bool verbose = (verboseCount > 0), debug = (verboseCount > 1);
    if (verbose) setbuf(stdout, NULL);

    if ((ubound > 0.) && (lbound > ubound)) isVelocity = true; // apparently velocities are meant
    if (isVelocity) { // velocity specified instead of slowness
        double dummy = lbound;
        if (ubound > 0.0) { // both bounds are given
            lbound = 1.0 / ubound;
            ubound = 1.0 / dummy;
        } else if (lbound > 0.0) {
            lbound = 1.0 / dummy;
        }
    }
    vcout << "Lower/upper slowness bound = " << lbound << "/" << ubound << std::endl;

    //!** load data file
    DataContainer dataIn(dataFileName, "s g"); //! predefine sensor indices and load
    if (verbose) dataIn.showInfos();

    //!** apply error model if not defined;
    if (!dataIn.allNonZero("err")) {
        vcout << "Estimate error: " << errPerc << "% + " 
              << errTime << "s" << std::endl;
        dataIn.set("err", dataIn("t")*  errPerc / 100.0 + errTime); // always relative error
    }
    vcout << "Data error: min = " << min(dataIn("err")) * 1000. 
          << "ms max = " << max(dataIn("err")) * 1000. << "ms" << std::endl;

    //!** Load mesh
    Mesh paraMesh;
    if (paraMeshFilename != NOT_DEFINED) {
        paraMesh.load(paraMeshFilename);
    } else {
        std::cerr << "No mesh given!" << std::endl;
        return 1;
    }
    vcout << "Paramesh: " << paraMesh << std::endl;

    //!** set up TT modeling class;
    TravelTimeDijkstraModelling f(paraMesh, dataIn, verbose);
    RVector appSlowness(f.getApparentSlowness());
    vcout << "min/max apparent velocity = " << 1.0 / max(appSlowness) 
          << " / " << 1.0 / min(appSlowness) << " m/s" << std::endl;

    //!** get mesh from region manager (BERT convention: first/smallest region is background
    if (isBertMesh && f.regionManager().regionCount() > 1) {
        f.regionManager().regions()->begin()->second->setBackground(true);
    }

    Mesh paraDomain(f.regionManager().paraDomain());
    vcout << "ParaDomain:\t"<< paraDomain << std::endl;
    DEBUG paraDomain.save("meshParaDomain.bms");
    DEBUG paraDomain.exportVTK("meshParaDomain");

    //!** necessary step to reflect changes due to possible region changes
    f.createRefinedForwardMesh(refineFwdMesh, false);

    size_t nModel = f.regionManager().parameterCount();
    vcout << "model size = " << nModel << std::endl;
    RVector startModel(nModel, median(appSlowness));

    if (startModelFilename != NOT_DEFINED) {
        load(startModel, startModelFilename);
    }

    if (createGradientModel) { // auslagern in ttdijkstramodel
        vcout << "Creating Gradient model ..." << std::endl;
        double smi = median(appSlowness);
        if (smi < lbound) smi = lbound * 1.1;
        
        double sma = max(appSlowness) / 2.0;
        if (ubound > 0.0 && sma > ubound) sma = ubound * 0.9;

        RVector zmid(nModel);
        int di = paraDomain.dim() - 1;
        for (size_t i = 0; i < paraDomain.cellCount() ; i++) {
            for (size_t j = 0; j < paraDomain.cell(i).nodeCount(); j++) {
                zmid[ i ] += paraDomain.cell(i).node(j).pos()[ di ];
            }
            zmid[ i ] /= double(paraDomain.cell(i).nodeCount());
        }
        double zmi = min(zmid);
        double zma = max(zmid);
        for (size_t i = 0; i < startModel.size(); i++) {
            startModel[ i ] = smi * std::exp((zmid[ i ] - zmi) / (zma - zmi) * std::log(sma / smi));
        }
        DEBUG save(startModel, "startModel.vector");
    }

    f.setStartModel(startModel);

    double sloRef = max(appSlowness);
    if (sloRef < lbound || (ubound > 0.0 && sloRef > ubound)) {
        sloRef = std::max(std::sqrt(lbound * ubound), (lbound + ubound) / 2);
    }

    vcout << "Start model size = " << startModel.size() 
            << " min/max = " << min(startModel) << "/" 
            << max(startModel) << std::endl;

    RTransLogLU transModel(lbound, ubound);
    RTrans transData;

    //!** Model transformation: log slowness with lower and upper bound
    Inversion < double > inv(dataIn("t"), f, transData, transModel, verbose, debug);
    inv.setModel(startModel);
    inv.setLambda(lambda);
    inv.setMaxIter(maxIter);
    inv.setAbsoluteError(dataIn("err"));
    inv.setOptimizeLambda(lambdaOpt);
    inv.setRobustData(isRobust);
    inv.setBlockyModel(isBlocky);
    inv.stopAtChi1(false);

    vcout << "setting z weight...";
    f.regionManager().setZWeight(zWeight);
    vcout << "ready. Start inversion." << std::endl;

    //!** Start Inversion
    RVector model(inv.run());
    vcout << "inversion ready." << std::endl;

    //!** Save velocity, model response and coverage
    RVector velocity(1.0 / (model + TOLERANCE));
    save(velocity, "velocity.vector");
    save(inv.response(), "response.vector");

    RVector one(dataIn.size(), 1.0);
    RVector coverage(transMult(*f.jacobian(), one));
    save(coverage, "coverage.vector");

    RVector scoverage(transMult(inv.constraintMatrix(), inv.constraintMatrix() * coverage));
    for (Index i = 0 ; i < scoverage.size() ; i++){
        scoverage[ i ] = sign(std::abs(scoverage[ i ]));
    }
    save(scoverage, "scoverage.vector");

    //!** Save VTK file with the vectors
    paraDomain.addExportData("slowness" , model);
    paraDomain.addExportData("velocity" , velocity);
    paraDomain.addExportData("coverage" , coverage);
    paraDomain.addExportData("scoverage" , scoverage);
    paraDomain.exportVTK("ttinv.result");

    //!** Cleanup and exit
    return EXIT_SUCCESS;
}
