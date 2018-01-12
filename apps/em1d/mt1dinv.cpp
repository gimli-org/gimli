/***************************************************************************
 *   Copyright (C) 2008-2009 by the resistivity.net development team       *
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
#include <em1dmodelling.h>
#include <inversion.h>

#include <string>
using namespace GIMLI;

#define vcout if (verbose) std::cout
#define dcout if (debug) std::cout
#define DEBUG if (debug)


int main(int argc, char *argv []) {
    bool lambdaOpt = false, doResolution = false, isRobust = false;
    double lambda = 10.0, lbound = 0.0, ubound = 0.0, errPerc = 5.0, errPhase = 5.;
    int maxIter = 10, nlay = 3, verboseCount = 0;
    std::string modelFile(NOT_DEFINED), dataFileName;
    GIMLI::Placeholder x__;

    OptionMap oMap;
    oMap.setDescription("Description. MT1dInv - 1D block or smooth inversion of dc resistivity data\n");
    oMap.addLastArg(dataFileName, "Data file");
    oMap.add(verboseCount, "v" , "verbose", "Verbose/debug/dosave mode (multiple use).");
    oMap.add(lambdaOpt,    "O" , "OptimizeLambda", "Optimize model smoothness using L-curve.");
    oMap.add(doResolution, "D" , "doResolution", "Do resolution analysis.");
    oMap.add(isRobust,     "R" , "RobustData"    , "Robust (L1) data weighting.");
    oMap.add(lambda,       "l:", "lambda", "Regularization strength lambda[100].");
    oMap.add(errPerc,      "e:", "error", "Error percentage[0]");
    oMap.add(errPhase,     "p:", "PhaseError", "Phase error/degree[1]");
    oMap.add(nlay,         "n:", "nlay", "Number of layers[3]");
    oMap.add(maxIter,      "i:", "maxIter", "Number of Iterations[10]");
    oMap.add(lbound,       "b:", "lbound", "Lower Resistivity bound[0]");
    oMap.add(ubound,       "u:", "ubound", "Upper Resistivity bound[0-inactive]");
    oMap.add(modelFile,    "m:", "modelFile", "Model file for pure modelling");
    oMap.parse(argc, argv);

    bool verbose = (verboseCount > 0), debug = (verboseCount > 1), dosave = (verboseCount > 2);

    RMatrix TRP;
    loadMatrixCol(TRP, dataFileName);
    size_t nperiods = TRP.cols();
    std::cout << "nperiods " << nperiods << std::endl;

    Mesh mesh;
    RVector model;
    double medrhoa = median(TRP[1]);
    std::cout << "medrhoa " << medrhoa << std::endl;
    /*! Transformations: log for app. resisivity and thickness, logLU for resistivity */
    RTransLogLU transThk(0.1, 1e5);
    RTransLogLU transRho(lbound, ubound);
    RTransLog transRhoa; //! log apparent resistivity
    RTrans transPhi; //! linear phases
    RTransCumulative transData; //! combination of two trans functions
    transData.add(transRhoa, nperiods);
    transData.add(transPhi, nperiods);
    RVector error(cat(RVector(TRP[1] * errPerc / 100.0), RVector(nperiods, errPhase)));
    save(error, "error.vec");
//     size_t nModel = 2 * nlay - 1;
    MT1dModelling f(TRP[0], nlay, debug);
    double medskindepth = sqrt(median(TRP[0]) * medrhoa) *503.0;
    model = cat(RVector(nlay - 1, medskindepth / nlay), RVector(nlay, medrhoa));
    model[nlay] = model[nlay - 1] * 2; // inhomogeneous model
    save(model, "start.vec");
    save(f.response(model), "response.vec");

    f.region(0)->setTransModel(transThk);
    f.region(1)->setTransModel(transRho);

    /*! Set up inversion with full matrix, data and forward operator */
    Inversion< double > inv(cat(TRP[1], TRP[2]), f, verbose, dosave); //! only rhoa
    inv.setTransData(transRhoa);
    inv.setMarquardtScheme(0.8); //! Marquardt method
    inv.setLambda(lambda);
    inv.setOptimizeLambda(lambdaOpt);
    inv.setRobustData(isRobust);
    inv.setMaxIter(maxIter);
    inv.setAbsoluteError(error); //! error model
    inv.setModel(model);         //! starting model

    /*! actual computation: run the inversion */
    model = inv.run();
    save(model, "model.vec");
    save(inv.response(), "response.vec");
    
    RVector thk(nlay - 1);
    RVector res(nlay);
    for (int i = 0 ; i < nlay - 1 ; i++) thk[i] = model[i];
    for (int i = 0 ; i < nlay ; i++) res[i] = model[nlay - 1 + i];
    save(res, "resistivity.vec");
    save(thk, "thickness.vec");
    
    if (verbose) {
        RVector cumz(thk);
        for (size_t i = 1 ; i < thk.size() ; i++) cumz[i] = cumz[i-1] + thk[i];
		cumz.round(0.1);        
		res.round(0.1);        

        std::cout << "Res = " << res << std::endl;
        std::cout << "  z =  " << cumz << std::endl;
    }
    
    if (doResolution) {
        RVector resolution(model.size());
        RVector resMDiag (model.size());
        RMatrix resM;
        for (size_t iModel = 0; iModel < model.size(); iModel++) {
            resolution = inv.modelCellResolution(iModel);
            resM.push_back(resolution);
            resMDiag[iModel] = resolution[iModel];
        }
        save(resMDiag, "resMDiag.vec");
        save(resM,     "resM");
    }

    return EXIT_SUCCESS;
}
