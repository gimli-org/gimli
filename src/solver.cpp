/******************************************************************************
 *   Copyright (C) 2006-2020 by the GIMLi development team                    *
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
#include "gimli.h"
#include "matrix.h"
#include "vectortemplates.h"
#include "solver.h"

namespace GIMLI{

int solveCGLSCDWWhtrans(const MatrixBase & S, const MatrixBase & C,
                        const Vec & dWeight, 
                        const Vec & b, // deltaData
                        Vec & x,       // deltaModel
                        const Vec & wc, const Vec & wm,
                        const Vec & tm, // m/d mod
                        const Vec & td, // d/d resp
                        double lambda, const Vec & roughness,
                        int maxIter, double tol,
                        bool verbose){ //ALLOW_PYTHON_THREADS
    
    uint nData = b.size();
    uint nModel = x.size();
    uint nConst = C.rows();

//     __MS(S.rtti())
//     __MS(C.rtti())

    if (S.rows() != nData)  std::cerr << "J.rows != nData " << S.rows() << " / " << nData << std::endl;
    //if (S.cols() != nModel) std::cerr << "J.cols != nModel " << S.cols() << " / " << nModel << std::endl;
    if (C.cols() != nModel) std::cerr << "C.cols != nModel " << C.cols() << " / " << nModel << std::endl;
    if (dWeight.size() != nData) std::cerr << "dWeight.size() != nData" << dWeight.size() << " != " << nData << std::endl;
    if (wc.size() != nConst) std::cerr << "wc.size() != nConst " << wc.size() << " / " << nConst << std::endl;
    if (wm.size() != nModel) std::cerr << "wm.size() != nModel" << wm.size() << " / " << nModel << std::endl;
    if (tm.size() != nModel) std::cerr << "tm.size() != nModel " << tm.size() << " / " << nModel << std::endl;
    if (td.size() != nData) std::cerr << "td.size() != nData " << td.size() << " / " << nData << std::endl;
    if (roughness.size() != nConst) std::cerr << "roughness.size != nConst " << roughness.size() << " / " << nConst << std::endl;

//Ch  Vec cdx(transMult(C, Vec(wc * wc * (C * Vec(wm * deltaX)))) * wm * lambda); // nModel
    Vec cdx(transMult(C, Vec(wc * roughness)) * wm * lambda); // nModel
    Vec z((b - S * Vec(x / tm) * td) * dWeight); // nData
    Vec p(transMult(S, Vec(z * dWeight * td)) / tm   
          - transMult(C, Vec(wc * wc * (C * Vec(wm * x)))) * wm * lambda
          - cdx );// nModel
    Vec r(transMult(S, Vec(b * dWeight * dWeight * td)) / tm - cdx); // nModel

// __MS(min(dWeight) << " " << max(dWeight) << " " << mean(dWeight))
// __MS(min(b) << " " << max(b) << " " << mean(b))
// __MS(min(wc) << " " << max(wc) << " " << mean(wc))
// __MS(min(wm) << " " << max(wm) << " " << mean(wm))
// __MS(min(tm) << " " << max(tm) << " " << mean(tm))
// __MS(min(td) << " " << max(td) << " " << mean(td))
// __MS(min(cdx) << " " << max(cdx) << " " << mean(cdx))
// __MS(min(z) << " " << max(z) << " " << mean(z))
// __MS(min(p) << " " << max(p) << " " << mean(p))
// __MS(min(r) << " " << max(r) << " " << mean(r))
// __MS(min(Vec(b * dWeight * dWeight * td)) << " " 
//      << max(Vec(b * dWeight * dWeight * td)) << " " 
//      << mean(Vec(b * dWeight * dWeight * td)))
// __MS(min(transMult(S,Vec(b.size(), 1))) << " " 
//      << max(transMult(S,Vec(b.size(), 1))) << " " 
//      << mean(transMult(S,Vec(b.size(), 1))))
// if (z.size() > 100 )exit(1);

    double accuracy = tol;
    if (accuracy < 0.0) accuracy = max(TOLERANCE, 1e-08 * dot(r, r));
    r = p;

    double normR2 = dot(r, r), normR2old = 0.0;
    double alpha = 0.0, beta = 0.0;

    int count = 0;

    Vec q(nData);
    Vec wcp(nConst); // nBounds

//     std::cout.precision(14);
//     std::cout << 0 << "  " << accuracy << std::endl;

    while (count < maxIter && normR2 > accuracy){
        count ++;
        q = S * Vec(p / tm) * dWeight * td;
        wcp = wc * (C * Vec(p * wm));

        alpha = normR2 / (dot(q, q) + lambda * dot(wcp, wcp));
        x += p * alpha;

        if((count % 10) == -1) { // TOM
            z = b - S * Vec(x / tm) * td; // z exakt durch extra Mult.
            z *= dWeight;
        } else {
            z -= q * alpha;
        }
        r = transMult(S, Vec(z * dWeight * td)) / tm 
            - transMult(C, Vec(wc * wc * (C * Vec(wm * x)))) * wm * lambda 
            - cdx;

        normR2old = normR2;
        normR2 = dot(r, r);
        beta = normR2 / normR2old;
        p = r + p * beta;

// __MS("---" << count << "-------" << normR2 << "-------------------------")
// __MS(min(q) << " " << max(q) << " " << mean(q))
// __MS(min(wcp) << " " << max(wcp) << " " << mean(wcp))
// __MS(alpha)

//         if((count < 200)) {
//             std::cout << count << "  " << normR2 << std::endl;
//         }

#ifndef _WIN32
        if (verbose) std::cout << "\r[ " << count << "/" << normR2 << "]\t";
#endif
    }
#ifdef _WIN32
    if (verbose) std::cout << "[ " << count << "/" << normR2 << "]\t" << std::endl;
#endif

    // if (z.size() > 100 )exit(1);
    return 1;
}

int solveCGLSCDWWtrans(const MatrixBase & S, const MatrixBase & C,
                       const Vec & dWeight,  const Vec & b, Vec & x,
                       const Vec & wc, const Vec & mc, const Vec & tm,
                       const Vec & td, double lambda, const Vec & deltaX,
                       int maxIter, bool verbose){ //ALLOW_PYTHON_THREADS

  uint nData = b.size();
  uint nModel = x.size();
  uint nConst = C.rows();

  if (S.rows() != nData)  std::cerr << "J.rows != nData " << S.rows() << " / " << nData << std::endl;
  // not needed
  // if (S.cols() != nModel) std::cerr << "J.cols != nModel " << S.cols() << " / " << nModel << std::endl;
  if (C.rows() != nConst) std::cerr << "C.rows != nConst " << C.rows() << " / " << nConst << std::endl;
  if (C.cols() != nModel) std::cerr << "C.cols != nModel " << C.cols() << " / " << nModel << std::endl;
  if (dWeight.size() != nData) std::cerr << "dWeight.size() != nData" << dWeight.size() << " != " << nData << std::endl;
  if (wc.size() != nConst) std::cerr << "wc.size() != nConst " << wc.size() << " / " << nConst << std::endl;
  if (mc.size() != nModel) std::cerr << "mc.size() != nModel" << mc.size() << " / " << nModel << std::endl;
  if (tm.size() != nModel) std::cerr << "tm.size() != nModel " << tm.size() << " / " << nModel << std::endl;
  if (td.size() != nData) std::cerr << "td.size() != nData " << td.size() << " / " << nData << std::endl;

  Vec cdx(transMult(C, Vec(wc * wc * (C * Vec(mc * deltaX)))) * mc * lambda); // nModel
  Vec z((b - S * Vec(x / tm) * td) * dWeight); // nData
  Vec p(transMult(S, Vec(z * dWeight * td)) / tm - cdx - transMult(C, Vec(wc * wc * (C * Vec(mc * x)))) * mc * lambda);// nModel
  Vec r(transMult(S, Vec(b * dWeight * dWeight * td)) / tm - cdx); // nModel

  double accuracy = max(TOLERANCE, 1e-08 * dot(r, r));
  r = p;

  double normR2 = dot(r, r), normR2old = 0.0;
  double alpha = 0.0, beta = 0.0;

  int count = 0;

  Vec q(nData);
  Vec wcp(nConst); // nBounds

  while (count < maxIter && normR2 > accuracy){

    count ++;
    q = S * Vec(p / tm) * dWeight * td;
    wcp = wc * (C * Vec(p * mc));

    alpha = normR2 / (dot(q, q) + lambda * dot(wcp, wcp));
    x += p * alpha;
    if((count % 10) == -1) { // TOM
      z = b - S * Vec(x / tm) * td; // z exakt durch extra Mult.
      z *= dWeight;
    } else {
      z -= q * alpha;
    }
    r = transMult(S, Vec(z * dWeight * td)) / tm - transMult(C, Vec(wc * wc * (C * Vec(mc * x)))) * mc * lambda - cdx;

    normR2old = normR2;
    normR2 = dot(r, r);
    beta = normR2 / normR2old;
    p = r + p * beta;
    #ifndef _WIN32
    if (verbose) std::cout << "\r[ " << count << "/" << normR2 << "]\t";
    #endif
  }
  #ifdef _WIN32
  if (verbose) std::cout << "[ " << count << "/" << normR2 << "]\t" << std::endl;
  #endif
  return 1;
}

} //namespace GIMLI

