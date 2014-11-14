/***************************************************************************
 *   Copyright (C) 2006-2014 by the resistivity.net development team       *
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

#ifndef _GIMLI_SOLVER__H
#define _GIMLI_SOLVER__H

#include "gimli.h"
#include "matrix.h"
#include "vectortemplates.h"

namespace GIMLI{

template < class Vec >
int solveCGLSCDWWhtrans(const MatrixBase & S, const MatrixBase & C,
                        const Vec & dWeight, const Vec & b, Vec & x,
                        const Vec & wc, const Vec & wm, 
                        const Vec & tm, const Vec & td,
                        double lambda, const Vec & roughness, 
                        int maxIter=200, double tol = -1.0,
                        bool verbose=false){ //ALLOW_PYTHON_THREADS

    uint nData = b.size();
    uint nModel = x.size();
    uint nConst = C.rows();

    if (S.rows() != nData)  std::cerr << "J.rows != nData " << S.rows() << " / " << nData << std::endl;
    if (S.cols() != nModel) std::cerr << "J.cols != nModel " << S.cols() << " / " << nModel << std::endl;
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
    Vec p(transMult(S, Vec(z * dWeight * td)) / tm - cdx - transMult(C, Vec(wc * wc * (C * Vec(wm * x)))) * wm * lambda);// nModel
    Vec r(transMult(S, Vec(b * dWeight * dWeight * td)) / tm - cdx); // nModel

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
        r = transMult(S, Vec(z * dWeight * td)) / tm - transMult(C, Vec(wc * wc * (C * Vec(wm * x)))) * wm * lambda - cdx;

        normR2old = normR2;
        normR2 = dot(r, r);
        beta = normR2 / normR2old;
        p = r + p * beta;
        
//         if((count < 200)) { 
//             std::cout << count << "  " << normR2 << std::endl;
//         }

#ifndef MINGW
        if (verbose) std::cout << "\r[ " << count << "/" << normR2 << "]\t";
#endif
    }
#ifdef MINGW
    if (verbose) std::cout << "[ " << count << "/" << normR2 << "]\t" << std::endl;
#endif
    return 1;
}

template < class Vec >
int solveCGLSCDWWtrans(const MatrixBase & S, const MatrixBase & C, const Vec & dWeight,
		  const Vec & b, Vec & x, const Vec & wc, const Vec & mc, const Vec & tm, const Vec & td,
		  double lambda, const Vec & deltaX, int maxIter = 200, bool verbose = false){ //ALLOW_PYTHON_THREADS

  uint nData = b.size();
  uint nModel = x.size();
  uint nConst = C.rows();

  if (S.rows() != nData)  std::cerr << "J.rows != nData " << S.rows() << " / " << nData << std::endl;
  if (S.cols() != nModel) std::cerr << "J.cols != nModel " << S.cols() << " / " << nModel << std::endl;
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
    #ifndef MINGW
    if (verbose) std::cout << "\r[ " << count << "/" << normR2 << "]\t";
    #endif
  }
  #ifdef MINGW
  if (verbose) std::cout << "[ " << count << "/" << normR2 << "]\t" << std::endl;
  #endif
  return 1;
}

template < class Mat, class CMatrix, class Vec >
int solveCGLSCDWW(const Mat & S, const CMatrix & C, const Vec & dWeight,
		  const Vec & b, Vec & x, const Vec & wc, const Vec & mc,
		  double lambda, const Vec & deltaX, int maxIter = 200, bool verbose = false){ ALLOW_PYTHON_THREADS

  uint nData = b.size();
  uint nModel = x.size();
  uint nConst = C.rows();
  if (nData != S.rows()) std::cerr << "size(ddata) != J.rows " << nData << " != " << S.rows() << std::endl;
  if (nModel != S.cols()) std::cerr << "size(model) != J.cols " << nModel << " != " << S.cols() << std::endl;
  if (nModel != C.cols()) std::cerr << "size(model) != C.cols " << nModel << " / " << C.cols() << std::endl;
  if (nConst != wc.size()) std::cerr << "size(wmodel) != C.rows " << nConst << " / " << wc.size() << std::endl;
  if (nModel != mc.size()) std::cerr << "size(wconst) != C.cols " << nModel << " / " << mc.size() << std::endl;

//   if (nData != S.rows()) std::cerr << "Data-dimension differ SensM" << nData << " / " << S.rows() << std::endl;
//   if (nModel != S.cols()) std::cerr << "Modell-dimension differ SensM" << nModel << " / " << S.cols() << std::endl;
  //  if (nModel != C.size()) std::cerr << "Modell-dimension differ ContraintsM (critcal)" << nModel << " / " << C.size() << std::endl;

  Vec cdx(transMult(C, Vec((wc * wc) * (C * Vec(deltaX * mc)))) * mc * lambda); // nModel
  Vec z((b - S * x) * dWeight); // nData
  Vec p(transMult(S, Vec(z * dWeight)) - cdx - transMult(C, Vec(wc * wc * (C * Vec(mc * x)))) * mc * lambda);// nModel
  Vec r(transMult(S, Vec(b * dWeight * dWeight)) - cdx); // nModel
  double accuracy = 1e-8 * dot(r, r);
  r = p;

  double normR2 = dot(r, r), normR2old = 0.0;
  double alpha = 0.0, beta = 0.0;

  int count = 0;

  Vec q(nData);
  Vec wcp(nConst); // nBounds

  while(count < maxIter && normR2 > accuracy){

    count ++;
    q = S * p * dWeight;
    wcp = wc * (C * Vec(p * mc));

    alpha = normR2 / (dot(q, q) + lambda * dot(wcp, wcp));
    x += p * alpha;
    if((count % 10) == -1) { // TOM
      z = b - S*x; // z exakt durch extra Mult.
      z *= dWeight;
    } else {
      z -= q * alpha;
    }
    r = transMult(S, Vec(z * dWeight)) - transMult(C, Vec(wc * wc * (C * Vec(mc * x)))) * mc * lambda - cdx;

    normR2old = normR2;
    normR2 = dot(r, r);
    beta = normR2 / normR2old;
    p = r + p * beta;
    #ifndef MINGW
    if (verbose)  std::cout << "\r[ " << count << "/" << normR2 << "]\t";
    #endif
  }
  if (verbose) std::cout << "\r[ " << count << "/" << normR2 << "]\t" << std::endl;
  return 1;
}


template < class Mat, class CMatrix, class Vec >
int solveCGLSCDWtrans(const Mat & S, const CMatrix & C, const Vec & dWeight,
		  const Vec & b, Vec & x, const Vec & wc, const Vec & tm, const Vec & td,
    		double lambda, const Vec & deltaX, int maxIter = 200, bool verbose = false){ ALLOW_PYTHON_THREADS


    uint nData = b.size();
    uint nModel = x.size();
    uint nBounds = C.rows();

  // if (nData != S.size()) cerr << "Data-dimension differ SensM" << nData << " / " << S.size() << endl;
  // if (nModel != S[ 0 ].size()) cerr << "Modell-dimension differ SensM" << nModel << " / " << S[ 0 ].size() << endl;
  // if (nModel != C.size()) cerr << "Modell-dimension differ ConstraintsM (critcal)" << nModel << " / " << C.size() << endl;

  Vec cdx(transMult(C, Vec(wc * wc * (C * deltaX) * lambda))); // nModel
  Vec z((b - S * Vec(x / tm) * td) * dWeight); // nData
  Vec p(transMult(S, Vec(z * dWeight * td)) / tm - cdx - transMult(C, Vec(wc * wc * (C * x) * lambda)));// nModel
  Vec r(transMult(S, Vec(b * dWeight * dWeight * td)) / tm - cdx); // nModel

  double accuracy = 1e-08 * dot(r, r);
//  std::cout << "init: " << dot(r, r)<<  " " << dot(p, p)<< " " << accuracy << std::endl;

  r = p;

  double normR2 = dot(r, r), normR2old = 0.0;
  double alpha = 0.0, beta = 0.0;

  int count = 0;
  
  Vec q(nData);
  Vec wcp(nBounds); // nBounds
  Vec pbytm(nModel);
  Vec zd(nData);
  Vec wwcx(nBounds);

  while(count < maxIter && normR2 > accuracy){

    count ++;
    pbytm = p / tm;
    q = S * pbytm * dWeight * td;
    wcp = wc * (C * p);

    alpha = normR2 / (dot(q, q) + lambda * dot(wcp, wcp));
    x += p * alpha;
    if((count % 10) == -1) { // TOM
      z = b - S * Vec(x / tm) * td; // z exakt durch extra Mult.
      z *= dWeight;
    } else {
      z -= q * alpha;
    }
    zd = z * dWeight * td;
    wwcx = wc * wc * (C * x) * lambda;
    r = transMult(S, zd) / tm - transMult(C, wwcx) - cdx;

    normR2old = normR2;
    normR2 = dot(r, r);
    beta = normR2 / normR2old;
    p = r + p * beta;
    #ifndef MINGW
    if (verbose) std::cout << "\r[ " << count << "/" << normR2 << "]\t";
    #endif
  }

  if (verbose) std::cout << "\r[ " << count << "/" << normR2 << "]\t" << std::endl;
  return 1;
}



template < class Mat, class CMatrix, class Vec >
int solveCGLSCDW(const Mat & S, const CMatrix & C, const Vec & dWeight,
		  const Vec & b, Vec & x, const Vec & wc,
		  double lambda, const Vec & deltaX, bool verbose = false){

  //if (verbose) save(x, "x.vec");
  uint nData = b.size();
  //  uint nModel = x.size();
  uint nBounds = C.rows();

//   if (nData != S.rows()) std::cerr << "Data-dimension differ SensM" << nData << " / " << S.rows() << std::endl;
//   if (nModel != S.cols()) std::cerr << "Modell-dimension differ SensM" << nModel << " / " << S.cols() << std::endl;
  //  if (nModel != C.size()) std::cerr << "Modell-dimension differ ContraintsM (critcal)" << nModel << " / " << C.size() << std::endl;

  Vec cdx(transMult(C, Vec((wc * wc) * (C * deltaX) * lambda))); // nModel
  Vec z((b - S * x) * dWeight); // nData
  //** version 1, should be better include starting vector (lambda optimization)
  Vec p(transMult(S, Vec(z * dWeight)) - cdx - transMult(C, Vec(wc * wc * (C * x) * lambda)));// nModel
  Vec r(transMult(S, Vec(b * dWeight * dWeight)) - cdx); // nModel
  double accuracy = 1e-8 * dot(r, r);
  r = p;
  //** version 2, strictly translated from dcmatlab
  //Vec p(transMult(S, Vec(z * dWeight)));
  //double accuracy = dot(p, p) * 1e-10;
  //p = p - cdx - transMult(C, Vec(wc * wc * (C * x) * lambda));// nModel
  //Vec r(p);

  double normR2 = dot(r, r), normR2old = 0.0;
  double alpha = 0.0, beta = 0.0;

  int count = 0;
  int maxcount = 200;

  Vec q(nData);
  Vec wcp(nBounds); // nBounds

  while(count < maxcount && normR2 > accuracy){

    count ++;
    q = S * p * dWeight;
    wcp = wc * (C * p);

    alpha = normR2 / (dot(q, q) + lambda * dot(wcp, wcp));
    x += p * alpha;
    if((count % 10) == -1) { // TOM
      z = b - S*x; // z exakt durch extra Mult.
      z *= dWeight;
    } else {
      z -= q * alpha;
    }
    r = transMult(S, Vec(z * dWeight)) - transMult(C, Vec(wc * wc * (C * x) * lambda)) - cdx;

    normR2old = normR2;
    normR2 = dot(r, r);
    beta = normR2 / normR2old;
    p = r + p * beta;
#ifndef MINGW
    if (verbose)  std::cout << "\r[ " << count << "/" << normR2 << "]\t";
#endif
  }
  if (verbose) std::cout << "\r[ " << count << "/" << normR2 << "]\t" << std::endl;
  return 1;
}

} //namespace _GIMLI_SOLVER__H

#endif
