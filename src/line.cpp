/******************************************************************************
 *   Copyright (C) 2008-2018 by the GIMLi development team                    *
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

#include "line.h"

namespace GIMLI{

std::ostream & operator << (std::ostream & str, const Line & l){
    if (l.valid()){
        str << "Line: " << l.p0() << "-- " << l.p1() << " ";
    } else {
        str << "Line: invalid ";
    }
    return str;
}

Line::Line()
    : valid_(false) {
}

Line::Line(const RVector3 & p)
    : p1_(p), valid_(false){
    p0_[0] = 0.0; p0_[1] = 0.0; p0_[2] = 0.0;
    valid_ = false;
    checkValidity();
}

Line::Line(const RVector3 & p0, const RVector3 & p1)
    : p0_(p0), p1_(p1), valid_(false){
    checkValidity();
}

Line::Line(const Line & line){
    copy_(line);
}

Line::~Line(){
}

Line & Line::operator = (const Line & line){
    if (this != & line){
        copy_(line);
    }
    return *this;
}

void Line::copy_(const Line & line){
    p0_ = line.p0();
    p1_ = line.p1();
    valid_ = line.valid();
}

bool Line::checkValidity(double tol){
    if (p0_.distance(p1_) > tol) valid_ = true; else valid_ = false;
    return valid_;
}

bool Line::compare(const Line & line, double epsilon) const {
    if (this->touch(line.p0()) && this->touch(line.p1())) return true;
    else return false;
}

RVector3 Line::intersect(const Line & line) const {
    THROW_TO_IMPL
    //check paralell and equal

  //** in 3D: wenn sie sich ber�hren spannen sie eine Ebene auf.

// x0 + s * x1 = x2 + t * x3;
// x0 - x2 = t * x3 - s * x1;
//  RVector3 x0(this->p0_);
//  RVector3 x1(this->p1_ - this->p0_);
//
//  RVector3 x2(line.p0());
//  RVector3 x3(line.p1() - line.p0());
//
//  RVector3 b(x0 - x2);

//  STLMatrix A(3, 2);
//  A.setColumn(0, -1.0 * x1);
//  A.setColumn(1, x3);
 //RVector3 st(3);

 //** simple check if 2D
/* if (x1[2] == 0.0 && x3[2] == 0.0){
   MyVec::STLVector b2(2);
   b2[0] = b[0];
   b2[1] = b[1];
   MyVec::STLMatrix A2(2, 2);
   A2[0][0] = A[0][0]; A2[1][0] = A[1][0];
   A2[0][1] = A[0][1]; A2[1][1] = A[1][1];

   solveLU(A2, st, b2);
//    cout << A2 << endl;
//    cout << b2[0] << " " << b2[1] << endl;

 } else {
   // L1 = P1 + a V1
   // L2 = P2 + b V2
   //   a (V1 X V2) = (P2 - P1) X V2

   RVector3 P1(this->p0());
   RVector3 V1(this->p1() - this->p0());
   RVector3 P2(line.p0());
   RVector3 V2(line.p1() - line.p0());

   double a = (P2-P1).cross(V2).abs() / V1.cross(V2).abs();

   if (!isnan(a) && !isnan(a)){
     return RVector3(P1 + V1 * a);
   }

   solveLU(A, st, b);*/
//     cout << A << endl;
//     cout << st[0] <<  " " << st[1] << endl;
//     cout << b[0] << " " << b[1] << " " << b[3] << endl;
// }

 // double s = st[0];
 // double t = st[1];
//
//  if (isnan(s) || isinf(s) || isnan(t) || isinf(t)) {
//    //   cout << WHERE_AM_I << " s: " << s << " t: " << t << endl;
//    return RVector3();
//  }
//

    //RVector3 iPos(x0 + x1 * s);
//
//  if (line.touch(iPos) < 0){
//    //cout << WHERE_AM_I << " touch: " << line.touch(iPos) << endl;
    return RVector3();
//  }

//    return iPos;
}

double Line::nearest(const RVector3 & p) const{
/*    double t = ((p.x - p0_.x) * (p1_.x - p0_.x)) +
               ((p.y - p0_.y) * (p1_.y - p0_.y)) +
               ((p.z - p0_.z) * (p1_.z - p2_.z))*/
    return (p-p0_).dot(p1_-p0_)/p1_.distSquared(p0_);
}

double Line::distance(const RVector3 & pos) const {
    //  cout << p1_ << p0_ << pos << endl;
    //cout << p1_ - p0_ << endl;
    //cout << (p1_ - p0_).cross(p0_ -pos) << endl;
    //cout << "ret " << ((p1_ - p0_).cross(p0_ - pos)).abs() / (p1_ - p0_).abs() << endl;
    return ((p1_ - p0_).cross(p0_ - pos)).abs() / (p1_ - p0_).abs();
}

double Line::t(const RVector3 & pos, double tol) const {
    RVector3 t(RVector3((pos - p0_).round(tol) / (p1_ - p0_).round(tol)));
//   cout << tol << endl;
//   cout << (pos - p0_).round(tol) << (p1_ - p0_).round(tol) << endl;
//   cout << RVector3((pos - p0_) / (p1_ - p0_)) << endl;
//   cout << WHERE_AM_I << t << endl;

    if (!isnan(t[0]) && !isinf(t[0])) return t[0];
    else if (!isnan(t[1]) && !isinf(t[1])) return t[1];
    else if (!isnan(t[2]) && !isinf(t[2])) return t[2];

    throwError(1, WHERE_AM_I + " pos is not at this line ");

    return 0.0;
}

bool Line::touch(const RVector3 & pos, double tol) const {
    if (this->distance(pos) > tol) return false;
    return true;
}

bool Line::touch1(const RVector3 & pos, int & pFunIdx, double tol) const{
    bool verbose = false;
    if (verbose) std::cout << pos << tol << std::endl;
    if (verbose) std::cout << "Dist = " << std::fabs(this->distance(pos)) << std::endl;
    double length = p0_.distance(p1_);
    double dist = std::fabs(this->distance(pos));
    if (length > 1) tol*= length;
    if (dist > (tol) * 10) {
        if (verbose) std::cout << "dist:" << dist << " length: " << length << " > " << "tol: " << (tol) * 10 << std::endl;
        pFunIdx = -1;
        return false;
    }

    //  cout << "next" <<endl;
    if ((dist) < tol) dist = tol;
    double tsol = this->t(pos, dist);
    if (verbose) std::cout << p0_ << pos << p1_ << "Tsol: " << tsol << std::endl;

    if (std::fabs(tsol) < tol) pFunIdx =  2;
    if (tsol < 0) pFunIdx = 1;
    if (std::fabs(1 - tsol) < tol) pFunIdx = 4;
    if (tsol > 1) pFunIdx = 5;

    pFunIdx = 3;
    return true;
}

} // namespace GIMLI

