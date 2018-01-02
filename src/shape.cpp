/******************************************************************************
 *   Copyright (C) 2006-2018 by the GIMLi development team                    *
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

#include "shape.h"
#include "node.h"
#include "line.h"
#include "plane.h"
#include "vectortemplates.h"
#include "meshentities.h"

#include "inversion.h"

#if USE_BOOST_THREAD
        boost::mutex ShapeFunctionWriteCacheMutex__;
#else
        std::mutex ShapeFunctionWriteCacheMutex__;
#endif


namespace GIMLI{

template < > DLLEXPORT ShapeFunctionCache * Singleton < ShapeFunctionCache >::pInstance_ = NULL;

std::vector < PolynomialFunction < double > >
createPolynomialShapeFunctions(const std::vector < RVector3 > & pnts,
                               uint dim, uint nCoeff, bool pascale,
                               bool serendipity, const RVector & startVector){
// __M
    bool verbose = false;

    if (verbose){
        std::cout << "dim: " << dim << " pascale: " << pascale << " serendipity: " << serendipity << std::endl;
        std::cout << "start: " << startVector << std::endl;
        for (Index i = 0; i < pnts.size(); i ++) {
            std::cout << "P" << i << ": " << pnts[i] << std::endl;
        }
    }

    PolynomialModelling fop(dim, nCoeff, pnts, startVector);
    fop.setPascalsStyle(pascale);
    fop.setSerendipityStyle(serendipity);
//
//     if (verbose){
//         PolynomialFunction < double >tmp(RVector(fop.polynomialFunction().size()));
//         tmp.fill(fop.startModel());
//         std::cout << "base: " << tmp << std::endl;
//     }


    std::vector < PolynomialFunction < double > > ret;

    for (Index i = 0; i < pnts.size(); i ++){
        fop.clearJacobian();
        if (verbose) std::cout << "-pre-N" << i << ": " << fop.polynomialFunction() << std::endl;

        RVector N(pnts.size(), 0.0);
        N[i] = 1.0;
        RInversion inv(N, fop, verbose, verbose);
        inv.setRelativeError(0.0);
        inv.stopAtChi1(false);
        inv.setCGLSTolerance(1e-40);
        inv.setLambda(0);
        inv.setMaxIter(20);
        inv.run();

        if (verbose) std::cout << "N" << i << ": " << fop.polynomialFunction() << std::endl;
        ret.push_back(fop.polynomialFunction());

    }

    return ret;
}

std::ostream & operator << (std::ostream & str, const Shape & c){
    str << c.name() << " " << std::endl;
    for (uint i = 0; i < c.nodes().size(); i ++) {
        str << c.nodes()[i]->pos() << " ";
    }
    return str;
}

Shape::Shape(){
    domSize_ = 0.0;
    hasDomSize_ = false;
}

Shape::~Shape(){
}

void Shape::changed(){
    invJacobian_.clear();
    invJacobian_.setValid(false);
    hasDomSize_ = false;
}

Node & Shape::node(Index i) {
    if (i > nodeCount() - 1){
        std::cerr << WHERE_AM_I << " requested shape node: " << i << " does not exist." << std::endl;
        exit(EXIT_MESH_NO_NODE);
    }
    return *nodeVector_[i];
}

const Node & Shape::node(Index i) const {
    if (i > nodeCount() - 1){
        std::cerr << WHERE_AM_I << " requested shape node: " << i << " does not exist." << std::endl;
        exit(EXIT_MESH_NO_NODE);
    }
    return *nodeVector_[i];
}

void Shape::setNode(Index i, Node & n) {
    if (i > nodeCount() - 1){
        std::cerr << WHERE_AM_I << " requested shape node: " << i << " does not exist." << std::endl;
        exit(EXIT_MESH_NO_NODE);
    }
    nodeVector_[i] = &n;
    this->changed();
}

bool Shape::enforcePositiveDirection(){
    if (createJacobian().det() < 0){
        std::reverse(nodeVector_.begin(), nodeVector_.end());
        this->changed();
        return true;
    }
    return false;
}

RVector3 Shape::center() const {
    RVector3 center(0.0, 0.0, 0.0);
    for (uint i = 0; i < nodeVector_.size(); i ++) {
        center += nodeVector_[i]->pos();
    }
    center /= nodeVector_.size();
    return center;
}

RVector3 Shape::norm() const {
    THROW_TO_IMPL
    return RVector3();
}

std::vector < PolynomialFunction < double > > Shape::createShapeFunctions() const {
    uint nCoeff = 2;
    bool pascale = false;
    bool serendipity = false;

    switch (this->rtti()){
        case MESH_SHAPE_EDGE_RTTI:
        case MESH_SHAPE_TRIANGLE_RTTI:
        case MESH_SHAPE_TETRAHEDRON_RTTI:
            pascale = true;
            break;
        case MESH_SHAPE_QUADRANGLE_RTTI:
        case MESH_SHAPE_HEXAHEDRON_RTTI:
            pascale = true;
            serendipity = true;
            break;
    }

    return createPolynomialShapeFunctions(*this, nCoeff, pascale, serendipity);
}

RVector Shape::N(const RVector3 & rst) const {
    RVector n(nodeCount(), 0.0);
    N(rst, n);
    return n;
}

void Shape::N(const RVector3 & rst, RVector & n) const {
    const std::vector< PolynomialFunction < double > > &N = ShapeFunctionCache::instance().shapeFunctions(*this);

    for (Index i = 0; i < N.size(); i ++) {
        n[i] = N[i](rst);
    }
}

RMatrix Shape::dNdrst(const RVector3 & rst) const {
    RMatrix MdNdrst(3, nodeCount());
    dNdrst(rst, MdNdrst);
    return MdNdrst;
}

void Shape::dNdrst(const RVector3 & rst, RMatrix & MdNdrst) const {
    MdNdrst *= 0.0;

    const std::vector< PolynomialFunction < double > > &dNx = ShapeFunctionCache::instance().deriveShapeFunctions(*this, 0);
    const std::vector< PolynomialFunction < double > > &dNy = ShapeFunctionCache::instance().deriveShapeFunctions(*this, 1);
    const std::vector< PolynomialFunction < double > > &dNz = ShapeFunctionCache::instance().deriveShapeFunctions(*this, 2);

    for (Index i = 0; i < dNx.size(); i ++) {
         MdNdrst[0][i] = dNx[i](rst);
         if (this->dim() > 1) MdNdrst[1][i] = dNy[i](rst);
         if (this->dim() > 2) MdNdrst[2][i] = dNz[i](rst);
    }
}

RVector crossN(const RVector & a, const RVector & b){
    RVector c(3, 0.0);
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    return c / norm(c);
}

RMatrix3 Shape::createJacobian() const {
    RMatrix3 J;
    createJacobian(J);
    return J;
}

void Shape::createJacobian(RMatrix3 & J) const {
    RVector x(nodeCount());
    RVector y(nodeCount());
    RVector z(nodeCount());

    for (uint i = 0; i < nodeCount(); i ++){
        x[i] = nodeVector_[i]->pos()[0];
        y[i] = nodeVector_[i]->pos()[1];
        z[i] = nodeVector_[i]->pos()[2];
    }

    if (ShapeFunctionCache::instance().RMatrixCache(rtti()).size() < 1){
        ShapeFunctionCache::instance().RMatrixCache(rtti()).push_back(RMatrix(3, nodeCount()));
    }

    RMatrix & MdNdrst = ShapeFunctionCache::instance().cachedRMatrix(rtti(), 0);
    this->dNdrst(RVector3(0.0, 0.0, 0.0), MdNdrst);
//     RMatrix MdNdrst(this->dNdrst(RVector3(0.0, 0.0, 0.0)));

    switch (this->dim()){
        case 1:{
            J.setVal(MdNdrst * x, 0);
            RVector J0(J.row(0));
            RVector J1(J.row(1));

            if ((abs(J0[0]) >= abs(J0[1]) && abs(J0[0]) >= abs(J0[2])) ||
                (abs(J0[1]) >= abs(J0[0]) && abs(J0[1]) >= abs(J0[2]))) {
                J1[0] =  J0[1];
                J1[1] = -J0[0];
                J1[2] =     0.;
            } else {
                J1[0] =     0.;
                J1[1] =  J0[2];
                J1[2] = -J0[1];
            }

            J.setVal(J1 / norml2(J1), 1);
            J.setVal(crossN(J.row(0), J.row(1)), 2);
// __M
//             std::cout << J1 << std::endl;
//             std::cout << J << std::endl;
         } break;
        case 2:
            J.setVal(MdNdrst * x, 0);
            J.setVal(MdNdrst * y, 1);
            J.setVal(crossN(J.row(0), J.row(1)), 2);
            break;
        case 3:
            J.setVal(MdNdrst * x, 0);
            J.setVal(MdNdrst * y, 1);
            J.setVal(MdNdrst * z, 2);
            break;
    }
}


RVector3 Shape::xyz(const RVector3 & rst) const{
    RVector3 xyz(0.0, 0.0, 0.0);
    rst2xyz(rst, xyz);
    return xyz;
}

void Shape::rst2xyz(const RVector3 & rst, RVector3 & xyz) const{
    RVector sf(this->N(rst));

    for (Index i = 0; i < nodeCount(); i ++){
        xyz += nodeVector_[i]->pos() * sf[i];
    }
}

const RMatrix3 & Shape::invJacobian() const {
    if (!invJacobian_.valid()){
       if (ShapeFunctionCache::instance().RMatrix3Cache().size() < 1){
           ShapeFunctionCache::instance().RMatrix3Cache().push_back(RMatrix3());
       }

       this->createJacobian(ShapeFunctionCache::instance().cachedRMatrix3(0));
       inv(ShapeFunctionCache::instance().cachedRMatrix3(0), invJacobian_);
       invJacobian_.setValid(true);
    }
//     if (invJacobian_.rows() != 3) {
//         invJacobian_.resize(3,3);
// //RMatrix J(3,3);
//         this->createJacobian(__tmp__J__);
//         inv(__tmp__J__, invJacobian_ );
//     }
    return invJacobian_;
}

void Shape::xyz2rst(const RVector3 & xyz, RVector3 & rst) const{

    double err = 1., tol = 1e-10;
    uint maxiter = 200;
    uint iter = 0;
    RVector3 dxyz(0.0, 0.0, 0.0);
    RVector3 drst(0.0, 0.0, 0.0);

    double dErr = 10.0;
//     double lastdErr = 0.0;
    double damping = 1.0;

    while (abs(dErr) > tol && iter < maxiter){
        if (err > 1000 && iter > 1){
            damping *= 0.9;
            rst = RVector3(0.0, 0.0, 0.0);
            err = 1;
            iter = 0;
        }

        iter ++;
        dxyz = xyz - this->xyz(rst);
        drst = invJacobian() * dxyz;

        rst += drst * damping;
//         lastdErr = dErr;
        dErr = err - drst.abs();
        err = drst.abs();

//        __MS(iter << " " << err << " " << damping << " " << dErr)
    }

//     if (err> 1){
//         __MS(iter << " " << err << " " << damping << " " << dErr)
//         for (Index i = 0; i < nodeCount(); i ++){
//             __MS(nodeVector_[i]->pos());
//         }
//             __MS(this->createJacobian().det())
//         std::cout << xyz << " " << rst << std::endl;
// exit(0);
//     }
}

RVector3 Shape::rst(const RVector3 & xyz) const{
    RVector3 rst(0.0, 0.0, 0.0);
    xyz2rst(xyz, rst);
    return rst;
}

RVector3 Shape::rst(Index i) const {
    std::cout << "shape: " << rtti() << std::endl;
    THROW_TO_IMPL

    return RVector3(0.0, 0.0, 0.0);
}

bool Shape::isInside(const RVector3 & xyz, bool verbose) const {
     RVector sf; return isInside(xyz, sf, verbose);
}

bool Shape::isInside(const RVector3 & xyz, RVector & sf, bool verbose) const {

    sf = N(rst(xyz));
    double minsf = min(sf);

    if (verbose){
        std::cout << "rst: " << rst(xyz)<< std::endl;
        std::cout << "sf: " << sf << std::endl;
        std::cout << std::fabs(minsf) << " " << max(TOUCH_TOLERANCE, TOUCH_TOLERANCE * xyz.abs()) << std::endl;
    }

    if (std::fabs(minsf) < max(TOUCH_TOLERANCE, TOUCH_TOLERANCE * xyz.abs())) return true; //** on boundary
    if (minsf > 0.0) return true; //** inside
    return false;
}

double Shape::domainSize() const {
    if (!hasDomSize_) {
        domSize_ = domainSize_();
        //** gefaehrlich nach nem node::trans, rot, scale stimmt das nicht mehr (siehe unittest:testFEM)
        hasDomSize_ = true;
    }
    return domSize_;
}

//** Start Node specific implementation

RVector3 NodeShape::norm() const {
    return RVector3(1.0, 0.0, 0.0);
}

RVector3 NodeShape::rst(Index i) const {
    return RVector3(0.0, 0.0, 0.0);
}


//** Start EDGE specific implementation
double EdgeShape::length() const {
    return nodeVector_[0]->dist(*nodeVector_[1]);
}

RVector3 EdgeShape::norm() const {
    return nodeVector_[0]->pos().normXY(nodeVector_[1]->pos());
}

RVector3 EdgeShape::rst(Index i) const{
    if (i < nodeCount())
        return RVector3(EdgeCoordinates[i][0],
                        EdgeCoordinates[i][1],
                        EdgeCoordinates[i][2]);
    THROW_TO_IMPL; return RVector3(0.0, 0.0, 0.0);
}

void TriangleShape::setNodes(Node * n0, Node * n1, Node * n2){
    setNode(0, *n0); setNode(1, *n1); setNode(2, *n2);
}

//** Start TRIANGLE specific implementation
double TriangleShape::area() const {
//   RVector3 p1 = nodeVector_[0]->pos();
//   RVector3 p2 = nodeVector_[1]->pos();
//   RVector3 p3 = nodeVector_[2]->pos();

    RVector3 a(nodeVector_[1]->pos() - nodeVector_[0]->pos());
    RVector3 b(nodeVector_[2]->pos() - nodeVector_[0]->pos());

    return ((a).cross(b)).abs() * 0.5;

//   return (nodeVector_[1]->pos() - nodeVector_[0]->pos())
//     .cross(nodeVector_[2]->pos() - nodeVector_[0]->pos()).abs() * 0.5;
}

RVector3 TriangleShape::norm() const{
    RVector3 a(nodeVector_[1]->pos() - nodeVector_[0]->pos());
    RVector3 b(nodeVector_[2]->pos() - nodeVector_[0]->pos());
    RVector3 n((a).cross(b));
    return n.norm();
}

RVector3 TriangleShape::rst(Index i) const{
    if (i < nodeCount()) return RVector3(TriCoordinates[i][0], TriCoordinates[i][1], TriCoordinates[i][2]);
    THROW_TO_IMPL; return RVector3(0.0, 0.0, 0.0);
}

void TriangleShape::xyz2rst(const RVector3 & pos, RVector3 & rst ) const {
    //return Shape::xyz2rst(pos, rst);

 //** Coordinate transformation
//**     xp = x1 + (x2 - x1) * r + (x3 - x1) * s
//**     yp = y1 + (y2 - y1) * r + (y3 - y1) * s
//**     xp1 = x21 * r + x31 * s
//**     yp1 = y21 * r + y31 * s
//** maxima: eqns: [xp1 = x21 * r + x31 * s, yp1 = y21 * r + y31 * s]
//** maxima: linsolve(eqns, [r,s])
//**    [r=-(xp1*y31-x31*yp1)/(x31*y21-x21*y31),s=(xp1*y21-x21*yp1)/(x31*y21-x21*y31)]
//**     J = x21 * y31 - x31 * y21

    double x21 = nodeVector_[1]->pos()[0] - nodeVector_[0]->pos()[0];
    double x31 = nodeVector_[2]->pos()[0] - nodeVector_[0]->pos()[0];
    double y21 = nodeVector_[1]->pos()[1] - nodeVector_[0]->pos()[1];
    double y31 = nodeVector_[2]->pos()[1] - nodeVector_[0]->pos()[1];
    double xp1 = pos[0] - nodeVector_[0]->pos()[0];
    double yp1 = pos[1] - nodeVector_[0]->pos()[1];

    //** use here the local J instead of jacobianDeterminant_(), while it use the area-hack
    double J = x21 * y31 - x31 * y21;

    rst[0] = (y31 * xp1 - x31 * yp1) / J; // r
    rst[1] = (x21 * yp1 - y21 * xp1) / J; // s
}

RVector3 QuadrangleShape::rst(Index i) const{
    if (i < nodeCount()) return RVector3(QuadCoordinates[i][0], QuadCoordinates[i][1], QuadCoordinates[i][2]);
    THROW_TO_IMPL; return RVector3(0.0, 0.0, 0.0);
}

double QuadrangleShape::area() const {
    double sum = 0.0;
    TriangleShape tri;
    tri.setNodes(nodeVector_[0], nodeVector_[1], nodeVector_[2]);
    sum += tri.area();
    tri.setNodes(nodeVector_[0], nodeVector_[2], nodeVector_[3]);
    sum += tri.area();

    return sum;

    //** Gau sche Trapezformel fails for surface boundaries

/*
  double x13 = nodeVector_[0]->pos()[0]- nodeVector_[2]->pos()[0];
  double x42 = nodeVector_[3]->pos()[0]- nodeVector_[1]->pos()[0];

  double y13 = nodeVector_[0]->pos()[1]- nodeVector_[2]->pos()[1];
  double y24 = nodeVector_[1]->pos()[1]- nodeVector_[3]->pos()[1];
  return 0.5 * std::fabs(y13 * x42 + y24 * x13);*/
    //return std::fabs(this->jacobianDeterminant());
}

RVector3 QuadrangleShape::norm() const {
    TriangleShape tri;
    tri.setNodes(nodeVector_[0], nodeVector_[1], nodeVector_[2]);
    return tri.norm();
}

RVector3 TetrahedronShape::rst(Index i) const{
    if (i < nodeCount()) return RVector3(TetCoordinates[i][0], TetCoordinates[i][1], TetCoordinates[i][2]);
    THROW_TO_IMPL; return RVector3(0.0, 0.0, 0.0);
}

void TetrahedronShape::xyz2rst(const RVector3 & pos, RVector3 & rst ) const {
    //return Shape::xyz2rst(pos, rst);

//**     xp = x1 + (x2 - x1) * r + (x3 - x1) * s + (x4 - x1) * t
//**     yp = y1 + (y2 - y1) * r + (y3 - y1) * s + (y4 - y1) * t
//**     zp = z1 + (z2 - z1) * r + (z3 - z1) * s + (z4 - z1) * t
//**     xp1 = x21 * r + x31 * s + x41 * t
//**     yp1 = y21 * r + y31 * s + y41 * t
//**     zp1 = z21 * r + z31 * s + z41 * t
//** maxima: eqns: [xp1 = x21 * r + x31 * s + x41 * t, yp1 = y21 * r + y31 * s + y41 * t, zp1 = z21 * r + z31 * s + z41 * t]
//** maxima: linsolve(eqns, [r,s,t])
//     r=(x31*(yp1*z41-y41*zp1)+x41*(y31*zp1-yp1*z31)+xp1*(y41*z31-y31*z41))/
//             (x21*(y41*z31-y31*z41)+x31*(y21*z41-y41*z21)+x41*(y31*z21-y21*z31)),
//     s=-(x21*(yp1*z41-y41*zp1)+x41*(y21*zp1-yp1*z21)+xp1*(y41*z21-y21*z41))/
//             (x21*(y41*z31-y31*z41)+x31*(y21*z41-y41*z21)+x41*(y31*z21-y21*z31)),
//     t=(x21*(yp1*z31-y31*zp1)+x31*(y21*zp1-yp1*z21)+xp1*(y31*z21-y21*z31))/
//             (x21*(y41*z31-y31*z41)+x31*(y21*z41-y41*z21)+x41*(y31*z21-y21*z31))]
//               Jac = -(x21*(y41*z31-y31*z41)+x31*(y21*z41-y41*z21)+x41*(y31*z21-y21*z31))
    double x21 = nodeVector_[1]->pos()[0] - nodeVector_[0]->pos()[0];
    double x31 = nodeVector_[2]->pos()[0] - nodeVector_[0]->pos()[0];
    double x41 = nodeVector_[3]->pos()[0] - nodeVector_[0]->pos()[0];
    double y21 = nodeVector_[1]->pos()[1] - nodeVector_[0]->pos()[1];
    double y31 = nodeVector_[2]->pos()[1] - nodeVector_[0]->pos()[1];
    double y41 = nodeVector_[3]->pos()[1] - nodeVector_[0]->pos()[1];
    double z21 = nodeVector_[1]->pos()[2] - nodeVector_[0]->pos()[2];
    double z31 = nodeVector_[2]->pos()[2] - nodeVector_[0]->pos()[2];
    double z41 = nodeVector_[3]->pos()[2] - nodeVector_[0]->pos()[2];
    double xp1 = pos[0] - nodeVector_[0]->pos()[0];
    double yp1 = pos[1] - nodeVector_[0]->pos()[1];
    double zp1 = pos[2] - nodeVector_[0]->pos()[2];

    double J = x21 * (y31 * z41 - y41 * z31) +
               x31 * (y41 * z21 - y21 * z41) +
               x41 * (y21 * z31 - y31 * z21);

    rst[0] = (xp1 * (y31 * z41 - y41 * z31) +
                    x31 * (y41 * zp1 - yp1 * z41) +
                    x41 * (yp1 * z31 - y31 * zp1)) / J;

    rst[1] = (x21 * (yp1 * z41 - y41 * zp1) +
                    xp1 * (y41 * z21 - y21 * z41) +
                    x41 * (y21 * zp1 - yp1 * z21)) / J;

    rst[2] = (x21 * (y31 * zp1 - yp1 * z31) +
                    x31 * (yp1 * z21 - y21 * zp1) +
                    xp1 * (y21 * z31 - y31 * z21)) / J;
}

double TetrahedronShape::volume() const {
    RVector3 a(nodeVector_[1]->pos() - nodeVector_[0]->pos());
    RVector3 b(nodeVector_[2]->pos() - nodeVector_[0]->pos());
    RVector3 c(nodeVector_[3]->pos() - nodeVector_[0]->pos());
//     std::cout << a << " " << b << " " << c << " " << a.cross(b) << std::endl;
//     std::cout << std::fabs((a.cross(b).dot(c))) << std::endl;
    return 1.0 / 6.0 * std::fabs((a.cross(b).dot(c)));
    //** pls check whish way is faster, profile, with valgrind and count ticks
//     return fabs(jacobianDeterminant() / 6.0);
}

void TetrahedronShape::setNodes(Node * n0, Node * n1, Node * n2, Node * n3){
    setNode(0, *n0); setNode(1, *n1); setNode(2, *n2); setNode(3, *n3);
}

RVector3 HexahedronShape::rst(Index i) const{
    if (i < nodeCount()) return RVector3(HexCoordinates[i][0], HexCoordinates[i][1], HexCoordinates[i][2]);
    THROW_TO_IMPL; return RVector3(0.0, 0.0, 0.0);
}

double HexahedronShape::volume() const {
    double sum = 0.0;
    TetrahedronShape tet;
    for (uint i = 0; i < 5; i ++){
        tet.setNodes(nodeVector_[HexahedronSplit5TetID[i][0]], nodeVector_[HexahedronSplit5TetID[i][1]],
                      nodeVector_[HexahedronSplit5TetID[i][2]], nodeVector_[HexahedronSplit5TetID[i][3]]);
        sum += tet.volume();
    }

    return sum;
}

bool HexahedronShape::enforcePositiveDirection(){
    if (createJacobian().det() < 0){
//         std::reverse(nodeVector_.begin(), nodeVector_.end());

        //swap up and down side
        if (createJacobian().det() < 0){
            std::swap(nodeVector_[0], nodeVector_[4]);
            std::swap(nodeVector_[1], nodeVector_[5]);
            std::swap(nodeVector_[2], nodeVector_[6]);
            std::swap(nodeVector_[3], nodeVector_[7]);
        }
        this->changed();
        return true;
    }
    return false;
}

RVector3 TriPrismShape::rst(Index i) const{
    if (i < nodeCount()) return RVector3(PrismCoordinates[i][0], PrismCoordinates[i][1], PrismCoordinates[i][2]);
    THROW_TO_IMPL;
    return RVector3(0.0, 0.0, 0.0);
}

double TriPrismShape::volume() const{
    double sum = 0.0;
    TetrahedronShape tet;
    for (Index i = 0; i < 3; i ++){
        tet.setNodes(nodeVector_[TriPrimSplit3TetID[i][0]], nodeVector_[TriPrimSplit3TetID[i][1]],
                      nodeVector_[TriPrimSplit3TetID[i][2]], nodeVector_[TriPrimSplit3TetID[i][3]]);
        sum += tet.volume();
    }

    return sum;
    /*
    RVector3 a(nodeVector_[1]->pos() - nodeVector_[0]->pos());
    RVector3 b(nodeVector_[2]->pos() - nodeVector_[0]->pos());
    double z41 = nodeVector_[3]->pos()[2] - nodeVector_[0]->pos()[2];

    return ((a).cross(b)).abs()  z41;*/
}

std::vector < PolynomialFunction < double > > TriPrismShape::createShapeFunctions() const {
    RVector e2(2); e2[0] = 0; e2[1] =  1; // x

    RPolynomialFunction T3_2(e2, RVector(0));
    RPolynomialFunction T3_3(RVector(0), e2);
    RPolynomialFunction T3_1 = -(-1.0 + T3_2 + T3_3);

    RPolynomialFunction E2_2T(RVector(0), RVector(0), e2);
    RPolynomialFunction E2_1T = -(-1.0 + E2_2T);
    std::vector < PolynomialFunction < double > > ret;

    ret.push_back(T3_1 * E2_1T);
    ret.push_back(T3_2 * E2_1T);
    ret.push_back(T3_3 * E2_1T);
    ret.push_back(T3_1 * E2_2T);
    ret.push_back(T3_2 * E2_2T);
    ret.push_back(T3_3 * E2_2T);
    return ret;
}

RVector3 PyramidShape::rst(Index i) const{
    if (i < nodeCount()) return RVector3(PyramidCoordinates[i][0],
                                                      PyramidCoordinates[i][1],
                                                      PyramidCoordinates[i][2]);
    THROW_TO_IMPL; return RVector3(0.0, 0.0, 0.0);
}

} // namespace GIMLI

