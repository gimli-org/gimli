#ifndef PYTHON_PYGIMLI__H
#define PYTHON_PYGIMLI__H

//See best practices section in Py++ documentation

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <list>
#include <map>

// #ifndef PACKAGE_NAME
//         #define PACKAGE_NAME "gimli (python build)"
//         #define PACKAGE_VERSION "0.9.0"
//         #define PACKAGE_BUGREPORT "carsten@resistivity.net"
// #endif // PACKAGE_NAME

#ifdef PYTEST

#include "stopwatch.h"
#include "vector.h"

namespace GIMLI{

	template class Vector< double >;
    template class Vector< Complex >;
    template class Vector< int >;

    inline void ___instantiation___(){
        sizeof(::GIMLI::Index *);
		sizeof(::GIMLI::Index &);
		sizeof(::GIMLI::Index);
		sizeof(::GIMLI::IndexArray *);
		sizeof(::GIMLI::IndexArray &);
		sizeof(::GIMLI::IndexArray);
		sizeof(int *);
        sizeof(int &);
        sizeof(int);
        sizeof(long unsigned int *);
        sizeof(long unsigned int &);
        sizeof(long unsigned int);
		sizeof(long long unsigned int *);
        sizeof(long long unsigned int &);
        sizeof(long long unsigned int);
		sizeof(long long int *);
        sizeof(long long int &);
        sizeof(long long int);
        sizeof(double *);
        sizeof(double);
        sizeof(double &);
    }

} // namespace GIMLI

namespace pyplusplus{ namespace aliases{
    //typedef std::complex< double >                  Complex;

    //typedef GIMLI::Vector< double >                 RVector;
  //  typedef GIMLI::Vector< Complex >                CVector;
//    typedef GIMLI::Pos< double >                    RVector3;

//    typedef GIMLI::Matrix< double >                 RMatrix;

    //typedef std::vector< double >                   stdVectorR;
    typedef std::vector< int >                        stdVectorI;
    //typedef std::vector< std::complex< double > >   stdVectorC;

    //typedef std::vector< std::string > stdVectorString;
    //typedef std::map< float, float > stdMapF_F;
}
} //pyplusplus::aliases

#else // if not PYTEST

#include "baseentity.h"
#include "curvefitting.h"
#include "cholmodWrapper.h"
#include "datacontainer.h"
#include "gimli.h"
#include "gravimetry.h"
#include "dc1dmodelling.h"
#include "elementmatrix.h"
#include "em1dmodelling.h"
#include "exitcodes.h"
#include "expressions.h"
#include "interpolate.h"
#include "integration.h"
#include "inversion.h"
#include "ipcClient.h"
#include "ldlWrapper.h"
#include "line.h"
#include "linSolver.h"
#include "matrix.h"
#include "mesh.h"
#include "meshentities.h"
#include "meshgenerators.h"
#include "modellingbase.h"
#include "node.h"
#include "numericbase.h"
#include "plane.h"
#include "platform.h"
#include "pos.h"
#include "polynomial.h"
#include "quaternion.h"
#include "regionManager.h"
#include "shape.h"
#include "solver.h"
#include "solverWrapper.h"
#include "sparsematrix.h"
#include "spline.h"
#include "stopwatch.h"
#include "trans.h"
#include "triangleWrapper.h"
#include "ttdijkstramodelling.h"
#include "vector.h"
#include "vectortemplates.h"

// #include "inversionRollalong.h"
// #include "eameshwrapper.h"
// #include "matrixTemplates.h"

namespace GIMLI{
#define DEFINE_PY_VEC_OPERATOR__(OP)                      \
    inline RVector operator OP (const RVector & a, const RVector & b){ \
			RVector c(a);	c OP##=b; return c; } 				\
    inline RVector operator OP (const RVector & a, double v){ \
			RVector b(a);	b OP##=v; return b; } 				\
    inline RVector operator OP (double v, const RVector & a){ \
			RVector b(a);	for (uint i = 0; i < b.size(); i ++) b[ i ] = v OP b[i]; return b; } \

    DEFINE_PY_VEC_OPERATOR__(+)
    DEFINE_PY_VEC_OPERATOR__(-)
    DEFINE_PY_VEC_OPERATOR__(*)
    DEFINE_PY_VEC_OPERATOR__(/)
#undef DEFINE_PY_VEC_OPERATOR__

#define DEFINE_PY_VEC_UNARY_OPERATOR__(OP, FUNCT)                      \
    inline RVector OP (const RVector & a) { \
        RVector tmp(a.size()); for (uint i = 0; i < a.size(); i ++) tmp[ i ] = FUNCT()(a[ i ]); return tmp; } \

DEFINE_PY_VEC_UNARY_OPERATOR__(abs,   ABS_)
DEFINE_PY_VEC_UNARY_OPERATOR__(acot,  ACOT)
DEFINE_PY_VEC_UNARY_OPERATOR__(atan,  ATAN)
DEFINE_PY_VEC_UNARY_OPERATOR__(cos,   COS)
DEFINE_PY_VEC_UNARY_OPERATOR__(cot,   COT)
DEFINE_PY_VEC_UNARY_OPERATOR__(exp,   EXP)
DEFINE_PY_VEC_UNARY_OPERATOR__(exp10, EXP10)
DEFINE_PY_VEC_UNARY_OPERATOR__(fabs,  ABS_)
DEFINE_PY_VEC_UNARY_OPERATOR__(log,   LOG)
DEFINE_PY_VEC_UNARY_OPERATOR__(log10, LOG10)
DEFINE_PY_VEC_UNARY_OPERATOR__(sign,  SIGN)
DEFINE_PY_VEC_UNARY_OPERATOR__(sin,   SIN)
DEFINE_PY_VEC_UNARY_OPERATOR__(sqrt,  SQRT)
DEFINE_PY_VEC_UNARY_OPERATOR__(square, SQR)
DEFINE_PY_VEC_UNARY_OPERATOR__(tan,   TAN)
DEFINE_PY_VEC_UNARY_OPERATOR__(tanh,  TANH)
#undef DEFINE_PY_VEC_UNARY_OPERATOR__

#define DEFINE_COMPARE_OPERATOR__(OP) \
template BVector operator OP (const std::vector < int > & vec, const int & v); \
template BVector operator OP (const std::vector < uint > & vec, const uint & v); \

DEFINE_COMPARE_OPERATOR__(<)
DEFINE_COMPARE_OPERATOR__(<=)
DEFINE_COMPARE_OPERATOR__(>=)
DEFINE_COMPARE_OPERATOR__(==)
DEFINE_COMPARE_OPERATOR__(!=)
DEFINE_COMPARE_OPERATOR__(>)

#undef DEFINE_COMPARE_OPERATOR__


    //** maybe better to move these instantiation into libgimli

    template class Vector< double >;
    template class Vector< Complex >;
    template class Vector< int >;

    template class Matrix<  double >;
   // template class Matrix< Complex >;
    
    template RVector unique(const RVector & a);
    template RVector sort(const RVector & a);
    template RVector pow(const RVector & a, double power);
    template RVector pow(const RVector & a, int power);

    template class Pos< double >;
    template class Quaternion< double >;
    
    template class PolynomialElement< double >;
    template class PolynomialFunction< double >;
    template std::ostream & operator << (std::ostream & str, const PolynomialFunction < double > & p);
    template PolynomialFunction < double > operator - (const PolynomialFunction < double > & f);
    template PolynomialFunction < double > operator - (const PolynomialFunction < double > & f, const PolynomialFunction < double > & g);
    template PolynomialFunction < double > operator + (const PolynomialFunction < double > & f, const PolynomialFunction < double > & g);
    template PolynomialFunction < double > operator * (const PolynomialFunction < double > & f, const PolynomialFunction < double > & g);
    template PolynomialFunction < double > operator * (const PolynomialFunction < double > & f, const double & val);
    template PolynomialFunction < double > operator * (const double & val, const PolynomialFunction < double > & f);
    template PolynomialFunction < double > operator + (const PolynomialFunction < double > & f, const double & val);
    template PolynomialFunction < double > operator + (const double & val, const PolynomialFunction < double > & f);
    
    template std::vector < PolynomialFunction < double > > 
        createPolynomialShapeFunctions(const Shape & ent, uint nCoeff, bool pascale, bool serendipity, const RVector &);
    template std::vector < PolynomialFunction < double > > 
        createPolynomialShapeFunctions(const MeshEntity & ent, uint nCoeff, bool pascale, bool serendipity, const RVector &);
        
    template class Inversion< double >;

    template class Trans< RVector >;
    template class TransLinear< RVector >;
    template class TransPower< RVector >;
    template class TransLog< RVector >;
    template class TransLogLU< RVector >;
    template class TransCotLU< RVector >;

    template std::vector< Index > unique(const std::vector < Index > & a);
    template std::vector< SIndex > unique(const std::vector < SIndex > & a);
    
    template std::vector< Index > sort(const std::vector < Index > & a);
    template std::vector< SIndex > sort(const std::vector < SIndex > & a);
    
    template RVector fliplr(const RVector & a);

    template RVector increasingRange(const double & first, const double & last, Index n);

    template Mesh & Mesh::transform(const Matrix < double > & mat);
    template Pos< double > & Pos< double >::transform(const Matrix < double > & mat);

    template class SparseMatrix< double >;
    template Vector < double > operator * (const SparseMatrix < double > & A, const Vector < double > & a);
    template SparseMatrix< double > operator + (const SparseMatrix < double > & A, const SparseMatrix< double > & B);
    template SparseMatrix< double > operator * (const SparseMatrix < double > & A, const double & b);
    template SparseMatrix< double > operator * (const double & b, const SparseMatrix < double > & A);

   
    template class ElementMatrix< double >;
    template std::ostream & operator << (std::ostream & str, const ElementMatrix< double > & p);

    template RVector cat(const RVector & a, const RVector & b);

    template double sum(const RVector & v);
    template double min(const RVector & v);
    template double max(const RVector & v);
    //template Complex sum(const CVector & v);
    template double rms(const RVector & a);
    template double rms(const RVector & a, const RVector & b);
    template double rrms(const RVector & a, const RVector & b);
    template double norm(const RVector & a);
    template double normlp(const RVector & a, int p);
    template double norml1(const RVector & a);
    template double norml2(const RVector & a);
    template double normlInfinity(const RVector & a);
    template double euclideanNorm(const RVector & a);
    template double stdDev(const RVector & v);
    template double median(const RVector & a);
    template double mean(const RVector & a);
    template double arithmeticMean(const RVector & v);
    template double geometricMean(const RVector & v);
    template double harmonicMean(const RVector & a);

    template void rand(RVector & vec, double min = 0.0, double max = 1.0);
    template void randn(RVector & vec);

    template double degToRad(const double & deg);
    template double radToDeg(const double & rad);

    template RVector3 degToRad(const RVector3 & deg);
    template RVector3 radToDeg(const RVector3 & rad);

    template bool save(const RVector &v, const std::string & fname, IOFormat format = Ascii);
    template bool load(RVector &v, const std::string & fname, IOFormat format = Ascii,
                        bool verbose = true);

  //  template bool save(const RMatrix & A, const std::string & filename, IOFormat format = Binary);
    template bool load(RMatrix & A, const std::string & filename);

    template bool loadMatrixSingleBin(RMatrix & A, const std::string & filename);
    template bool loadMatrixVectorsBin(RMatrix & A, const std::string & filenameBody,
                                        uint kCount = 1);

    template bool saveMatrixCol(const RMatrix & A, const std::string & filename);
    template bool saveMatrixCol(const RMatrix & A, const std::string & filename,
                                 const std::string & comments);

    template bool loadMatrixCol(RMatrix & A, const std::string & filename);
    template bool loadMatrixCol(RMatrix & A, const std::string & filename,
                                 std::vector < std::string > & comments);

    template bool saveMatrixRow(const RMatrix & A, const std::string & filename);
    template bool saveMatrixRow(const RMatrix & A, const std::string & filename,
                                 const std::string & comments);

    template bool loadMatrixRow(RMatrix & A, const std::string & filename);
    template bool loadMatrixRow(RMatrix & A, const std::string & filename,
                                std::vector < std::string > & comments);

    template std::set< Node * > commonNodes(const std::set < Boundary * > &);
    template std::set< Node * > commonNodes(const std::set < Cell * > &);

    inline void ___instantiation___(){
        sizeof(bool);
        sizeof(bool *);
        sizeof(char);
        sizeof(char *);
        sizeof(char &);
        sizeof(false);
        sizeof(true);
        sizeof(int *);
        sizeof(int &);
        sizeof(int);
        sizeof(uint8 *);
        sizeof(uint8 &);
        sizeof(uint8);
        sizeof(unsigned char);
        sizeof(long unsigned int *);
        sizeof(long unsigned int &);
        sizeof(long unsigned int);
		sizeof(long long unsigned int *);
		sizeof(long long unsigned int &);
		sizeof(long long unsigned int);
        sizeof(double *);
        sizeof(double);
        sizeof(double &);
        sizeof(GIMLI::SparseMatrix<double> &);
        sizeof(GIMLI::SparseMatrix<double> *);
        sizeof(GIMLI::SparseMatrix<double>);
		sizeof(::GIMLI::Index *);
		sizeof(::GIMLI::Index &);
		sizeof(::GIMLI::Index);
		sizeof(::GIMLI::IndexArray *);
		sizeof(::GIMLI::IndexArray &);
		sizeof(::GIMLI::IndexArray);
    }

    /*! Temporary workaround until there is as solution for
        http://www.mail-archive.com/cplusplus-sig@python.org/msg00333.html/

        declare and init new function, inject them to original class within __init__.py
    */
    // template void Quaternion< double >::rotMatrix(Matrix < double > & mat) const;
    inline void getRotMatrix__(const Quaternion < double > & q, Matrix < double > & mat){
        q.rotMatrix(mat);
    }

    inline void interpolate_GILsave__(const Mesh & mesh, const RMatrix & data,
                                       const std::vector< RVector3 > & pos, RMatrix & iData,
                                       bool verbose = false){
        std::cout << "interpolate_GILsave__ 1" << std::endl;
        ALLOW_PYTHON_THREADS
        return interpolate(mesh, data, pos, iData, verbose);
    }
    inline RVector interpolate_GILsave__(const Mesh & mesh, const RVector & data,
                                          const RVector & x, const RVector & y, bool verbose = false){
        ALLOW_PYTHON_THREADS
        std::cout << "interpolate_GILsave__ 2" << std::endl;
        return interpolate(mesh, data, x, y, verbose);
    }

    inline void interpolate_GILsave__(const Mesh & mesh, const RVector & data,
                                       const Mesh & pos, RVector & iData, bool verbose = false){
        ALLOW_PYTHON_THREADS
        std::cout << "interpolate_GILsave__ 3" << std::endl;
        return interpolate(mesh, data, pos, iData, verbose);
    }

} // namespace GIMLI

//** define some aliases to avoid insane method names
namespace pyplusplus{ namespace aliases{
    typedef std::complex< double >                       Complex;

    typedef GIMLI::Vector< double >                      RVector;
    typedef GIMLI::Vector< Complex >                     CVector;
    typedef GIMLI::Vector< int >                         BVector;

    typedef GIMLI::Matrix< double >                      RMatrix;
    //typedef GIMLI::Matrix< Complex >                     CMatrix;

    typedef GIMLI::VectorIterator< double >              RVectorIter;
    typedef GIMLI::Pos< double >                         RVector3;

    typedef GIMLI::Quaternion< double >                  RQuaternion;
    
    typedef GIMLI::PolynomialFunction< double >          RPolynomialFunction;
    typedef GIMLI::PolynomialElement< double >           RPolynomialElement;

    typedef GIMLI::Inversion< double >                   RInversion;
//     typedef GIMLI::RollalongInSpace< double >            RRollalongInSpace;

    typedef GIMLI::ElementMatrix< double >                   DElementMatrix;
    typedef GIMLI::SparseMapMatrix< int, GIMLI::Index >      ISparseMapMatrix;
    typedef GIMLI::SparseMapMatrix< double, GIMLI::Index >   DSparseMapMatrix;
    typedef GIMLI::SparseMatrix< int >                       ISparseMatrix;
    typedef GIMLI::SparseMatrix< double >                    DSparseMatrix;

    typedef GIMLI::Trans< GIMLI::RVector >        RTrans;
    typedef GIMLI::TransLinear< GIMLI::RVector >  RTransLinear;
    typedef GIMLI::TransPower< GIMLI::RVector >   RTransPower;
    typedef GIMLI::TransLog< GIMLI::RVector >     RTransLog;
    typedef GIMLI::TransLogLU< GIMLI::RVector >   RTransLogLU;
    typedef GIMLI::TransCotLU< GIMLI::RVector >   RTransCotLU;
    typedef GIMLI::CumulativeTrans < GIMLI::RVector > RTransCumulative;

    typedef GIMLI::BoundingBox<double> RBoundingBox;

    typedef std::vector< GIMLI::RVector3 >           stdVectorRVector3;
    typedef std::vector< GIMLI::Vector< std::complex< double > > > stdVectorCVector;

    typedef std::vector< GIMLI::Boundary * >         stdVectorBounds;
    typedef std::vector< GIMLI::Cell * >             stdVectorCells;
    typedef std::vector< GIMLI::Node * >             stdVectorNodes;
//     typedef std::vector< GIMLI::ElectrodeShape * >   stdVectorElectrodeShape;
//     typedef std::vector< GIMLI::Electrode >          stdVectorElectrode;
//     typedef std::vector< GIMLI::Electrode  * >       stdVectorpElectrode;
    typedef std::vector< GIMLI::MeshEntity * >       stdVectorMeshEntity;
    typedef std::vector< GIMLI::CubicFunct >         stdVectorCubicFunct;

    typedef std::set< GIMLI::Node * >                 stdSetNodes;
    typedef std::set< GIMLI::Boundary * >             stdSetBoundary;
    typedef std::set< GIMLI::Cell * >                 stdSetCell;

    typedef std::map< std::string, GIMLI::Vector< double > > stdMapStringRVector;
    typedef std::map< int, GIMLI::Region * >            stdMapRegion;
    typedef std::map< float, float >                    stdMapF_F;
    typedef std::map< int, double >                     stdMapI_D;
    typedef std::map< long, long >                      stdMapL_L;
    typedef std::map< long, unsigned long >             stdMapL_UL;
    typedef std::map< unsigned long, double >           stdMapUL_D;

    typedef std::vector< std::string >                  stdVectorString;
    typedef std::vector< int >                          stdVectorI;
    typedef std::vector< GIMLI::Index >                 stdVectorUL;
    typedef std::vector< double >                       stdVectorR;
    typedef std::vector< std::complex < double > >      stdVectorC;
    typedef std::vector< std::pair< GIMLI::Index, GIMLI::Index > > stdVectorPairLongLong;
    typedef std::vector< GIMLI::Vector< double > >      stdVectorRVector;

    typedef std::vector< GIMLI::Trans < GIMLI::Vector < double > > * > stdVectorTrans;
    typedef std::vector< GIMLI::Vector< double > >      stdVectorRVector;
    typedef std::vector< GIMLI::Vector< std::complex< double > > > stdVectorCVector;

    typedef std::list< long unsigned int >              stdListUL;

    typedef std::set< long int >                        stdSetL;

}} //pyplusplus::aliases

#endif // else not PYTEST

#endif // PYTHON_PYGIMLI__H
