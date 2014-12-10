#ifndef _GIMLI_TESTEXTERNALS__H
#define _GIMLI_TESTEXTERNALS__H

#include <cppunit/extensions/HelperMacros.h>

#include <gimli.h>
#include <triangleWrapper.h>
#include <cholmodWrapper.h>
#include <sparsematrix.h>
#include <mesh.h>
#include <pos.h>

class GIMLIExternalTest : public CppUnit::TestFixture  {
    CPPUNIT_TEST_SUITE(GIMLIExternalTest);
    CPPUNIT_TEST(testTriangle);
    CPPUNIT_TEST(testCHOLMOD);
    CPPUNIT_TEST_SUITE_END();
    
public:    
    void testTriangle(){
        GIMLI::Mesh poly(2);
        GIMLI::Node *n1 = poly.createNode( GIMLI::RVector3( 0.0, 0.0 ) );
        GIMLI::Node *n2 = poly.createNode( GIMLI::RVector3( 1.0, 0.0 ) );
        GIMLI::Node *n3 = poly.createNode( GIMLI::RVector3( 1.0, 1.0 ) );
        poly.createEdge( *n1, *n2 );
        poly.createEdge( *n2, *n3 );
        poly.createEdge( *n3, *n1 );
        GIMLI::TriangleWrapper tri( poly );
        tri.setSwitches("-pzeQ");
        GIMLI::Mesh outMesh(2);
        tri.generate( outMesh );
        CPPUNIT_ASSERT( outMesh.cellCount() == 1 );
        CPPUNIT_ASSERT( outMesh.nodeCount() == 3 );
        CPPUNIT_ASSERT( outMesh.boundaryCount() == 3 );
    }
    
    template < class Matrix, class ValueType > 
        void testCHOLMODSolve(const Matrix & Sm){
        
        GIMLI::SparseMatrix< ValueType > S(Sm);
        GIMLI::CHOLMODWrapper solver(S, false);
        GIMLI::Vector < ValueType > b(S.rows(), ValueType(1));
        GIMLI::Vector < ValueType > x(S.rows());
        solver.solve(b, x);
        solver.solve(b, x);
        
//         std::cout << b - S * x << std::endl;
//         std::cout << GIMLI::Vector < ValueType >(b - S * x)<< std::endl;
//         std::cout << GIMLI::norm(b - S * x)<< std::endl;
        
        CPPUNIT_ASSERT(GIMLI::norm(b - S * x) < TOLERANCE);
        CPPUNIT_ASSERT(GIMLI::norm(b - Sm * x) < TOLERANCE);
    }
        
    void testCHOLMOD(){
        GIMLI::RSparseMapMatrix Sm(3,3);
        Sm.setVal(0,0,1.0);
        Sm.setVal(1,1,1.0);
        Sm.setVal(2,2,1.0);
        GIMLI::RSparseMatrix S(Sm);
        GIMLI::CHOLMODWrapper solver(S, true);

        GIMLI::RVector ones(Sm.rows(), 1.0);
        GIMLI::RVector x(Sm.rows());
        solver.solve(ones, x);
        CPPUNIT_ASSERT(GIMLI::norm(ones - S * x) == 0.0 );
        
        GIMLI::RSparseMapMatrix Rm(3, 3, 0);
        Rm.setVal(0, 0, 1.0);
        Rm.setVal(0, 1, 0.0);
        Rm.setVal(0, 2, 2.0);
        Rm.setVal(1, 0, 0.0);
        Rm.setVal(1, 1, 1.0);
        Rm.setVal(1, 2, 3.0);
        Rm.setVal(2, 0, 2.0);
        Rm.setVal(2, 1, 3.0);
        Rm.setVal(2, 2, 42.0);
        testCHOLMODSolve< GIMLI::RSparseMapMatrix, double>(Rm);
        
        GIMLI::RSparseMapMatrix RmU(3, 3, 1);
        RmU.setVal(0, 0, 1.0);
        RmU.setVal(0, 1, 0.0);
        RmU.setVal(0, 2, 2.0);
        RmU.setVal(1, 0, 0.0);
        RmU.setVal(1, 1, 1.0);
        RmU.setVal(1, 2, 3.0);
        RmU.setVal(2, 0, 2.0);
        RmU.setVal(2, 1, 3.0);
        RmU.setVal(2, 2, 42.0);
        testCHOLMODSolve< GIMLI::RSparseMapMatrix, double>(RmU);
        
        GIMLI::RSparseMapMatrix RmL(3, 3, -1);
        RmL.setVal(0, 0, 1.0);
        RmL.setVal(0, 1, 0.0);
        RmL.setVal(0, 2, 2.0);
        RmL.setVal(1, 0, 0.0);
        RmL.setVal(1, 1, 1.0);
        RmL.setVal(1, 2, 3.0);
        RmL.setVal(2, 0, 2.0);
        RmL.setVal(2, 1, 3.0);
        RmL.setVal(2, 2, 42.0);
        testCHOLMODSolve< GIMLI::RSparseMapMatrix, double>(RmL);
        
        //hermetian
        GIMLI::CSparseMapMatrix SmC(3, 3, 0);
        SmC.setVal(0, 0, GIMLI::Complex( 1.0,  0.0));
        SmC.setVal(0, 1, GIMLI::Complex( 1e-5, 0.0));
        SmC.setVal(0, 2, GIMLI::Complex( 2.0,  1.0));
        SmC.setVal(1, 0, GIMLI::Complex( 1e-5, 0.0));
        SmC.setVal(1, 1, GIMLI::Complex( 1.0,  0.0));
        SmC.setVal(1, 2, GIMLI::Complex( 3.0,  0.0));
        SmC.setVal(2, 0, GIMLI::Complex( 2.0,  -1.0));
        SmC.setVal(2, 1, GIMLI::Complex( 3.0,  0.0));
        SmC.setVal(2, 2, GIMLI::Complex(42.0,  0.0));
        testCHOLMODSolve< GIMLI::CSparseMapMatrix, GIMLI::Complex >(SmC);

        GIMLI::CSparseMapMatrix SmCL(3, 3, -1);
        SmCL.setVal(0, 0, GIMLI::Complex( 1.0,  0.0));
        SmCL.setVal(0, 1, GIMLI::Complex( 1e-5, 0.0));
        SmCL.setVal(0, 2, GIMLI::Complex( 2.0,  1.0));
        SmCL.setVal(1, 0, GIMLI::Complex( 1e-5, 0.0));
        SmCL.setVal(1, 1, GIMLI::Complex( 1.0,  0.0));
        SmCL.setVal(1, 2, GIMLI::Complex( 3.0,  0.0));
        SmCL.setVal(2, 0, GIMLI::Complex( 2.0,  -1.0));
        SmCL.setVal(2, 1, GIMLI::Complex( 3.0,  0.0));
        SmCL.setVal(2, 2, GIMLI::Complex(42.0,  0.0));
        testCHOLMODSolve< GIMLI::CSparseMapMatrix, GIMLI::Complex >(SmCL);

        GIMLI::CSparseMapMatrix SmCU(3, 3, 1);
        SmCU.setVal(0, 0, GIMLI::Complex( 1.0,  0.0));
        SmCU.setVal(0, 1, GIMLI::Complex( 1e-5, 0.0));
        SmCU.setVal(0, 2, GIMLI::Complex( 2.0,  1.0));
        SmCU.setVal(1, 0, GIMLI::Complex( 1e-5, 0.0));
        SmCU.setVal(1, 1, GIMLI::Complex( 1.0,  0.0));
        SmCU.setVal(1, 2, GIMLI::Complex( 3.0,  0.0));
        SmCU.setVal(2, 0, GIMLI::Complex( 2.0,  -1.0));
        SmCU.setVal(2, 1, GIMLI::Complex( 3.0,  0.0));
        SmCU.setVal(2, 2, GIMLI::Complex(42.0,  0.0));
        testCHOLMODSolve< GIMLI::CSparseMapMatrix, GIMLI::Complex >(SmCU);

        //non-hermetian symmetric
        GIMLI::CSparseMapMatrix SmCN(3, 3, 0);
        SmCN.setVal(0, 0, GIMLI::Complex( 1.0,  0.0));
        SmCN.setVal(0, 1, GIMLI::Complex( 1e-5, 0.0));
        SmCN.setVal(0, 2, GIMLI::Complex( -2.0,  1.0));
        SmCN.setVal(1, 0, GIMLI::Complex( 1e-5, 0.0));
        SmCN.setVal(1, 1, GIMLI::Complex( 1.0,  0.0));
        SmCN.setVal(1, 2, GIMLI::Complex( 3.0,  0.0));
        SmCN.setVal(2, 0, GIMLI::Complex( -2.0,  1.0));
        SmCN.setVal(2, 1, GIMLI::Complex( 3.0,  0.0));
        SmCN.setVal(2, 2, GIMLI::Complex(42.0,  0.0));
        testCHOLMODSolve< GIMLI::CSparseMapMatrix, GIMLI::Complex >(SmCN);
        
    }
    
};

CPPUNIT_TEST_SUITE_REGISTRATION(GIMLIExternalTest);

#endif