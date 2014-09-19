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
    
    void testCHOLMOD(){
        GIMLI::RSparseMapMatrix Sm(3,3);
        Sm.addVal(0,0,1.0);
        Sm.addVal(1,1,1.0);
        Sm.addVal(2,2,1.0);
        GIMLI::RSparseMatrix S(Sm);
        GIMLI::CHOLMODWrapper solver(S, true);

        GIMLI::RVector ones(Sm.rows(), 1.0);
        GIMLI::RVector x(Sm.rows());
        solver.solve(ones, x);
        std::cout << GIMLI::norm(ones - S * x) << std::endl;
//         3 3 5 -1
//         1 1  1.  0.
//         3 1  2. -1.
//         2 2  1.  0.
//         3 2  3.  0.
//         3 3 42.  0.

        return;
        GIMLI::CSparseMapMatrix Smc(3,3);
        Smc.addVal(0, 0, Complex( 1.0,  0.0));
        Smc.addVal(2, 0, Complex( 2.0, -1.0));
        Smc.addVal(1, 1, Complex( 1.0,  0.0));
        Smc.addVal(2, 1, Complex( 3.0,  0.0));
        Smc.addVal(2, 2, Complex(42.0,  0.0));
        
    
        GIMLI::CSparseMatrix Sc(Smc);
        GIMLI::CHOLMODWrapper solverC(Sc, true);
        GIMLI::CVector bC(Sm.rows(), Complex(1.0, 0.0));
        GIMLI::CVector xC(Sm.rows());
        solver.solve(bC, xC);
        std::cout << GIMLI::norm(xC) << std::endl;
        std::cout << GIMLI::norm(Sc * xC) << std::endl;
        std::cout << GIMLI::norm(bC - Sc * xC) << std::endl;
        
//std::cout << norm(ones - S * solverC.solve(toComplex(ones))) << std::endl;
        
                

    }
    
};

CPPUNIT_TEST_SUITE_REGISTRATION(GIMLIExternalTest);

#endif