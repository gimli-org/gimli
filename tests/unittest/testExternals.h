#ifndef _GIMLI_TESTEXTERNALS__H
#define _GIMLI_TESTEXTERNALS__H

#include <cppunit/extensions/HelperMacros.h>

#include <gimli.h>
#include <triangleWrapper.h>
#include <mesh.h>
#include <pos.h>

class GIMLIExternalTest : public CppUnit::TestFixture  {
    CPPUNIT_TEST_SUITE( GIMLIExternalTest );
    CPPUNIT_TEST( testTriangle );
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
        GIMLI::Mesh outMesh(2);
        tri.generate( outMesh );
        CPPUNIT_ASSERT( outMesh.cellCount() == 1 );
        CPPUNIT_ASSERT( outMesh.nodeCount() == 3 );
        CPPUNIT_ASSERT( outMesh.boundaryCount() == 3 );
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( GIMLIExternalTest );

#endif