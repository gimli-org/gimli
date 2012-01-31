#include <cppunit/extensions/HelperMacros.h>

#include <gimli.h>
#include <mesh.h>

#include <stdexcept>

using namespace GIMLI;

class MeshTest : public CppUnit::TestFixture  {
    CPPUNIT_TEST_SUITE( MeshTest );
    CPPUNIT_TEST( testRefine2d );
    CPPUNIT_TEST( testRefine3d );
        
    //CPPUNIT_TEST_EXCEPTION( funct, exception );
    CPPUNIT_TEST_SUITE_END();
    
public:    
    
    void testRefine2d(){
                
        Mesh mesh( 2 );
        
        Node *n0 = mesh.createNode( RVector3( 0.0, 0.0 ) );
        Node *n1 = mesh.createNode( RVector3( 1.0, 0.0 ) );
        Node *n2 = mesh.createNode( RVector3( 0.0, 1.0 ) );
        Node *n3 = mesh.createNode( RVector3( 0.0, -1.0 ) );
        Node *n4 = mesh.createNode( RVector3( 1.0, -1.0 ) );
        
        mesh.createTriangle( *n0, *n1, *n2 );
        mesh.createQuadrangle( *n0, *n3, *n4, *n1 );
        
        mesh.createNeighbourInfos();
        
        mesh.save("mesh.bms");
        mesh.showInfos();
        
        Mesh tmp;
        tmp.createH2Mesh( mesh );
//         tmp.save("tmph2.bms");
//         tmp.showInfos();
        
        // split tri into 4 tris, quad into 8 tris
        CPPUNIT_ASSERT( tmp.cellCount() == 12 );
        CPPUNIT_ASSERT( tmp.nodeCount() == 12 );
        CPPUNIT_ASSERT( tmp.boundaryCount() == 23 );
                 
        tmp.createP2Mesh( mesh );
        
//         tmp.save("tmpp2.bms");
//         tmp.showInfos();
        // split tri into 1 tri6, quad into 1 quad8
        CPPUNIT_ASSERT( tmp.cellCount() == 2 );
        CPPUNIT_ASSERT( tmp.nodeCount() == 11 );
        CPPUNIT_ASSERT( tmp.boundaryCount() == 6 );
    }
    
     void testRefine3d(){
                
        Mesh mesh( 3 );
        
        Node *n0 = mesh.createNode( RVector3( 0.0, 0.0 ) );
        Node *n1 = mesh.createNode( RVector3( 1.0, 0.0 ) );
        Node *n2 = mesh.createNode( RVector3( 0.0, 1.0 ) );
        Node *n3 = mesh.createNode( RVector3( 0.0, 0.0, 1.0) );
        
        mesh.createTetrahedron( *n0, *n1, *n2, *n3 );
        mesh.createNeighbourInfos();
        
        Mesh tmp;
        tmp.createH2Mesh( mesh );
       
        CPPUNIT_ASSERT( tmp.cellCount() == 8 );
        CPPUNIT_ASSERT( tmp.nodeCount() == 10 );
        // boundaries with marker 0 will not yet refined
        CPPUNIT_ASSERT( tmp.boundaryCount() == 0 );
                 
        tmp.createP2Mesh( mesh );
        CPPUNIT_ASSERT( tmp.cellCount() == 1 );
        CPPUNIT_ASSERT( tmp.nodeCount() == 10 );
        CPPUNIT_ASSERT( tmp.boundaryCount() == 4 );
    }
    
};

CPPUNIT_TEST_SUITE_REGISTRATION( MeshTest );