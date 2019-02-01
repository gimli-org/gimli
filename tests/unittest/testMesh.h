#include <cppunit/extensions/HelperMacros.h>

#include <gimli.h>
#include <mesh.h>
#include <meshgenerators.h>

#include <stdexcept>

using namespace GIMLI;

class MeshTest : public CppUnit::TestFixture{
    CPPUNIT_TEST_SUITE(MeshTest);
    CPPUNIT_TEST(testSimple);
    CPPUNIT_TEST(testRefine2d);
    CPPUNIT_TEST(testRefine3d);
        
    CPPUNIT_TEST(testPolygonInsertion);
    
    //CPPUNIT_TEST_EXCEPTION(funct, exception);
    CPPUNIT_TEST_SUITE_END();
    
public:

    void testSimple(){
        Node n0(RVector3(0.0, 0.0, 0.0));
        Node n1(RVector3(1.0, 0.0, 0.0));
        Node n2(RVector3(0.0, 1.0, 0.0));

        Triangle *tri = new Triangle (n0, n1, n2);
        delete tri;
    }
    
    void testRefine2d(){
                
        Mesh mesh(2);
        
        Node *n0 = mesh.createNode(RVector3(0.0, 0.0));
        Node *n1 = mesh.createNode(RVector3(1.0, 0.0));
        Node *n2 = mesh.createNode(RVector3(0.0, 1.0));
        Node *n3 = mesh.createNode(RVector3(0.0, -1.0));
        Node *n4 = mesh.createNode(RVector3(1.0, -1.0));
        
        mesh.createTriangle(*n0, *n1, *n2);
        mesh.createQuadrangle(*n0, *n3, *n4, *n1);
        
        mesh.createNeighbourInfos();
        
        // split tri into 4 tris, quad into 4 quads
        CPPUNIT_ASSERT(mesh.cellCount() == 2);
        CPPUNIT_ASSERT(mesh.nodeCount() == 5);
        CPPUNIT_ASSERT(mesh.boundaryCount() == 6);
                
//         mesh.save("mesh.bms");
//         mesh.showInfos();
        
        Mesh tmp(mesh.createH2());
        tmp.createNeighbourInfos();
        
//         tmp.save("tmph2.bms");
//         tmp.showInfos();
        
//         
        // split tri into 4 tris, quad into 4 quads
        CPPUNIT_ASSERT(tmp.cellCount() == 8);
        CPPUNIT_ASSERT(tmp.nodeCount() == 12);
        CPPUNIT_ASSERT(tmp.boundaryCount() == 19);
                 
        tmp.findCell(RVector3(0.5, 0.5));
        
        
        tmp = mesh.createP2();
        
//         tmp.save("tmpp2.bms");
//         tmp.showInfos();
        // split tri into 1 tri6, quad into 1 quad8
        CPPUNIT_ASSERT(tmp.cellCount() == 2);
        CPPUNIT_ASSERT(tmp.nodeCount() == 11);
        CPPUNIT_ASSERT(tmp.boundaryCount() == 6);
    }
    
     void testRefine3d(){
                
        Mesh mesh(3);
        
        Node *n0 = mesh.createNode(RVector3(0.0, 0.0));
        Node *n1 = mesh.createNode(RVector3(1.0, 0.0));
        Node *n2 = mesh.createNode(RVector3(0.0, 1.0));
        Node *n3 = mesh.createNode(RVector3(0.0, 0.0, 1.0));
        
        mesh.createTetrahedron(*n0, *n1, *n2, *n3);
        mesh.createNeighbourInfos();
        
        CPPUNIT_ASSERT(mesh.cellCount() == 1);
        CPPUNIT_ASSERT(mesh.nodeCount() == 4);
        // boundaries with marker 0 will not yet refined
        CPPUNIT_ASSERT(mesh.boundaryCount() == 4);
                
        Mesh tmp(mesh.createH2());
       
        CPPUNIT_ASSERT(tmp.cellCount() == 8);
        CPPUNIT_ASSERT(tmp.nodeCount() == 10);
        CPPUNIT_ASSERT(tmp.boundaryCount() == 16);
        
        tmp.createNeighbourInfos();
        CPPUNIT_ASSERT(tmp.boundaryCount() == 24);
        
        tmp = mesh.createP2();
        CPPUNIT_ASSERT(tmp.cellCount() == 1);
        CPPUNIT_ASSERT(tmp.nodeCount() == 10);
        CPPUNIT_ASSERT(tmp.boundaryCount() == 4);
        
        Mesh q(3);
        std::vector < Node * > nodes;
        nodes.push_back(q.createNode(RVector3(0.0, 0.0, 0.0)));
        nodes.push_back(q.createNode(RVector3(1.0, 0.0, 0.0)));
        nodes.push_back(q.createNode(RVector3(1.0, 1.0, 0.0)));
        nodes.push_back(q.createNode(RVector3(0.0, 1.0, 0.0)));
        
        nodes.push_back(q.createNode(RVector3(0.0, 0.0, 1.0)));
        nodes.push_back(q.createNode(RVector3(1.0, 0.0, 1.0)));
        nodes.push_back(q.createNode(RVector3(1.0, 1.0, 1.0)));
        nodes.push_back(q.createNode(RVector3(0.0, 1.0, 1.0)));
        
        q.createCell(nodes);
        q = q.createH2();
                
        CPPUNIT_ASSERT(q.cellCount() == 8);
        CPPUNIT_ASSERT(q.nodeCount() == 27);
        
        q = createMesh3D(1u, 1u, 1u, 0);
        q = q.createH2();

        CPPUNIT_ASSERT(q.cellCount() == 8);
        CPPUNIT_ASSERT(q.nodeCount() == 27);
    }
    
    void testPolygonInsertion(){
        GIMLI::Node n0(0.0, 0.0, -0.5); n0.setId(0);
        GIMLI::Node n1(1.0, 0.0, -0.5); n1.setId(1);
        GIMLI::Node n2(1.0, 0.0,  0.5); n2.setId(2);
        GIMLI::Node n3(0.0, 0.0,  0.5); n3.setId(3);
        GIMLI::Node n4(0.0, 0.5, -0.5); n4.setId(4);
        GIMLI::Node n5(1.0, 0.5, -0.5); n5.setId(5);
        GIMLI::Node n6(1.0, 0.5,  0.5); n6.setId(6);
        GIMLI::Node n7(0.0, 0.5,  0.5); n7.setId(7);
        GIMLI::Node n8(0.5, 0.0, -0.5); n8.setId(8);

        std::vector< Node * > nodes;
        nodes.push_back(&n0);
        nodes.push_back(&n4);
        nodes.push_back(&n5);
        nodes.push_back(&n1);

        GIMLI::PolygonFace q1(nodes);


        q1.insertNode(&n8);
        CPPUNIT_ASSERT(q1.secondaryNodes().size() == 0);
        CPPUNIT_ASSERT(q1.allNodeCount() == 5);
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION(MeshTest);