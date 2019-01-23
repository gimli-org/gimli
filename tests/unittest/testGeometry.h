#include <cppunit/extensions/HelperMacros.h>

#include <gimli.h>
#include <pos.h>
#include <line.h>
#include <plane.h>

class GeometryTest : public CppUnit::TestFixture{
    CPPUNIT_TEST_SUITE(GeometryTest);
    CPPUNIT_TEST(testEquality);
    CPPUNIT_TEST(testTouch);
    CPPUNIT_TEST(testIntersection2d);
    CPPUNIT_TEST(testIntersection3d);

    //CPPUNIT_TEST_EXCEPTION(funct, exception);
    CPPUNIT_TEST_SUITE_END();

public:
    void testEquality(){
    }

    void testTouch(){
        GIMLI::Line l(GIMLI::RVector3(0.0, 0.0), GIMLI::RVector3(1.0, 0.0));

        int p = 0;
        CPPUNIT_ASSERT(l.touch1(GIMLI::RVector3(-1.0, 0.0), p) == true);
        CPPUNIT_ASSERT(p == 1);
        CPPUNIT_ASSERT(l.touch1(GIMLI::RVector3(0.0, 0.0), p) == true);
        CPPUNIT_ASSERT(p == 2);
        CPPUNIT_ASSERT(l.touch1(GIMLI::RVector3(0.5, 0.0), p) == true);
        CPPUNIT_ASSERT(p == 3);
        CPPUNIT_ASSERT(l.touch1(GIMLI::RVector3(1.0, 0.0), p) == true);
        CPPUNIT_ASSERT(p == 4);
        CPPUNIT_ASSERT(l.touch1(GIMLI::RVector3(2.0, 0.0), p) == true);
        CPPUNIT_ASSERT(p == 5);
    }

    void testIntersection2d(){
        GIMLI::Line l(GIMLI::RVector3(0.0, 0.0), GIMLI::RVector3(10.0, 10.0));

        CPPUNIT_ASSERT(l.intersect(GIMLI::RVector3(0.0, 5.0),
                                   GIMLI::RVector3(1.0, 0.0)) == GIMLI::RVector3(5.0, 5.0));

        CPPUNIT_ASSERT(l.intersect(GIMLI::RVector3(10.0, 5.0),
                                   GIMLI::RVector3(-1.0, 0.0)) == GIMLI::RVector3(5.0, 5.0));

        CPPUNIT_ASSERT(l.intersect(GIMLI::RVector3(5.0, 10.0),
                                   GIMLI::RVector3(0.0, -1.0)) == GIMLI::RVector3(5.0, 5.0));
        CPPUNIT_ASSERT(l.intersect(GIMLI::RVector3(5.0, 0.0),
                                   GIMLI::RVector3(0.0, 1.0)) == GIMLI::RVector3(5.0, 5.0));

        CPPUNIT_ASSERT(l.intersect(GIMLI::RVector3(0.0, 10.0),
                                   GIMLI::RVector3(1.0, -1.0)) == GIMLI::RVector3(5.0, 5.0));
        CPPUNIT_ASSERT(l.intersect(GIMLI::RVector3(10.0, 0.0),
                                   GIMLI::RVector3(-1.0, 1.0)) == GIMLI::RVector3(5.0, 5.0));

        CPPUNIT_ASSERT(l.intersect(GIMLI::RVector3(0.0, 10.0),
                                   GIMLI::RVector3(0.0, -1.0)) == GIMLI::RVector3(0.0, 0.0));

        CPPUNIT_ASSERT(l.intersect(GIMLI::RVector3(0.0, 0.0),
                                   GIMLI::RVector3(0.0, 1.0)) == GIMLI::RVector3(0.0, 0.0));

        CPPUNIT_ASSERT(l.intersect(GIMLI::RVector3(10.0, 0.0),
                                   GIMLI::RVector3(0.0,  1.0)) == GIMLI::RVector3(10.0, 10.0));
        CPPUNIT_ASSERT(l.intersect(GIMLI::RVector3(10.0, 10.0),
                                   GIMLI::RVector3(0.0,  1.0)) == GIMLI::RVector3(10.0, 10.0));


        // parallel and distance 0.0, return true but t is invalid 
        GIMLI::RVector3 t;
        CPPUNIT_ASSERT(l.intersectRay(GIMLI::RVector3(0.0, 0.0), GIMLI::RVector3(1.0, 1.0), t) == true);
        
        t = l.intersect(GIMLI::RVector3(0.0, 0.0), GIMLI::RVector3(1.0, 1.0));
        CPPUNIT_ASSERT(t.valid() == false);

        t = l.intersect(GIMLI::RVector3(10.0, 10.0), GIMLI::RVector3(-1.0, -1.0));
        CPPUNIT_ASSERT(t.valid() == false);

        // parallel and distance > 0.0, return false and t is invalid 
        l.intersectRay(GIMLI::RVector3(10.0, 11.0), GIMLI::RVector3(-1.0, -1.0), t);
        CPPUNIT_ASSERT(l.intersectRay(GIMLI::RVector3(10.0, 11.0), GIMLI::RVector3(-1.0, -1.0), t) == false);

        t = l.intersect(GIMLI::RVector3(10.0, 11.0), GIMLI::RVector3(-1.0, -1.0));
        CPPUNIT_ASSERT(t.valid() == false);

        CPPUNIT_ASSERT(l.intersectRay(GIMLI::RVector3(1.0, 0.0), GIMLI::RVector3(1.0, 1.0), t) == false);
        
        t = l.intersect(GIMLI::RVector3(1.0, 0.0), GIMLI::RVector3(1.0, 1.0));
        CPPUNIT_ASSERT(t.valid() == false);
    }

    void testIntersection3d(){
        GIMLI::RVector3 p0(0.0, 0.0, 0.0);
        GIMLI::RVector3 p1(1.0, 0.0, 0.0);
        GIMLI::RVector3 t;

        GIMLI::Line l(p0, p1);
        GIMLI::RVector3 tp((p1-p0)/2.0);
        
        GIMLI::RVector3 pq(0.0, 1.0, 0.0);
        CPPUNIT_ASSERT(l.intersect(pq, (tp-pq)) == tp);
        
        pq = GIMLI::RVector3(0.0, -1.0, 0.0);
        CPPUNIT_ASSERT(l.intersect(pq, (tp-pq)) == tp);

        pq = GIMLI::RVector3(0.0, -1.0, -1.0);
        CPPUNIT_ASSERT(l.intersect(pq, (tp-pq)) == tp);

        pq = GIMLI::RVector3(-14.0, -1.0, -1.0);
        CPPUNIT_ASSERT(l.intersect(pq, 3.14*(tp-pq)) == tp);

        // wrong ray direction
        pq = GIMLI::RVector3(-14.0, -1.0, -1.0);
        CPPUNIT_ASSERT(l.intersectRay(pq, -3.14 * (tp-pq), t) == false);
        CPPUNIT_ASSERT(t == tp);

    }

private:

};

CPPUNIT_TEST_SUITE_REGISTRATION(GeometryTest);
