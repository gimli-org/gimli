#include <cppunit/extensions/HelperMacros.h>

#include <gimli.h>
#include <pos.h>
#include <line.h>
#include <plane.h>

class GeometryTest : public CppUnit::TestFixture{
    CPPUNIT_TEST_SUITE(GeometryTest);
    CPPUNIT_TEST(testEquality);
    CPPUNIT_TEST(testTouch);
    CPPUNIT_TEST(testIntersection);

    //CPPUNIT_TEST_EXCEPTION(funct, exception);
    CPPUNIT_TEST_SUITE_END();

public:
    void testEquality(){

    }

    void testTouch(){

    }

    void testIntersection(){
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


        GIMLI::RVector3 t = l.intersect(GIMLI::RVector3(0.0, 0.0),
                                        GIMLI::RVector3(1.0, 1.0));
        CPPUNIT_ASSERT(t.valid() == false);
        t = l.intersect(GIMLI::RVector3(10.0, 10.0),
                        GIMLI::RVector3(-1.0, -1.0));
        CPPUNIT_ASSERT(t.valid() == false);

        t = l.intersect(GIMLI::RVector3(10.0, 10.0),
                        GIMLI::RVector3(-1.0, -1.0));
        CPPUNIT_ASSERT(t.valid() == false);

        t = l.intersect(GIMLI::RVector3(10.0, 11.0),
                        GIMLI::RVector3(-1.0, -1.0));
        CPPUNIT_ASSERT(t.valid() == false);


    }

private:

};

CPPUNIT_TEST_SUITE_REGISTRATION(GeometryTest);
