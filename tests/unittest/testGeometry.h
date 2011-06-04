#include <cppunit/extensions/HelperMacros.h>

#include <gimli.h>
#include <pos.h>
#include <line.h>
#include <plane.h>

class GeometryTest : public CppUnit::TestFixture  {
    CPPUNIT_TEST_SUITE( GeometryTest );
    CPPUNIT_TEST( testEquality );
    CPPUNIT_TEST( testTouch );
    CPPUNIT_TEST( testIntersection );
    
    //CPPUNIT_TEST_EXCEPTION( funct, exception );
    CPPUNIT_TEST_SUITE_END();
    
public:    
    void testEquality(){
    
    }

    void testTouch(){
    
    }

    void testIntersection(){
    }
    
private:
       
};

CPPUNIT_TEST_SUITE_REGISTRATION( GeometryTest );