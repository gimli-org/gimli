#include <cppunit/extensions/HelperMacros.h>

#include <gimli.h>
#include <datacontainer.h>
#include <pos.h>

#include <stdexcept>

using namespace GIMLI;
Placeholder x__;

class DataContainerTest : public CppUnit::TestFixture  {
    CPPUNIT_TEST_SUITE( DataContainerTest );
    CPPUNIT_TEST( testCreate );
    CPPUNIT_TEST( testIO );
    CPPUNIT_TEST( testEdit );
    
    CPPUNIT_TEST_SUITE_END();
    
public:    
    
    void testCreate(){
        DataContainer data;
        uint nProbes = 10;
        
        for ( uint i = 0; i < nProbes; i++ ){
            data.createSensor( RVector3( double(i), 0.0, 0.0 ) );
        }
        CPPUNIT_ASSERT( data.sensorCount() == nProbes );
        
        long id = data.createSensor( data.sensorPositions()[( nProbes - 1 )] + 0.01, 0.1 );
        
        CPPUNIT_ASSERT( id == ((long)nProbes - 1) );
        CPPUNIT_ASSERT( data.sensorCount() == nProbes );
        
        data.registerSensorIndex( "S1" );
        data.registerSensorIndex( "S2" );
        
        data.resize( nProbes );
        
        RVector tmp( data.size() ); tmp.fill( x__ );
        
        data.set("S1", tmp );
        data.set("S2", tmp+1.0);
        
        data.set("vals", tmp );
        data.set("valid", RVector( nProbes, 1.0 ) );
        data.save( "test.dat" );
        data.remove( find( ( data("S1") == 1 ) & ( data("S2") == 1 ) ) );
        
        
        CPPUNIT_ASSERT( data.size() == nProbes-2 );
        
    }
    
    void testIO(){
        
    }
    
    void testEdit(){
        
    }
    
    

private:
   
    
};

CPPUNIT_TEST_SUITE_REGISTRATION( DataContainerTest );