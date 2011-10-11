#include <cppunit/extensions/HelperMacros.h>

#include <gimli.h>
#include <datacontainer.h>
#include <pos.h>

#include <stdexcept>

using namespace GIMLI;

class DataContainerTest : public CppUnit::TestFixture  {
    CPPUNIT_TEST_SUITE( DataContainerTest );
    CPPUNIT_TEST( testCreate );
    
    CPPUNIT_TEST( testIO );
    CPPUNIT_TEST( testEdit );
    
    CPPUNIT_TEST_SUITE_END();
    
public:    
    
    void testCreate(){
        DataContainer data;
        uint nSensors = 10;
        
        for ( uint i = 0; i < nSensors; i++ ){
            data.createSensor( RVector3( double(i), 0.0, 0.0 ) );
        }
        CPPUNIT_ASSERT( data.sensorCount() == nSensors);
        
        long id = data.createSensor( data.sensorPositions()[ ( nSensors - 1 ) ] + 0.01, 0.1 );
        
        CPPUNIT_ASSERT( id == ((long)nSensors - 1) );
        CPPUNIT_ASSERT( data.sensorCount() == nSensors );
        
        data.registerSensorIndex( "S1" );
        data.registerSensorIndex( "S2" );
        
        data.resize( nSensors );
        
        RVector tmp( data.size() ); tmp.fill( x__ );
        
        // create sensor index data [0,9]
        data.set("S1", tmp );
        // create sensor index data [1,10]
        data.set("S2", tmp+1.0);
        
        // create values [1,10]
        data.set("vals", tmp );
        
        // set all valid
        data.set("valid", RVector( nSensors, 1.0 ) );
                
        data.markInvalidSensorIndices();
        data.removeInvalid();
        CPPUNIT_ASSERT( data.size() == nSensors - 1 );
        
        data.remove( find( ( data("S1") == 3 ) | ( data("S2") == 3 ) ) );
        CPPUNIT_ASSERT( data.size() == nSensors - 3 );
        
        data.removeUnusedSensors();
        CPPUNIT_ASSERT( data.sensorCount() == nSensors - 1 );
        
        data.removeSensorIdx( 2 );
                
        CPPUNIT_ASSERT( data.sensorCount() == nSensors - 2 );
        CPPUNIT_ASSERT( data.size()        == nSensors - 4 );
        
        data.save( "test0.dat" );        
        
        DataContainer data2( data );
        data.add( data2 );
        CPPUNIT_ASSERT( data.sensorCount() == data2.sensorCount() );
        CPPUNIT_ASSERT( data.size()        == data2.size()*2 );
        
        DataContainer data3( data2 );
        data3.setSensorPosition(0, RVector3( 0.5, .0, .0 ) );
        data.add( data3 );
        CPPUNIT_ASSERT( data.sensorCount() == data2.sensorCount()+1 );
        CPPUNIT_ASSERT( data.size()        == data2.size()*3 );
        
        data.save( "test1.dat" );        
        
    }   
    
    void testIO(){
        
    }
    
    void testEdit(){
        
    }
    
    

private:
   
    
};

CPPUNIT_TEST_SUITE_REGISTRATION( DataContainerTest );