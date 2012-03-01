#ifndef _GIMLI_TESTGIMLIMISC__H
#define _GIMLI_TESTGIMLIMISC__H

#include <cppunit/extensions/HelperMacros.h>

#include <gimli.h>
#include <ipcClient.h>
#include <memwatch.h>
#include <matrix.h>

class GIMLIMiscTest : public CppUnit::TestFixture  {
    CPPUNIT_TEST_SUITE( GIMLIMiscTest );
    CPPUNIT_TEST( testGimliMisc );
    CPPUNIT_TEST( testBooleanLogic );
    CPPUNIT_TEST( testFunctorTemplates );
    CPPUNIT_TEST( testStringFunctions );
    CPPUNIT_TEST( testIPCSHM );
    CPPUNIT_TEST( testMemWatch );
//     CPPUNIT_TEST( testRotationByQuaternion );
    
	//CPPUNIT_TEST_EXCEPTION( funct, exception );
    CPPUNIT_TEST_SUITE_END();
    
public:    
	void testGimliMisc(){
		std::cout << "Hello, this is: " << GIMLI::versionStr() << std::endl;
		CPPUNIT_ASSERT( GIMLI::fileExist( "unittest.sh" ) == true );
        std::cout << "number of CPU: " << GIMLI::numberOfCPU() << std::endl;
        std::cout << "sizes: int" << " " << sizeof( int ) 
                  << " long" << " " << sizeof( long )
                  << " size_t" << " " << sizeof( size_t )
                  << " ptr"<< " " << sizeof( void *)<< std::endl;
    }

    void testBooleanLogic(){
       
        std::set< int > s1; s1.insert( 1 ); s1.insert( 2 ); s1.insert( 3 ); s1.insert( 4 );
        std::set< int > s2; s2.insert( 0 ); s2.insert( 2 ); s2.insert( 3 ); s2.insert( 5 );  
        std::set< int > s3; s3.insert( 1 ); s3.insert( 2 ); s3.insert( 6 ); s3.insert( 5 );
        std::set< int > s4; s4.insert( 0 ); s4.insert( 2 ); s4.insert( 3 ); s4.insert( 4 );    
        
        std::set< int > stest;
        GIMLI::intersectionSet( stest, s1, s2 );
        CPPUNIT_ASSERT( stest.size() == 2 );
        CPPUNIT_ASSERT( *stest.begin() == 2 );
        GIMLI::intersectionSet( stest, s1, s2, s3 );
        CPPUNIT_ASSERT( stest.size() == 1 );
        CPPUNIT_ASSERT( *stest.begin() == 2 );
        
        std::vector < std::set < int > > setVec;
        setVec.push_back( s1 );
        setVec.push_back( s2 );
        setVec.push_back( s3 );
        setVec.push_back( s4 );
        GIMLI::intersectionSet( stest, setVec );
        CPPUNIT_ASSERT( stest.size() == 1 );
        CPPUNIT_ASSERT( *stest.begin() == 2 );
        
    }
    
    void testFunctorTemplates(){
//         int * t = new int[4];
//         std::cout << t << std::endl;
//         GIMLI::deletePtr()(t);
//         std::cout << t << std::endl;
//         CPPUNIT_ASSERT( t == NULL );
    }
    
    void testStringFunctions(){
        std::string t1( "a:bb:ccc:ddd" );
        GIMLI::split( t1, ':' );
        CPPUNIT_ASSERT( GIMLI::split( t1, ':' ).size() == 4 );
        CPPUNIT_ASSERT( GIMLI::split( t1, ':' )[ 0 ] == "a" );
        CPPUNIT_ASSERT( GIMLI::split( t1, ':' )[ 3 ] == "ddd" );
    }
    
    void testRotationByQuaternion(){
        
    }

    void testIPCSHM(){
        
        // Init shared memory
        GIMLI::IPCClientSHM ipc;

        // Init shared memory
        ipc.setSegmentName( "unittest" );
        ipc.setInt( "testInt", 1 );
        ipc.setBool( "testBool", false );
        ipc.setDouble( "testDouble", 3.14 );
        
        // this can by done by any other process or program on the same machine
        GIMLI::IPCClientSHM ipc2;
        ipc2.setSegmentName( "unittest" );
        ipc2.info();
        CPPUNIT_ASSERT( ipc2.getInt("testInt" ) == 1 );
        CPPUNIT_ASSERT( ipc2.getBool("testBool" ) == false );
        CPPUNIT_ASSERT( ipc2.getDouble("testDouble" ) == 3.14 );
        
        // free the shared memory
        ipc.free( "unittest" );
    }
    
    void testMemWatch(){
        GIMLI::__GIMLI_DEBUG__ = true;
        GIMLI::MemWatch::singleton().info( WHERE );
        double * mat2 = new double[ 10000 * 10000 ];
        GIMLI::MemWatch::singleton().info( WHERE );
        delete [] mat2;
        GIMLI::MemWatch::singleton().info( WHERE );
        GIMLI::RMatrix mat( 10000, 10000 );
        GIMLI::MemWatch::singleton().info( WHERE );
        mat.clear();
        GIMLI::MemWatch::singleton().info( WHERE );
    }

};


CPPUNIT_TEST_SUITE_REGISTRATION( GIMLIMiscTest );

#endif