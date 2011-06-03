#ifndef _GIMLI_TESTGIMLIMISC__H
#define _GIMLI_TESTGIMLIMISC__H

#include <cppunit/extensions/HelperMacros.h>

#include <gimli.h>

class GIMLIMiscTest : public CppUnit::TestFixture  {
    CPPUNIT_TEST_SUITE( GIMLIMiscTest );
    CPPUNIT_TEST( testGimliMisc );
	CPPUNIT_TEST( testBooleanLogic );
    CPPUNIT_TEST( testFunctorTemplates );
    CPPUNIT_TEST( testStringFunctions );
    CPPUNIT_TEST( testRotationByQuaternion );
    
	//CPPUNIT_TEST_EXCEPTION( funct, exception );
    CPPUNIT_TEST_SUITE_END();
    
public:    
	void testGimliMisc(){
		std::cout << "Hello, this is: " << GIMLI::versionStr() << std::endl;
		CPPUNIT_ASSERT( GIMLI::fileExist( "unittest.sh" ) == true );
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
        std::string t1( "a:bb:ccc:ddd");
        GIMLI::split( t1, ':' );
        //CPPUNIT_ASSERT( GIMLI::split( t1, ':' ).size() == 4 );
	
        //CPPUNIT_ASSERT( GIMLI::split( t1, ':' )[ 0 ] == "a" );
        //CPPUNIT_ASSERT( GIMLI::split( t1, ':' )[ 3 ] == "ddd" );
    }
    
    void testRotationByQuaternion(){
        
    }

};


CPPUNIT_TEST_SUITE_REGISTRATION( GIMLIMiscTest );

#endif