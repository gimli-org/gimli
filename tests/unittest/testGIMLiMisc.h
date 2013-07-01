#ifndef _GIMLI_TESTGIMLIMISC__H
#define _GIMLI_TESTGIMLIMISC__H

#include <cppunit/extensions/HelperMacros.h>

#include <gimli.h>
#include <ipcClient.h>
#include <memwatch.h>

#include <matrix.h>

#include <polynomial.h>
#include <pos.h>

class GIMLIMiscTest : public CppUnit::TestFixture  {
    CPPUNIT_TEST_SUITE( GIMLIMiscTest );
    CPPUNIT_TEST( testGimliMisc );
    CPPUNIT_TEST( testBooleanLogic );
    CPPUNIT_TEST( testFunctorTemplates );
    CPPUNIT_TEST( testStringFunctions );
    CPPUNIT_TEST( testIPCSHM );
    CPPUNIT_TEST( testMemWatch );
    CPPUNIT_TEST( testPolynomialFunction );
//     CPPUNIT_TEST( testRotationByQuaternion );
    
	//CPPUNIT_TEST_EXCEPTION( funct, exception );
    CPPUNIT_TEST_SUITE_END();
    
public:    
	void testGimliMisc(){
		std::cout << "Hello, this is: " << GIMLI::versionStr() << std::endl;
        std::cout << "I'm calling from" << WHERE_AM_I << std::endl;
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
        std::cout << "MemWatch" << std::endl;
        GIMLI::__GIMLI_DEBUG__ = true;
        GIMLI::MemWatch::instance().info(WHERE);
        double * mat2 = new double[10000 * 10000];
        GIMLI::MemWatch::instance().info(WHERE);
        delete [] mat2;
        GIMLI::MemWatch::instance().info(WHERE);
        GIMLI::RMatrix mat(10000, 10000);
        GIMLI::MemWatch::pInstance()->info(WHERE);
        mat.clear();
        GIMLI::MemWatch::pInstance()->info(WHERE);
    }
    
    void testPolynomialFunction(){
        CPPUNIT_ASSERT( GIMLI::PolynomialFunction< double > ( GIMLI::RVector( 0.0 ) )( GIMLI::RVector3( 3.14, 0.0, 0.0 ) ) == 0.0 );
        
        GIMLI::RVector ax( 3, 0 ); ax[ 1 ] = 1.0;
        GIMLI::RVector ay( 3, 0 ); ay[ 1 ] = 2.0;
        GIMLI::RVector az( 3, 0 ); az[ 1 ] = 3.0;
        
        GIMLI::PolynomialFunction< double > f( ax );
        CPPUNIT_ASSERT( f( GIMLI::RVector3( 3.14, 0.0, 0.0 ) ) == ax[ 1 ] * 3.14 );
        CPPUNIT_ASSERT( f.derive( 0 )( GIMLI::RVector3( 3.14, 0.0, 0.0 ) ) == ax[ 1 ] );
        CPPUNIT_ASSERT( f.derive( 1 )( GIMLI::RVector3( 3.14, 0.0, 0.0 ) ) == 0.0 );
        CPPUNIT_ASSERT( f.derive( 2 )( GIMLI::RVector3( 3.14, 0.0, 0.0 ) ) == 0.0 );
        
        GIMLI::PolynomialFunction< double > g( GIMLI::RVector(0), ay );
        CPPUNIT_ASSERT( g( GIMLI::RVector3( 3.14, 0.0, 0.0 ) ) == 0.0 );
        CPPUNIT_ASSERT( g( GIMLI::RVector3( 0.0, 3.14, 0.0 ) ) == ay[ 1 ] * 3.14 );
        CPPUNIT_ASSERT( g( GIMLI::RVector3( 0.0, 0.0, 3.14 ) ) == 0.0 );
        CPPUNIT_ASSERT( g.derive( 0 )( GIMLI::RVector3( 3.14, 0.0, 0.0 ) ) == 0.0 );
        CPPUNIT_ASSERT( g.derive( 1 )( GIMLI::RVector3( 3.14, 0.0, 0.0 ) ) == ay[ 1 ] );
        CPPUNIT_ASSERT( g.derive( 2 )( GIMLI::RVector3( 3.14, 0.0, 0.0 ) ) == 0.0 );
        
        GIMLI::PolynomialFunction< double > h( GIMLI::RVector(0), GIMLI::RVector(0), az );
        CPPUNIT_ASSERT( h( GIMLI::RVector3( 3.14, 0.0, 0.0 ) ) == 0.0 );
        CPPUNIT_ASSERT( h( GIMLI::RVector3( 0.0, 3.14, 0.0 ) ) == 0.0 );
        CPPUNIT_ASSERT( h( GIMLI::RVector3( 0.0, 0.0, 3.14 ) ) == az[ 1 ] * 3.14 );
        CPPUNIT_ASSERT( h.derive( 0 )( GIMLI::RVector3( 3.14, 0.0, 0.0 ) ) == 0.0 );
        CPPUNIT_ASSERT( h.derive( 1 )( GIMLI::RVector3( 3.14, 0.0, 0.0 ) ) == 0.0 );
        CPPUNIT_ASSERT( h.derive( 2 )( GIMLI::RVector3( 3.14, 0.0, 0.0 ) ) == az[ 1 ] );
        
//        std::cout << f << std::endl;
//         std::cout << g << std::endl;
//         std::cout << h << std::endl;
        ax[ 1 ] = -1;
        ay[ 2 ] = 0;
        f = GIMLI::PolynomialFunction< double >( ax );
        g = GIMLI::PolynomialFunction< double >( GIMLI::RVector(0), ay );
        

        GIMLI::PolynomialFunction< double > t( 2 );
        GIMLI::RVector tmp( GIMLI::powInt( 2, 3 ), 0.0 );
        tmp[ 3 ] = 1.0; t.fill( tmp );//t= +1xy
        CPPUNIT_ASSERT( t( GIMLI::RVector3( 3.0, 3.0, 0.0 ) ) == 3.0*3.0 );
        CPPUNIT_ASSERT( t.derive(0)( GIMLI::RVector3( 3.0, 3.0, 0.0 ) ) == 3.0 );
        CPPUNIT_ASSERT( t.derive(1)( GIMLI::RVector3( 3.0, 3.0, 0.0 ) ) == 3.0 );
        CPPUNIT_ASSERT( t.derive(2)( GIMLI::RVector3( 3.0, 3.0, 0.0 ) ) == 0.0 );
        
//         std::cout << t << std::endl;
//         std::cout << t.derive(0) << std::endl;
//         std::cout << t.derive(1) << std::endl;
//         std::cout << t.derive(2) << std::endl;
// 
//         exit(0);
        //         
//         std::cout << f * g << std::endl;
//         std::cout << (f * g).derive( 0 ) << std::endl;
        
    }
    
};

CPPUNIT_TEST_SUITE_REGISTRATION( GIMLIMiscTest );

#endif