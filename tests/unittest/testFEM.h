#include <cppunit/extensions/HelperMacros.h>

#include <gimli.h>
#include <node.h>
#include <shape.h>
#include <pos.h>
#include <meshentities.h>
#include <elementmatrix.h>
#include <integration.h>

class FEMTest : public CppUnit::TestFixture  {
    CPPUNIT_TEST_SUITE( FEMTest );

    CPPUNIT_TEST( testFEMBasics );

    CPPUNIT_TEST( testFEM1D );
    CPPUNIT_TEST( testFEM2D );
    CPPUNIT_TEST( testFEM3D );

    CPPUNIT_TEST_SUITE_END();

public:

    void testFEMBasics(){
        for ( uint i = 1; i < 10; i ++ ){
//             std::cout << "n = " << i << " " << sum( IntegrationRules::instance().gauWeights( i ) )
//                       << " " << sum( IntegrationRules::instance().edgWeights( i ) )
//                       << " " << sum( IntegrationRules::instance().triGLWeights( i ) )
//                       << " " << sum( IntegrationRules::instance().quaWeights( i ) ) << std::endl;

            CPPUNIT_ASSERT( ::fabs( sum( IntegrationRules::instance().gauWeights( i ) ) - 2.0 ) <  TOLERANCE );

//             for ( uint j = 0; j < rules.triGLAbscissa( i ).size(); j ++ ){
//                 std::cout << rules.triGLAbscissa( i )[ j ] << " " <<  rules.triGLWeights( i )[ j ] << std::endl;
//             }

            
        }
    }
    
    void testFEM1D(){
        resetNodes( );
        testStiffness1D();
        for ( size_t i = 0; i < nodes_.size(); i ++ ) {
            nodes_[ i ]->scale( RVector3( 42.42, 0.0 ) );
            nodes_[ i ]->translate( RVector3( -42.42, 0.0 ) );
        }
        testStiffness1D();
    }
    void testFEM2D(){
        resetNodes( );
        testStiffness2D();
        for ( size_t i = 0; i < nodes_.size(); i ++ ) {
            nodes_[ i ]->scale( RVector3( 42.42, 14.0 ) );
            nodes_[ i ]->rotate( RVector3( 0.0, 0.0, -42.42 ) );
            nodes_[ i ]->translate( RVector3( -42.42, 14.0 ) );
        }
        testStiffness2D();
    }
    void testFEM3D(){
        resetNodes( );
        testStiffness3D();
        for ( size_t i = 0; i < nodes_.size(); i ++ ) {
            nodes_[ i ]->scale( RVector3( 42.42, 14.0, 0.4) );
            nodes_[ i ]->rotate( RVector3( 33.0, 13.0, -42.42 ) );
            nodes_[ i ]->translate( RVector3( -42.42, 14.0, 24.23 ) );
        }
        testStiffness3D();
    }

    void testStiffness1D(){
        std::vector < Node * > n( 2 );
        n[ 0 ] = nodes_[ 0 ]; n[ 1 ] = nodes_[ 1 ];
        
        RVector e2( 2 ); e2[ 0 ] = 0; e2[ 1 ] =  1; // x
        GIMLI::EdgeCell edg( n );
        
        GIMLI::RPolynomialFunction E2_2( e2 );
        GIMLI::RPolynomialFunction E2_1 = -( -1.0 + E2_2 );
        
        CPPUNIT_ASSERT( edg.createShapeFunctions()[ 0 ] == E2_1 );
        CPPUNIT_ASSERT( edg.createShapeFunctions()[ 1 ] == E2_2 );
        
        testStiffness( edg );

        n.resize( 3 );
        n[ 0 ] = nodes_[ 0 ]; n[ 1 ] = nodes_[ 1 ]; n[ 2 ] = nodes_[ 8 ];
        GIMLI::Edge3Cell edg3( n );
        
        CPPUNIT_ASSERT( edg3.createShapeFunctions()[ 0 ] == E2_1 * ( 2.0 * E2_1 + -1.0 ) );
        CPPUNIT_ASSERT( edg3.createShapeFunctions()[ 1 ] == E2_2 * ( 2.0 * E2_2 + -1.0 ) );
        CPPUNIT_ASSERT( edg3.createShapeFunctions()[ 2 ] == E2_1 * E2_2 * 4.0 );
        
        testStiffness( edg3 );
    }

    void testStiffness2D(){
        std::vector < Node * > n( 2 );

        n.resize( 3 );
        n[ 0 ] = nodes_[ 0 ]; n[ 1 ] = nodes_[ 1 ]; n[ 2 ] = nodes_[ 2 ];
        GIMLI::Triangle tri( n );

        RVector e2( 2 ); e2[ 0 ] = 0; e2[ 1 ] =  1; // x
        
        GIMLI::RPolynomialFunction T3_2( e2, RVector( 0 ) );
        GIMLI::RPolynomialFunction T3_3( RVector( 0 ), e2 );
        GIMLI::RPolynomialFunction T3_1 = -( -1.0 + T3_2 + T3_3 );
        
        CPPUNIT_ASSERT( tri.createShapeFunctions()[ 0 ] == T3_1 );
        CPPUNIT_ASSERT( tri.createShapeFunctions()[ 1 ] == T3_2 );
        CPPUNIT_ASSERT( tri.createShapeFunctions()[ 2 ] == T3_3 );
        testStiffness( tri );
        
        n.resize( 6 );
        n[ 0 ] = nodes_[ 0 ]; n[ 1 ] = nodes_[ 1 ];  n[ 2 ]  = nodes_[ 3 ];
        n[ 3 ] = nodes_[ 8 ]; n[ 4 ] = nodes_[ 20 ]; n[ 5 ] = nodes_[ 11 ];
        GIMLI::Triangle6 tri6( n );

        CPPUNIT_ASSERT( tri6.createShapeFunctions()[ 0 ] == T3_1 * ( 2.0 * T3_1 + -1.0 ) );
        CPPUNIT_ASSERT( tri6.createShapeFunctions()[ 1 ] == T3_2 * ( 2.0 * T3_2 + -1.0 ) );
        CPPUNIT_ASSERT( tri6.createShapeFunctions()[ 2 ] == T3_3 * ( 2.0 * T3_3 + -1.0 ) );
        CPPUNIT_ASSERT( tri6.createShapeFunctions()[ 3 ] == T3_1 * T3_2 * 4.0 );
        CPPUNIT_ASSERT( tri6.createShapeFunctions()[ 4 ] == T3_2 * T3_3 * 4.0 );
        CPPUNIT_ASSERT( tri6.createShapeFunctions()[ 5 ] == T3_3 * T3_1 * 4.0 );
        
        testStiffness( tri6 );
        
        n.resize( 4 );
        n[ 0 ] = nodes_[ 0 ]; n[ 1 ] = nodes_[ 1 ]; n[ 2 ] = nodes_[ 2 ];n[ 3 ] = nodes_[ 3 ];
        GIMLI::Quadrangle quad( n );
        
        GIMLI::RPolynomialFunction E2_2R( e2, RVector( 0 ) );
        GIMLI::RPolynomialFunction E2_1R = -( -1.0 + E2_2R );
        GIMLI::RPolynomialFunction E2_2S( RVector( 0 ), e2 );
        GIMLI::RPolynomialFunction E2_1S = -( -1.0 + E2_2S );
        
        GIMLI::RPolynomialFunction Q4_1 = E2_1R * E2_1S;
        GIMLI::RPolynomialFunction Q4_2 = E2_2R * E2_1S;
        GIMLI::RPolynomialFunction Q4_3 = E2_2R * E2_2S;
        GIMLI::RPolynomialFunction Q4_4 = E2_1R * E2_2S;
        
        CPPUNIT_ASSERT( quad.createShapeFunctions()[ 0 ] == Q4_1  );
        CPPUNIT_ASSERT( quad.createShapeFunctions()[ 1 ] == Q4_2  );
        CPPUNIT_ASSERT( quad.createShapeFunctions()[ 2 ] == Q4_3  );
        CPPUNIT_ASSERT( quad.createShapeFunctions()[ 3 ] == Q4_4  );
        
        testStiffness( quad );

        n.resize( 8 );
        n[ 0 ] = nodes_[ 0 ]; n[ 1 ] = nodes_[ 1 ]; n[ 2 ] = nodes_[ 2 ];  n[ 3 ] = nodes_[ 3 ];
        n[ 4 ] = nodes_[ 8 ]; n[ 5 ] = nodes_[ 9 ]; n[ 6 ] = nodes_[ 10 ]; n[ 7 ] = nodes_[ 11 ];
        GIMLI::Quadrangle8 quad8( n );
        
        GIMLI::RPolynomialFunction E3_3R = E2_1R * E2_2R * 4.0;
        GIMLI::RPolynomialFunction E3_3S = E2_1S * E2_2S * 4.0;
        
        GIMLI::RPolynomialFunction Q8_5 = E3_3R * E2_1S;
        GIMLI::RPolynomialFunction Q8_6 = E2_2R * E3_3S;
        GIMLI::RPolynomialFunction Q8_7 = E3_3R * E2_2S;
        GIMLI::RPolynomialFunction Q8_8 = E2_1R * E3_3S;

        GIMLI::RPolynomialFunction Q8_1 = Q4_1 - Q8_8*0.5 - Q8_5*0.5;
        GIMLI::RPolynomialFunction Q8_2 = Q4_2 - Q8_5*0.5 - Q8_6*0.5;
        GIMLI::RPolynomialFunction Q8_3 = Q4_3 - Q8_6*0.5 - Q8_7*0.5;
        GIMLI::RPolynomialFunction Q8_4 = Q4_4 - Q8_7*0.5 - Q8_8*0.5;
        
        CPPUNIT_ASSERT( quad8.createShapeFunctions()[ 0 ] == Q8_1 );
        CPPUNIT_ASSERT( quad8.createShapeFunctions()[ 1 ] == Q8_2 );
        CPPUNIT_ASSERT( quad8.createShapeFunctions()[ 2 ] == Q8_3 );
        CPPUNIT_ASSERT( quad8.createShapeFunctions()[ 3 ] == Q8_4 );
        CPPUNIT_ASSERT( quad8.createShapeFunctions()[ 4 ] == Q8_5 );
        CPPUNIT_ASSERT( quad8.createShapeFunctions()[ 5 ] == Q8_6 );
        CPPUNIT_ASSERT( quad8.createShapeFunctions()[ 6 ] == Q8_7 );
        CPPUNIT_ASSERT( quad8.createShapeFunctions()[ 7 ] == Q8_8 );
                
        testStiffness( quad8 );
    }

    void testStiffness3D(){
        std::vector < Node * > n( 4 );

        n.resize( 4 );
        n[ 0 ] = nodes_[ 0 ]; n[ 1 ] = nodes_[ 1 ]; n[ 2 ] = nodes_[ 3 ];n[ 3 ] = nodes_[ 4 ];
        GIMLI::Tetrahedron tet( n );

        RVector e2( 2 ); e2[ 0 ] = 0; e2[ 1 ] =  1; // x
        
        GIMLI::RPolynomialFunction T4_2( e2, RVector( 0 ), RVector( 0 ) );
        GIMLI::RPolynomialFunction T4_3( RVector( 0 ), e2, RVector( 0 ) );
        GIMLI::RPolynomialFunction T4_4( RVector( 0 ), RVector( 0 ), e2 );
        GIMLI::RPolynomialFunction T4_1 = - ( -1. + T4_2 + T4_3 + T4_4 );
        
        CPPUNIT_ASSERT( tet.createShapeFunctions()[ 0 ] == T4_1 );
        CPPUNIT_ASSERT( tet.createShapeFunctions()[ 1 ] == T4_2 );
        CPPUNIT_ASSERT( tet.createShapeFunctions()[ 2 ] == T4_3 );
        CPPUNIT_ASSERT( tet.createShapeFunctions()[ 3 ] == T4_4 );
       
        testStiffness( tet );

        /*    *------*  \n
             /|     /|  \n
            4------* |  \n
            | 3----|-*  \n
            |/     |/   \n
            0---8---1    \n
        */
        
        n.resize( 10 );
        n[ 0 ] = nodes_[ 0 ]; n[ 1 ] = nodes_[ 1 ]; n[ 2 ] = nodes_[ 3 ]; n[ 3 ] = nodes_[ 4 ];
        
        //*! VTK,Flaherty,Gimli count: 1-2-3-4, 5(1-2), 6(2-3), 7(3-1), 8(1-4), 9(2-4), 10(3-4)* //
        n[ 4 ] = nodes_[ 8 ];  n[ 5 ] = nodes_[ 20 ]; n[ 6 ] = nodes_[ 11 ];
        n[ 7 ] = nodes_[ 16 ]; n[ 8 ] = nodes_[ 22 ]; n[ 9 ] = nodes_[ 23 ];
        
        GIMLI::Tetrahedron10 tet10_b( n );
        testStiffness( tet10_b );
        
        CPPUNIT_ASSERT( tet10_b.createShapeFunctions()[ 0 ] == T4_1 * ( 2.0 * T4_1 + -1.0 ) );
        CPPUNIT_ASSERT( tet10_b.createShapeFunctions()[ 1 ] == T4_2 * ( 2.0 * T4_2 + -1.0 ) );
        CPPUNIT_ASSERT( tet10_b.createShapeFunctions()[ 2 ] == T4_3 * ( 2.0 * T4_3 + -1.0 ) );
        CPPUNIT_ASSERT( tet10_b.createShapeFunctions()[ 3 ] == T4_4 * ( 2.0 * T4_4 + -1.0 ) );
        CPPUNIT_ASSERT( tet10_b.createShapeFunctions()[ 4 ] == 4.0 * T4_1 * T4_2 );
        CPPUNIT_ASSERT( tet10_b.createShapeFunctions()[ 5 ] == 4.0 * T4_2 * T4_3 );
        CPPUNIT_ASSERT( tet10_b.createShapeFunctions()[ 6 ] == 4.0 * T4_3 * T4_1 );
        CPPUNIT_ASSERT( tet10_b.createShapeFunctions()[ 7 ] == 4.0 * T4_1 * T4_4 );
        CPPUNIT_ASSERT( tet10_b.createShapeFunctions()[ 8 ] == 4.0 * T4_2 * T4_4 );
        CPPUNIT_ASSERT( tet10_b.createShapeFunctions()[ 9 ] == 4.0 * T4_3 * T4_4 );
        
        //*! Zienkiewicz count: 1-2-3-4, 5(1-2), 6(1-3), 7(1-4), 8(2-3), 9(3-4), 10(4-2)* //
        n[ 4 ] = nodes_[ 8 ];  n[ 5 ] = nodes_[ 11 ]; n[ 6 ] = nodes_[ 16 ];
        n[ 7 ] = nodes_[ 20 ]; n[ 8 ] = nodes_[ 23 ]; n[ 9 ] = nodes_[ 22 ];

        GIMLI::Tetrahedron10 tet10( n );

        /*! Clear the cache to rebuild shapefunctions for the alternative tet10 numbering. */
        ShapeFunctionCache::instance().clear( );

        CPPUNIT_ASSERT( tet10.createShapeFunctions()[ 0 ] == T4_1 * ( 2.0 * T4_1 + -1.0 ) );
        CPPUNIT_ASSERT( tet10.createShapeFunctions()[ 1 ] == T4_2 * ( 2.0 * T4_2 + -1.0 ) );
        CPPUNIT_ASSERT( tet10.createShapeFunctions()[ 2 ] == T4_3 * ( 2.0 * T4_3 + -1.0 ) );
        CPPUNIT_ASSERT( tet10.createShapeFunctions()[ 3 ] == T4_4 * ( 2.0 * T4_4 + -1.0 ) );
        CPPUNIT_ASSERT( tet10.createShapeFunctions()[ 4 ] == 4.0 * T4_1 * T4_2 );
        CPPUNIT_ASSERT( tet10.createShapeFunctions()[ 5 ] == 4.0 * T4_1 * T4_3 );
        CPPUNIT_ASSERT( tet10.createShapeFunctions()[ 6 ] == 4.0 * T4_1 * T4_4 );
        CPPUNIT_ASSERT( tet10.createShapeFunctions()[ 7 ] == 4.0 * T4_2 * T4_3 );
        CPPUNIT_ASSERT( tet10.createShapeFunctions()[ 8 ] == 4.0 * T4_3 * T4_4 );
        CPPUNIT_ASSERT( tet10.createShapeFunctions()[ 9 ] == 4.0 * T4_4 * T4_2 );

        testStiffness( tet10 );
        
        n.resize( 8 );

        n[ 0 ] = nodes_[ 0 ];  n[ 1 ] = nodes_[ 1 ]; n[ 2 ] = nodes_[ 2 ]; n[ 3 ] = nodes_[ 3 ];
        n[ 4 ] = nodes_[ 4 ];  n[ 5 ] = nodes_[ 5 ]; n[ 6 ] = nodes_[ 6 ]; n[ 7 ] = nodes_[ 7 ];
        
        GIMLI::Hexahedron hex( n );
        
        GIMLI::RPolynomialFunction E2_2R( e2, RVector( 0 ), RVector( 0 ) );
        GIMLI::RPolynomialFunction E2_1R = -( -1.0 + E2_2R );
        GIMLI::RPolynomialFunction E2_2S( RVector( 0 ), e2, RVector( 0 ) );
        GIMLI::RPolynomialFunction E2_1S = -( -1.0 + E2_2S );
        GIMLI::RPolynomialFunction E2_2T( RVector( 0 ), RVector( 0 ), e2 );
        GIMLI::RPolynomialFunction E2_1T = -( -1.0 + E2_2T );
        
        GIMLI::RPolynomialFunction Q4_1 = E2_1R * E2_1S;
        GIMLI::RPolynomialFunction Q4_2 = E2_2R * E2_1S;
        GIMLI::RPolynomialFunction Q4_3 = E2_2R * E2_2S;
        GIMLI::RPolynomialFunction Q4_4 = E2_1R * E2_2S;
        
        CPPUNIT_ASSERT( E2_1R * E2_1S * E2_1T == Q4_1 * E2_1T );
        CPPUNIT_ASSERT( E2_2R * E2_1S * E2_1T == Q4_2 * E2_1T );
        CPPUNIT_ASSERT( E2_2R * E2_2S * E2_1T == Q4_3 * E2_1T );
        CPPUNIT_ASSERT( E2_1R * E2_2S * E2_1T == Q4_4 * E2_1T );
        CPPUNIT_ASSERT( E2_1R * E2_1S * E2_2T == Q4_1 * E2_2T );
        CPPUNIT_ASSERT( E2_2R * E2_1S * E2_2T == Q4_2 * E2_2T );
        CPPUNIT_ASSERT( E2_2R * E2_2S * E2_2T == Q4_3 * E2_2T );
        CPPUNIT_ASSERT( E2_1R * E2_2S * E2_2T == Q4_4 * E2_2T );
        
        CPPUNIT_ASSERT( hex.createShapeFunctions()[ 0 ] == Q4_1 * E2_1T );
        CPPUNIT_ASSERT( hex.createShapeFunctions()[ 1 ] == Q4_2 * E2_1T );
        CPPUNIT_ASSERT( hex.createShapeFunctions()[ 2 ] == Q4_3 * E2_1T );
        CPPUNIT_ASSERT( hex.createShapeFunctions()[ 3 ] == Q4_4 * E2_1T );
        CPPUNIT_ASSERT( hex.createShapeFunctions()[ 4 ] == Q4_1 * E2_2T );
        CPPUNIT_ASSERT( hex.createShapeFunctions()[ 5 ] == Q4_2 * E2_2T );
        CPPUNIT_ASSERT( hex.createShapeFunctions()[ 6 ] == Q4_3 * E2_2T );
        CPPUNIT_ASSERT( hex.createShapeFunctions()[ 7 ] == Q4_4 * E2_2T );
        
        testStiffness( hex );

        /*        7------*  \n
                 /|     /|  \n
                4------5 |  \n
                | 3----|-*  \n
                |/     |/   \n
                0------1    \n
            */
        n.resize( 6 );

        n[ 0 ] = nodes_[ 0 ];  n[ 1 ] = nodes_[ 1 ]; n[ 2 ] = nodes_[ 3 ]; 
        n[ 3 ] = nodes_[ 4 ];  n[ 4 ] = nodes_[ 5 ]; n[ 5 ] = nodes_[ 7 ];
        
        GIMLI::TriPrism pri( n );
        
        GIMLI::RPolynomialFunction T3_2( e2, RVector( 0 ) );
        GIMLI::RPolynomialFunction T3_3( RVector( 0 ), e2 );
        GIMLI::RPolynomialFunction T3_1 = -( -1.0 + T3_2 + T3_3 );
        
        GIMLI::RPolynomialFunction P6_1 = T3_1 * E2_1T;
        GIMLI::RPolynomialFunction P6_2 = T3_2 * E2_1T;
        GIMLI::RPolynomialFunction P6_3 = T3_3 * E2_1T;
        GIMLI::RPolynomialFunction P6_4 = T3_1 * E2_2T;
        GIMLI::RPolynomialFunction P6_5 = T3_2 * E2_2T;
        GIMLI::RPolynomialFunction P6_6 = T3_3 * E2_2T;
        
        
        CPPUNIT_ASSERT( pri.createShapeFunctions()[ 0 ] == P6_1 );
        CPPUNIT_ASSERT( pri.createShapeFunctions()[ 1 ] == P6_2 );
        CPPUNIT_ASSERT( pri.createShapeFunctions()[ 2 ] == P6_3 );
        CPPUNIT_ASSERT( pri.createShapeFunctions()[ 3 ] == P6_4 );
        CPPUNIT_ASSERT( pri.createShapeFunctions()[ 4 ] == P6_5 );
        CPPUNIT_ASSERT( pri.createShapeFunctions()[ 5 ] == P6_6 );
                
        testStiffness( pri );
        
        n.resize( 15 );        
        n[ 0 ] = nodes_[ 0 ];  n[ 1 ] = nodes_[ 1 ]; n[ 2 ] = nodes_[ 3 ]; 
        n[ 3 ] = nodes_[ 4 ];  n[ 4 ] = nodes_[ 5 ]; n[ 5 ] = nodes_[ 7 ];
        
        n[ 6 ] = nodes_[ 8  ];  n[ 7 ]  = nodes_[ 20 ]; n[ 8 ]  = nodes_[ 11 ];
        n[ 9 ] = nodes_[ 12 ];  n[ 10 ] = nodes_[ 21 ]; n[ 11 ] = nodes_[ 15 ];
        n[ 12 ] = nodes_[ 16 ]; n[ 13 ] = nodes_[ 17 ]; n[ 14 ] = nodes_[ 19 ]; 
        
        GIMLI::TriPrism15 pri2( n );

//         for (uint i = 0; i < pri2.nodeCount(); i ++ ){
//             std::cout << pri2.rst(i) << " " << pri2.node(i).pos() << std::endl;
//         }
//         std::cout << pri2.createShapeFunctions()[ 0 ] << std::endl;
//         std::cout << T3_1 * E2_1T * ( 2. * ( T3_1 + E2_1T ) + -3. ) << std::endl;
//         
        
        CPPUNIT_ASSERT( pri2.createShapeFunctions()[ 0 ] == T3_1 * E2_1T * ( 2. * ( T3_1 + E2_1T ) + -3. ) );
        CPPUNIT_ASSERT( pri2.createShapeFunctions()[ 1 ] == T3_2 * E2_1T * ( 2. * ( T3_2 + E2_1T ) + -3. ) );
        CPPUNIT_ASSERT( pri2.createShapeFunctions()[ 2 ] == T3_3 * E2_1T * ( 2. * ( T3_3 + E2_1T ) + -3. ) );
        CPPUNIT_ASSERT( pri2.createShapeFunctions()[ 3 ] == T3_1 * E2_2T * ( 2. * ( T3_1 + E2_2T ) + -3. ) );
        CPPUNIT_ASSERT( pri2.createShapeFunctions()[ 4 ] == T3_2 * E2_2T * ( 2. * ( T3_2 + E2_2T ) + -3. ) );
        CPPUNIT_ASSERT( pri2.createShapeFunctions()[ 5 ] == T3_3 * E2_2T * ( 2. * ( T3_3 + E2_2T ) + -3. ) );
        
        CPPUNIT_ASSERT( pri2.createShapeFunctions()[ 6 ] == 4.0 * T3_1 * T3_2 * E2_1T ); // T6_4 * E2_1T
        CPPUNIT_ASSERT( pri2.createShapeFunctions()[ 7 ] == 4.0 * T3_2 * T3_3 * E2_1T ); // T6_5 * E2_1T
        CPPUNIT_ASSERT( pri2.createShapeFunctions()[ 8 ] == 4.0 * T3_3 * T3_1 * E2_1T ); // T6_6 * E2_1T
        CPPUNIT_ASSERT( pri2.createShapeFunctions()[ 9 ] == 4.0 * T3_1 * T3_2 * E2_2T ); // T6_4 * E2_2T
        CPPUNIT_ASSERT( pri2.createShapeFunctions()[ 10 ] == 4.0 * T3_2 * T3_3 * E2_2T ); // T6_5 * E2_2T
        CPPUNIT_ASSERT( pri2.createShapeFunctions()[ 11 ] == 4.0 * T3_3 * T3_1 * E2_2T ); // T6_6 * E2_2T
        
        CPPUNIT_ASSERT( pri2.createShapeFunctions()[ 12 ] == 4.0 * T3_1 * E2_1T * E2_2T ); // T3_1 * E3_3T
        CPPUNIT_ASSERT( pri2.createShapeFunctions()[ 13 ] == 4.0 * T3_2 * E2_1T * E2_2T ); // T3_2 * E3_3T
        CPPUNIT_ASSERT( pri2.createShapeFunctions()[ 14 ] == 4.0 * T3_3 * E2_1T * E2_2T ); // T3_3 * E3_3T
        
        testStiffness( pri2 );
    }

    void testStiffness( const Cell & ent ){
//         std::cout.precision( 14 );
//         std::cout << ent.rtti() << " " << ent.shape().name() << std::endl;
        
        //** check if shapefunctions forms unity matrix
        for ( size_t i = 0; i < ent.nodeCount(); i ++ ){
                       
//             std::cout << ent.rst( i ) << std::endl;
                       
            RVector N( ent.N( ent.rst( i ) ) );
                        
            for ( size_t j = 0; j < ent.nodeCount(); j ++ ){
//                 std::cout << ent.createShapeFunctions()[j] << " " << ent.createShapeFunctions()[j]( ent.rst( i ) ) << std::endl;
//                 std::cout << std::setw( 5 ) << N[ j ];
                if ( i == j ) CPPUNIT_ASSERT( ( N[ j ] - 1.0) <  TOLERANCE );
                //else CPPUNIT_ASSERT( ::fabs( N[ j ] ) <  TOLERANCE );
            }
//              std::cout <<  std::endl;

            
        }
        
        ElementMatrix < double > S;
        S.ux2uy2uz2( ent );
        for ( size_t i = 0; i < S.size(); i ++ ){
            double s = sum( S.row( i ) );
            if ( ::fabs( s ) > TOLERANCE ) {
                std::cout << ent.rtti() << " " << ent.shape().name() << std::endl;
                std::cout << S.row( i ) << std::endl;
                std::cout << s << std::endl;
            }
            CPPUNIT_ASSERT( ::fabs( s ) < TOLERANCE );
        }
        // check derivatives sum to 0.\n";
        testMassElement( ent );
        testLoadElement( ent );
    }

    void testMassElement( const Cell & ent ){
        ElementMatrix < double > S;
        S.u2( ent );
        double s = 0.0;
        for ( size_t i = 0; i < S.size(); i ++ ){
            s += sum( S.row( i ) );
        }
        CPPUNIT_ASSERT( ::fabs( s - ent.shape().domainSize() ) < TOLERANCE );
    }
    
    void testLoadElement( const Cell & ent ){
        ElementMatrix < double > S;
        S.u( ent );
        double s = sum( S.row( 0 ) );
        //std::cout << s << " " << ent.shape().domainSize() << std::endl;
        CPPUNIT_ASSERT( ::fabs( s - ent.shape().domainSize() ) < TOLERANCE );
        
        //std::cout << ent << " size:" << ent.shape().domainSize() << std::endl << S << std::endl;
    }

    void resetNodes(){
        tearDown();
        setUp();
    }

    void setUp( ){
/*    7------6  \n
     /|     /|  \n
    4------5 |  \n
    | 3----|-2  \n
    |/     |/   \n
    0------1    \n
*/
        nodes_.resize( 24 );
        nodes_[ 0 ] = new GIMLI::Node( 0.0, 0.0, 0.0 );
        nodes_[ 1 ] = new GIMLI::Node( 1.0, 0.0, 0.0 );
        nodes_[ 2 ] = new GIMLI::Node( 1.0, 1.0, 0.0 );
        nodes_[ 3 ] = new GIMLI::Node( 0.0, 1.0, 0.0 );
        nodes_[ 4 ] = new GIMLI::Node( 0.0, 0.0, 1.0 );
        nodes_[ 5 ] = new GIMLI::Node( 1.0, 0.0, 1.0 );
        nodes_[ 6 ] = new GIMLI::Node( 1.0, 1.0, 1.0 );
        nodes_[ 7 ] = new GIMLI::Node( 0.0, 1.0, 1.0 );

        nodes_[ 8  ] = new GIMLI::Node( ( nodes_[ 0 ]->pos() + nodes_[ 1 ]->pos() ) / 2.0 );
        nodes_[ 9  ] = new GIMLI::Node( ( nodes_[ 1 ]->pos() + nodes_[ 2 ]->pos() ) / 2.0 );
        nodes_[ 10 ] = new GIMLI::Node( ( nodes_[ 2 ]->pos() + nodes_[ 3 ]->pos() ) / 2.0 );
        nodes_[ 11 ] = new GIMLI::Node( ( nodes_[ 3 ]->pos() + nodes_[ 0 ]->pos() ) / 2.0 );
        
        nodes_[ 12 ] = new GIMLI::Node( ( nodes_[ 4 ]->pos() + nodes_[ 5 ]->pos() ) / 2.0 );
        nodes_[ 13 ] = new GIMLI::Node( ( nodes_[ 5 ]->pos() + nodes_[ 6 ]->pos() ) / 2.0 );
        nodes_[ 14 ] = new GIMLI::Node( ( nodes_[ 6 ]->pos() + nodes_[ 7 ]->pos() ) / 2.0 );
        nodes_[ 15 ] = new GIMLI::Node( ( nodes_[ 7 ]->pos() + nodes_[ 4 ]->pos() ) / 2.0 );

        nodes_[ 16 ] = new GIMLI::Node( ( nodes_[ 0 ]->pos() + nodes_[ 4 ]->pos() ) / 2.0 );
        nodes_[ 17 ] = new GIMLI::Node( ( nodes_[ 1 ]->pos() + nodes_[ 5 ]->pos() ) / 2.0 );
        nodes_[ 18 ] = new GIMLI::Node( ( nodes_[ 2 ]->pos() + nodes_[ 6 ]->pos() ) / 2.0 );
        nodes_[ 19 ] = new GIMLI::Node( ( nodes_[ 3 ]->pos() + nodes_[ 7 ]->pos() ) / 2.0 );

        nodes_[ 20 ] = new GIMLI::Node( ( nodes_[ 1 ]->pos() + nodes_[ 3 ]->pos() ) / 2.0 );
        nodes_[ 21 ] = new GIMLI::Node( ( nodes_[ 5 ]->pos() + nodes_[ 7 ]->pos() ) / 2.0 );
        nodes_[ 22 ] = new GIMLI::Node( ( nodes_[ 1 ]->pos() + nodes_[ 4 ]->pos() ) / 2.0 );
        nodes_[ 23 ] = new GIMLI::Node( ( nodes_[ 3 ]->pos() + nodes_[ 4 ]->pos() ) / 2.0 );
        
        for ( size_t i = 0; i < nodes_.size(); i ++ ) nodes_[ i ]->setId( i );
    }

    void tearDown(){
        for ( size_t i = 0; i < nodes_.size(); i ++ ) delete nodes_[ i ];
        nodes_.clear();
    }

private:
    std::vector < Node * > nodes_;
    //* n1_, * n2_, * n3_, * n4_, * n5_, * n6_, * n7_, * n8_;
};

CPPUNIT_TEST_SUITE_REGISTRATION( FEMTest );