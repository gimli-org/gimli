#include <cppunit/extensions/HelperMacros.h>

#include <gimli.h>
#include <node.h>
#include <shape.h>
#include <pos.h>
#include <meshentities.h>
#include <elementmatrix.h>

class FEMTest : public CppUnit::TestFixture  {
    CPPUNIT_TEST_SUITE( FEMTest );
    CPPUNIT_TEST( testFEM1D );
    CPPUNIT_TEST( testFEM2D );
    CPPUNIT_TEST( testFEM3D );

    CPPUNIT_TEST_SUITE_END();

public:

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
        GIMLI::EdgeCell edg( n );
        testStiffness( edg );

        n.resize( 3 );
        n[ 0 ] = nodes_[ 0 ]; n[ 1 ] = nodes_[ 1 ]; n[ 2 ] = nodes_[ 8 ];
        GIMLI::Edge3Cell edg3( n );
        testStiffness( edg3 );
    }

    void testStiffness2D(){
        std::vector < Node * > n( 2 );

        n.resize( 3 );
        n[ 0 ] = nodes_[ 0 ]; n[ 1 ] = nodes_[ 1 ]; n[ 2 ] = nodes_[ 2 ];
        GIMLI::Triangle tri( n );
        testStiffness( tri );

        n.resize( 6 );
        n[ 0 ] = nodes_[ 0 ]; n[ 1 ] = nodes_[ 1 ]; n[ 2 ] = nodes_[ 3 ];
        n[ 3 ] = nodes_[ 8 ]; n[ 4 ] = nodes_[ 12 ]; n[ 5 ] = nodes_[ 11 ];
        GIMLI::Triangle6 tri6( n );
        testStiffness( tri6 );

        n.resize( 4 );
        n[ 0 ] = nodes_[ 0 ]; n[ 1 ] = nodes_[ 1 ]; n[ 2 ] = nodes_[ 2 ];n[ 3 ] = nodes_[ 3 ];
        GIMLI::Quadrangle quad( n );
        testStiffness( quad );

        n.resize( 8 );
        n[ 0 ] = nodes_[ 0 ]; n[ 1 ] = nodes_[ 1 ]; n[ 2 ] = nodes_[ 2 ];  n[ 3 ] = nodes_[ 3 ];
        n[ 4 ] = nodes_[ 8 ]; n[ 5 ] = nodes_[ 9 ]; n[ 6 ] = nodes_[ 10 ]; n[ 7 ] = nodes_[ 11 ];
        GIMLI::Quadrangle8 quad8( n );
        testStiffness( quad8 );
    }

    void testStiffness3D(){
        std::vector < Node * > n( 4 );

        n.resize( 4 );
        n[ 0 ] = nodes_[ 0 ]; n[ 1 ] = nodes_[ 1 ]; n[ 2 ] = nodes_[ 3 ];n[ 3 ] = nodes_[ 4 ];
        GIMLI::Tetrahedron tet( n );
        testStiffness( tet );

        n.resize( 10 );
        n[ 0 ] = nodes_[ 0 ]; n[ 1 ] = nodes_[ 1 ]; n[ 2 ] = nodes_[ 3 ]; n[ 3 ] = nodes_[ 4 ];

        //*! Zienkiewicz count: 1-2-3-4, 5(1-2), 6(1-3), 7(1-4), 8(2-3), 9(3-4), 10(4-2)* //
        n[ 4 ] = nodes_[ 8 ];  n[ 5 ] = nodes_[ 11 ]; n[ 6 ] = nodes_[ 13 ];
        n[ 7 ] = nodes_[ 12 ]; n[ 8 ] = nodes_[ 17 ]; n[ 9 ] = nodes_[ 16 ];

        //*! VTK,Flaherty,Gimli count: 1-2-3-4, 5(1-2), 6(2-3), 7(3-1), 8(1-4), 9(2-4), 10(3-4)* //

        GIMLI::Tetrahedron10 tet10( n );
        testStiffness( tet10 );
    }

    void testStiffness( const Cell & ent ){
        ElementMatrix < double > S;
        S.ux2uy2uz2( ent );

        //** check if shapefunctions forms unity matrix
        for ( size_t i = 0; i < ent.nodeCount(); i ++ ){
            std::pair < std::vector < size_t >, RVector > N( ent.shapeFunctions( ent.node( i ).pos() ) );
            //std::cout << ent.node( i ).pos() << std::endl;
            for ( size_t j = 0; j < ent.nodeCount(); j ++ ){
                if ( i == j ) CPPUNIT_ASSERT( ( N.second[ j ] - 1.0) <  TOLERANCE );
                else CPPUNIT_ASSERT( ::fabs( N.second[ j ] ) <  TOLERANCE );
                //std::cout << std::setw( 5 ) << N.second[j];
            }
            //std::cout <<  std::endl;
        }

        for ( size_t i = 0; i < S.size(); i ++ ){
            double s = sum( S.row( i ) );
            if ( ::fabs( s ) > TOLERANCE ) {
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
        nodes_.resize( 18 );
        nodes_[ 0 ] = new GIMLI::Node( 0.0, 0.0 );
        nodes_[ 1 ] = new GIMLI::Node( 1.0, 0.0 );
        nodes_[ 2 ] = new GIMLI::Node( 1.0, 1.0 );
        nodes_[ 3 ] = new GIMLI::Node( 0.0, 1.0 );
        nodes_[ 4 ] = new GIMLI::Node( 0.0, 0.0, 1.0 );
        nodes_[ 5 ] = new GIMLI::Node( 1.0, 0.0, 1.0 );
        nodes_[ 6 ] = new GIMLI::Node( 1.0, 1.0, 1.0 );
        nodes_[ 7 ] = new GIMLI::Node( 0.0, 1.0, 1.0 );

        nodes_[ 8  ] = new GIMLI::Node( ( nodes_[ 0 ]->pos() + nodes_[ 1 ]->pos() ) / 2.0 );
        nodes_[ 9  ] = new GIMLI::Node( ( nodes_[ 1 ]->pos() + nodes_[ 2 ]->pos() ) / 2.0 );
        nodes_[ 10 ] = new GIMLI::Node( ( nodes_[ 2 ]->pos() + nodes_[ 3 ]->pos() ) / 2.0 );
        nodes_[ 11 ] = new GIMLI::Node( ( nodes_[ 3 ]->pos() + nodes_[ 0 ]->pos() ) / 2.0 );
        nodes_[ 12 ] = new GIMLI::Node( ( nodes_[ 1 ]->pos() + nodes_[ 3 ]->pos() ) / 2.0 );

        nodes_[ 13 ] = new GIMLI::Node( ( nodes_[ 0 ]->pos() + nodes_[ 4 ]->pos() ) / 2.0 );
        nodes_[ 14 ] = new GIMLI::Node( ( nodes_[ 4 ]->pos() + nodes_[ 7 ]->pos() ) / 2.0 );
        nodes_[ 15 ] = new GIMLI::Node( ( nodes_[ 7 ]->pos() + nodes_[ 3 ]->pos() ) / 2.0 );
        nodes_[ 16 ] = new GIMLI::Node( ( nodes_[ 1 ]->pos() + nodes_[ 4 ]->pos() ) / 2.0 );
        nodes_[ 17 ] = new GIMLI::Node( ( nodes_[ 3 ]->pos() + nodes_[ 4 ]->pos() ) / 2.0 );



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