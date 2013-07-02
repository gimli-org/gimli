#include <cppunit/extensions/HelperMacros.h>

#include <gimli.h>
#include <node.h>
#include <shape.h>
#include <pos.h>
#include <meshentities.h>
#include <mesh.h>
#include <meshgenerators.h>

class ShapeTest : public CppUnit::TestFixture  {
    CPPUNIT_TEST_SUITE( ShapeTest );
    // CPPUNIT_TEST( testShapeFunctions );
    // CPPUNIT_TEST( testEquality );
    // CPPUNIT_TEST( testRSTXYZ );
    // CPPUNIT_TEST( testDomainSizes );
    // CPPUNIT_TEST( testJacobiDeterminat );
    // CPPUNIT_TEST( testTouch );
    // CPPUNIT_TEST( testTouch1 );
    // CPPUNIT_TEST( testInterpolate );
    // CPPUNIT_TEST( testSplit );
    // CPPUNIT_TEST( testGridGen );
    CPPUNIT_TEST_SUITE_END();

public:
    void testEquality(){
        CPPUNIT_ASSERT( *n1_ == *n1_ );
        CPPUNIT_ASSERT( !( *n1_ == *n2_ ) );
    }
    
    void testShapeFunctions(){

        RVector e1( 2 ); e1[ 0 ] = 1; e1[ 1 ] = -1; // 1-x
        RVector e2( 2 ); e2[ 0 ] = 0; e2[ 1 ] =  1; // x
        
        GIMLI::RPolynomialFunction E2_1_R( e1 );// 1-x
        GIMLI::RPolynomialFunction E2_2_R( e2 );// x
        GIMLI::RPolynomialFunction E2_1_S( RVector(0), e1 );// 1-y
        GIMLI::RPolynomialFunction E2_2_S( RVector(0), e2 );// y

        GIMLI::RPolynomialFunction Q4_1 = E2_1_R * E2_1_S ; // (1-x) * (1-y)    = 1 - x - y + xy
        GIMLI::RPolynomialFunction Q4_2 = E2_2_R * E2_1_S ; // x * (1-y)        = x - xy
        GIMLI::RPolynomialFunction Q4_3 = E2_2_R * E2_2_S ; // x * y            = xy
        GIMLI::RPolynomialFunction Q4_4 = E2_1_R * E2_2_S ; // (1-x) * y        = y - xy
        
        CPPUNIT_ASSERT( q1_->shape().createShapeFunctions()[ 0 ] == Q4_1 );
        CPPUNIT_ASSERT( q1_->shape().createShapeFunctions()[ 1 ] == Q4_2 );
        CPPUNIT_ASSERT( q1_->shape().createShapeFunctions()[ 2 ] == Q4_3 );
        CPPUNIT_ASSERT( q1_->shape().createShapeFunctions()[ 3 ] == Q4_4 );
    }
    
    void testRSTXYZ(){
        CPPUNIT_ASSERT( e1_->shape().xyz( RVector3( 0.0, 0.0, 0.0 ) ) == e1_->node( 0 ).pos() );
        CPPUNIT_ASSERT( e1_->shape().xyz( RVector3( 1.0, 0.0, 0.0 ) ) == e1_->node( 1 ).pos() );
        CPPUNIT_ASSERT( e1_->shape().xyz( RVector3( 0.5, 0.0, 0.0 ) ) == e1_->center() );
        
        CPPUNIT_ASSERT( t1_->shape().xyz( RVector3( 0.0, 0.0, 0.0 ) ) == t1_->node( 0 ).pos() );
        CPPUNIT_ASSERT( t1_->shape().xyz( RVector3( 1.0, 0.0, 0.0 ) ) == t1_->node( 1 ).pos() );
        CPPUNIT_ASSERT( t1_->shape().xyz( RVector3( 0.0, 1.0, 0.0 ) ) == t1_->node( 2 ).pos() );
        CPPUNIT_ASSERT( t1_->shape().xyz( RVector3( 1./3., 1./3.0, 0.0 ) ) == t1_->center() );
    }
    
    void testDomainSizes(){
        CPPUNIT_ASSERT( t1_->shape().domainSize() == 0.5 );
        CPPUNIT_ASSERT( t2_->shape().domainSize() == 0.5 );
        CPPUNIT_ASSERT( t3_->shape().domainSize() == 0.5 );
        CPPUNIT_ASSERT( q1_->shape().domainSize() == 1.0 );
        CPPUNIT_ASSERT( tet1_->shape().domainSize() == 1.0 / 6.0 );
        CPPUNIT_ASSERT( tet2_->shape().domainSize() == 1.0 / 6.0 );
        CPPUNIT_ASSERT( tet3_->shape().domainSize() == 1.0 / 6.0 );
        CPPUNIT_ASSERT( hex1_->shape().domainSize() == 1.0 );
        CPPUNIT_ASSERT( pri1_->shape().domainSize() == 0.5 );
    }
    void testJacobiDeterminat(){
        //std::cout << "testJacobiDeterminat(): " << det( inv( t1_->shape().createJacobian() ) ) << std::endl;
        CPPUNIT_ASSERT( ::fabs( det( inv( t1_->shape().createJacobian() ) ) - 1.0 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( det( inv( t2_->shape().createJacobian() ) ) - 1.0 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( det( inv( t3_->shape().createJacobian() ) ) - 1.0 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( det( inv( tet1_->shape().createJacobian() ) ) - 1.0 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( det( inv( tet2_->shape().createJacobian() ) ) + 1.0 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( det( inv( tet3_->shape().createJacobian() ) ) + 1.0 ) < TOLERANCE );
    }
    void testTouch(){
        CPPUNIT_ASSERT( t1_->shape().isInside( t1_->shape().center(), false ) );
        CPPUNIT_ASSERT( t1_->shape().isInside( GIMLI::RVector3( 0.5, 0.25 ) ) );
        CPPUNIT_ASSERT( t1_->shape().isInside( GIMLI::RVector3( 0.5, 0.5 ) ) );
        CPPUNIT_ASSERT( t1_->shape().isInside( GIMLI::RVector3( 1.0, 1.0 ) ) );
        CPPUNIT_ASSERT( t1_->shape().isInside( GIMLI::RVector3( 0.0, 0.0 ) ) );
        CPPUNIT_ASSERT( !t1_->shape().isInside( GIMLI::RVector3( 0.5, 0.75 ) ) );
        CPPUNIT_ASSERT( !t1_->shape().isInside( GIMLI::RVector3( 2.0, -1.0 ) ) );

        CPPUNIT_ASSERT( tet1_->shape().isInside( tet1_->shape().center() ) );
    }
    void testTouch1(){
//         RVector sf;
//         t1_->shape().isInside( GIMLI::RVector3( -0.1,  0.0 ), sf, false ); CPPUNIT_ASSERT( pIndx == 1 );
//         t1_->shape().isInside( GIMLI::RVector3(  0.5,  1.0 ), sf, false ); CPPUNIT_ASSERT( pIndx == 1 );
//         t1_->shape().isInside( GIMLI::RVector3(  0.5, -1.0 ), sf, false ); CPPUNIT_ASSERT( pIndx == 2 );
// 
//         GIMLI::Node n1( 5.0, -1.0 ); n1.setId( 0 );
//         GIMLI::Node n2( 6.0, -1.0 ); n2.setId( 1 );
//         GIMLI::Node n3( 6.0,  0.0 ); n3.setId( 2 );
//         GIMLI::Node n4( 5.0,  0.0 ); n4.setId( 3 );
//         GIMLI::Quadrangle q( n1, n2, n3, n4 );
// 
//         q.shape().isInside( GIMLI::RVector3(  4.0,  0.0 ), pIndx, false ); CPPUNIT_ASSERT( pIndx == 2 );
//         q.shape().isInside( GIMLI::RVector3(  4.0, -0.5 ), pIndx, false ); CPPUNIT_ASSERT( pIndx == 1 );
//         q.shape().isInside( GIMLI::RVector3(  4.0, -1.0 ), pIndx, false ); CPPUNIT_ASSERT( pIndx == 1 );
//         q.shape().isInside( GIMLI::RVector3(  7.0, -0.5 ), pIndx, false ); CPPUNIT_ASSERT( pIndx == 0 );
//         q.shape().isInside( GIMLI::RVector3(  7.0, -1.0 ), pIndx, false ); CPPUNIT_ASSERT( pIndx == 0 );
//         q.shape().isInside( GIMLI::RVector3(  5.0, -2.0 ), pIndx, false ); CPPUNIT_ASSERT( pIndx == 3 );
//         q.shape().isInside( GIMLI::RVector3(  5.5, -2.0 ), pIndx, false ); CPPUNIT_ASSERT( pIndx == 2 );
//         q.shape().isInside( GIMLI::RVector3(  6.0, -2.0 ), pIndx, false ); CPPUNIT_ASSERT( pIndx == 2 );
//         q.shape().isInside( GIMLI::RVector3(  5.0, 1.0 ), pIndx, false ); CPPUNIT_ASSERT( pIndx == 0 );
//         q.shape().isInside( GIMLI::RVector3(  5.5, 1.0 ), pIndx, false ); CPPUNIT_ASSERT( pIndx == 0 );
//         q.shape().isInside( GIMLI::RVector3(  6.0, 1.0 ), pIndx, false ); CPPUNIT_ASSERT( pIndx == 1 );
    }

    void testInterpolate(){
        GIMLI::Triangle tri( *n1_, *n2_, *n4_ );

        RVector u( 5, 0.0 );
        u[ 0 ] = 1.0;
        CPPUNIT_ASSERT( ::fabs( tri.shape().rst( tri.center() )[ 0 ] - 1.0/3.0 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( tri.shape().rst( tri.center() )[ 1 ] - 1.0/3.0 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( tri.pot( GIMLI::RVector3( 0.0, 0.5 ), u) - 0.5 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( tri.pot( GIMLI::RVector3( 0.5, 0.0 ), u) - 0.5 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( tri.pot( GIMLI::RVector3( 0.5, 0.5 ), u) - 0.0 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( tri.grad( GIMLI::RVector3( 0, 0 ), u).dist( GIMLI::RVector3( -1.0, -1.0 ) ) ) < TOLERANCE );

        GIMLI::Triangle tri2( *n1_, *n3_, *n4_ );
        CPPUNIT_ASSERT( ::fabs( tri2.shape().rst( tri2.center() )[ 0 ] - 1.0/3.0 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( tri2.shape().rst( tri2.center() )[ 1 ] - 1.0/3.0 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( tri2.pot( GIMLI::RVector3( 0.0, 0.5 ), u) - 0.5 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( tri2.pot( GIMLI::RVector3( 1.0, 1.0 ), u) - 0.0 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( tri2.pot( GIMLI::RVector3( 0.5, 0.5 ), u) - 0.5 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( tri2.grad( GIMLI::RVector3( 0, 0 ), u).dist( GIMLI::RVector3( 0.0, -1.0 ) ) ) < TOLERANCE );

        GIMLI::Tetrahedron tet( *n1_, *n2_, *n4_, *n5_ );
        CPPUNIT_ASSERT( ::fabs( tet.shape().rst( tet.center() )[ 0 ] - 0.25 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( tet.shape().rst( tet.center() )[ 1 ] - 0.25 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( tet.shape().rst( tet.center() )[ 2 ] - 0.25 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( tet.pot( GIMLI::RVector3( 0.0, 0.5 ), u ) - 0.5 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( tet.pot( GIMLI::RVector3( 0.5, 0.5 ), u ) - 0.0 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( tet.pot( GIMLI::RVector3( 0.5, 0.0 ), u ) - 0.5 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( tet.pot( GIMLI::RVector3( 0.0, 0.0, 0.5 ), u ) - 0.5 ) < TOLERANCE );
        CPPUNIT_ASSERT( ::fabs( tet.grad( GIMLI::RVector3( 0, 0 ), u).dist( GIMLI::RVector3( -1.0, -1.0, -1.0 ) ) ) < TOLERANCE );
    }

    void testSplit(){
        
        double testSum = 0.0;
        for ( size_t j = 0; j < 5; j ++ ){
            GIMLI::Tetrahedron tet( hex1_->node( GIMLI::HexahedronSplit5TetID[ j ][ 0 ] ),
                                    hex1_->node( GIMLI::HexahedronSplit5TetID[ j ][ 1 ] ),
                                    hex1_->node( GIMLI::HexahedronSplit5TetID[ j ][ 2 ] ),
                                    hex1_->node( GIMLI::HexahedronSplit5TetID[ j ][ 3 ] ) );
            testSum += tet.shape().domainSize();
            CPPUNIT_ASSERT( det( inv( tet.shape().createJacobian() ) ) > 1e-12 );
        }
        CPPUNIT_ASSERT( std::fabs( testSum -1.0 ) < 1e-12 );

        testSum = 0.0;
        for ( size_t j = 0; j < 6; j ++ ){
            //std::cout << j << std::endl;
            GIMLI::Tetrahedron tet( hex1_->node( GIMLI::HexahedronSplit6TetID[ j ][ 0 ] ),
                                    hex1_->node( GIMLI::HexahedronSplit6TetID[ j ][ 1 ] ),
                                    hex1_->node( GIMLI::HexahedronSplit6TetID[ j ][ 2 ] ),
                                    hex1_->node( GIMLI::HexahedronSplit6TetID[ j ][ 3 ] ) );
//             std::cout << hex1_->node( GIMLI::HexahedronSplit6TetID[ j ][ 0 ])  << "\n " << hex1_->node( GIMLI::HexahedronSplit6TetID[ j ][ 1 ] ) <<"\n "
//                     << hex1_->node( GIMLI::HexahedronSplit6TetID[ j ][ 2 ] )<< "\n" << hex1_->node( GIMLI::HexahedronSplit6TetID[ j ][ 3 ] ) <<std::endl;
            testSum += tet.shape().domainSize();

//             std::cout << j << " " <<  tet << std::endl;
//             std::cout << j << " " <<  tet.shape().jacobianDeterminant() << std::endl;
            CPPUNIT_ASSERT( det( inv( tet.shape().createJacobian() ) ) > 1e-12 );
        }
        CPPUNIT_ASSERT( std::fabs( testSum -1.0 ) < 1e-12 );
    }

    void setUp( ){
/*    8------7  \n
     /|     /|  \n
    5------6 |  \n
    | 4----|-3  \n
    |/     |/   \n
    1------2    \n
*/
        n1_ = new GIMLI::Node( 0.0, 0.0 ); n1_->setId( 0 );
        n2_ = new GIMLI::Node( 1.0, 0.0 ); n2_->setId( 1 );
        n3_ = new GIMLI::Node( 1.0, 1.0 ); n3_->setId( 2 );
        n4_ = new GIMLI::Node( 0.0, 1.0 ); n4_->setId( 3 );
        n5_ = new GIMLI::Node( 0.0, 0.0, 1.0 ); n5_->setId( 4 );
        n6_ = new GIMLI::Node( 1.0, 0.0, 1.0 ); n6_->setId( 5 );
        n7_ = new GIMLI::Node( 1.0, 1.0, 1.0 ); n7_->setId( 6 );
        n8_ = new GIMLI::Node( 0.0, 1.0, 1.0 ); n8_->setId( 7 );

        e1_ = new GIMLI::Edge( *n1_, *n2_ );
        
        t1_ = new GIMLI::Triangle( *n1_, *n2_, *n3_ );
        t2_ = new GIMLI::Triangle( *n3_, *n1_, *n2_ );
        t3_ = new GIMLI::Triangle( *n1_, *n3_, *n4_ );

        q1_ = new GIMLI::Quadrangle( *n1_, *n2_, *n3_, *n4_);

        tet1_ = new GIMLI::Tetrahedron( *n1_, *n2_, *n4_, *n5_);
        tet2_ = new GIMLI::Tetrahedron( *n5_, *n1_, *n2_, *n4_ );
        tet3_ = new GIMLI::Tetrahedron( *n1_, *n5_, *n4_, *n2_ );

        GIMLI::Node *nodes[]={ n1_, n2_, n3_, n4_, n5_, n6_, n7_, n8_ };
        std::vector < GIMLI::Node* > n(8); std::copy( &nodes[0], &nodes[8], &n[0] );
        hex1_ = new GIMLI::Hexahedron( n );
        
        GIMLI::Node *nnpri[]={ n1_, n2_, n4_, n5_, n6_, n8_ };
        n.resize( 6 ); std::copy( &nnpri[0], &nnpri[6], &n[0] );
        pri1_ = new GIMLI::TriPrism( n );
    }

    void testGridGen(){
        RVector x( 2 ); x[0]=0; x[1]=1;
        RVector y( 2 ); y[0]=0; y[1]=1;
        Mesh mesh( createMesh2D( x, y ) );
        for ( size_t i = 0; i < mesh.cellCount(); i ++ ){
            CPPUNIT_ASSERT( det( inv( mesh.cell( i ).shape().createJacobian() ) ) > 1e-12 );
        }
        
        RVector x1( 2 ); x1[0]=-1; x1[1]=1;
        RVector y1( 2 ); y1[0]=-2; y1[1]=-1;
        Mesh mesh1( createMesh2D( x1, y1 ) );
        for ( size_t i = 0; i < mesh1.cellCount(); i ++ ){
            CPPUNIT_ASSERT( det( inv( mesh1.cell( i ).shape().createJacobian() ) ) > 1e-12 );
        }
        
    }
    
    void tearDown(){
        delete tet1_;
        delete tet2_;
        delete tet3_;

        delete hex1_;
        delete q1_;
        
        delete pri1_;

        delete e1_;
        
        delete t1_;
        delete t2_;
        delete t3_;

        delete n1_;
        delete n2_;
        delete n3_;
        delete n4_;
        delete n5_;
        delete n6_;
        delete n7_;
        delete n8_;
    }

private:
    GIMLI::Node         * n1_, * n2_, * n3_, * n4_, * n5_, * n6_, * n7_, * n8_;

    GIMLI::Edge         * e1_;
    
    GIMLI::Triangle     * t1_, * t2_, * t3_;
    GIMLI::Quadrangle   * q1_;
    
    GIMLI::Tetrahedron  * tet1_, * tet2_, * tet3_ ;
    GIMLI::Hexahedron   * hex1_;
    
    GIMLI::TriPrism     * pri1_;

};

CPPUNIT_TEST_SUITE_REGISTRATION( ShapeTest );