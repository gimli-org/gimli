#include <cppunit/extensions/HelperMacros.h>

#include <gimli.h>
#include <pos.h>
#include <vector.h>
#include <matrix.h>
#include <vectortemplates.h>
#include <vector>

#include <stdexcept>

using namespace GIMLI;

class VectorTest : public CppUnit::TestFixture  {
    CPPUNIT_TEST_SUITE(VectorTest);
    CPPUNIT_TEST(testEquality);
    CPPUNIT_TEST(testFill);
    CPPUNIT_TEST(testSetVal);
    CPPUNIT_TEST(testUnaryOperations);
    CPPUNIT_TEST(testBinaryOperations);
    CPPUNIT_TEST(testFunctions);
    CPPUNIT_TEST(testStdVectorTemplates);
    CPPUNIT_TEST(testCVector);
    CPPUNIT_TEST(testRVector3);
    CPPUNIT_TEST(testMatrix);
    CPPUNIT_TEST(testFind);
    CPPUNIT_TEST(testIO);
    
    //CPPUNIT_TEST_EXCEPTION(funct, exception);
    CPPUNIT_TEST_SUITE_END();
    
public:    
    
    void testEquality(){
        testEquality_< RVector > ();
        testEquality_< CVector > ();
    }
 
    template < class Vec > void testEquality_(){
        Vec v1 = Vec(10);
        
        try{ v1.getVal(0); CPPUNIT_ASSERT(1); } catch(...){}
        try{ v1.getVal(10); CPPUNIT_ASSERT(0); } catch(...){}
        
        Vec v2 = Vec(v1);
        CPPUNIT_ASSERT(v1 == v1);
        CPPUNIT_ASSERT(v1 == v2);
            
        Vec v;          CPPUNIT_ASSERT(v1 != v);
        v = v1;         CPPUNIT_ASSERT(v1 == v);
        v = v1 + 1.0;   CPPUNIT_ASSERT(v1 != v);
        v = v1 + 1.0;   CPPUNIT_ASSERT(v == v1 + 1.0);
    }
     
    void testFill(){
        double vec[5]={0.0, 1.0, 2.0, 3.0, 4.0};
        
        RVector v1(5); v1.fill(vec);
        Placeholder x__;
        RVector v2(5); v2.fill(x__);
        CPPUNIT_ASSERT(v1 == v2);
        CPPUNIT_ASSERT(fliplr(v1)[ 0 ] == 4.0);
        CPPUNIT_ASSERT(fliplr(v1)[ 4 ] == 0.0);
        CPPUNIT_ASSERT(fliplr(fliplr(v1)) == v1);
        
    }
     
    void testUnaryOperations(){
        testUnaryOperations_< double > (3.1415);
        testUnaryOperations_< Complex > (Complex(3.1415, 1.0));
        RVector v(10, 42.42);
        RVector v1(sign(square(v))); CPPUNIT_ASSERT(sum(v1) == 10);    
        randn(v);
        CPPUNIT_ASSERT(RVector(abs(v + 1e-6)) == RVector(exp10(log10(abs(v + 1e-6)))));
        CPPUNIT_ASSERT(RVector(abs(v + 1e-6)) == RVector(exp(log(abs(v + 1e-6)))));
        CPPUNIT_ASSERT(RVector(abs(v + 1e-6)) != RVector(exp10(log(abs(v + 1e-6)))));   
    }
    template < class ValueType > void testUnaryOperations_(const ValueType & fac){
        typedef Vector < ValueType > Vec;
        Vec v1 = Vec(10);
        Vec v(v1), v2(v1);     
        v.fill(fac); 
        v2.fill(fac);     
        v += fac; CPPUNIT_ASSERT(v == Vec(v1.size(), fac * 2.0));
        v -= fac; CPPUNIT_ASSERT(v == Vec(v1.size(), fac));
        v *= fac; CPPUNIT_ASSERT(v == Vec(v1.size(), fac * fac));
        v /= fac; CPPUNIT_ASSERT(v == Vec(v1.size(), fac));    
        v += v2;  CPPUNIT_ASSERT(v == Vec(v1.size(), fac * 2.0));
        v -= v2;  CPPUNIT_ASSERT(v == Vec(v1.size(), fac));
        v *= v2;  CPPUNIT_ASSERT(v == Vec(v1.size(), fac * fac));
        v /= v2;  CPPUNIT_ASSERT(v == Vec(v1.size(), fac));
        v = -v2;  CPPUNIT_ASSERT(v == v2 * -1.);   
    }
    void testSetVal(){
        typedef Vector < double > Vec;
        Vec v1(10);
        v1.fill(x__ + 1.0);
        Vec v2(10, 1.0);
        Vec v3(10);

        CPPUNIT_ASSERT(sum(v3.setVal(v1, 0, 10)) == 55);
        CPPUNIT_ASSERT(sum(v3.addVal(v2, 5, 10)) == 60);
        CPPUNIT_ASSERT(sum(v3.setVal(1.0)) == 10);
        CPPUNIT_ASSERT(sum(v3.setVal(2.0, find(v1 > 5))) == 15);
        
        RVector v4(5, 3.0);
        CPPUNIT_ASSERT(sum(v3.setVal(v4, 5, 10)) == 20);
        CPPUNIT_ASSERT(sum(v3.addVal(v4, 5, 10)) == 35);
        
        // set from 0 to end == 1
        CPPUNIT_ASSERT(sum(v1.setVal(1.0, 0, -1)) == 10);
        // set from 5 to end == 0
        CPPUNIT_ASSERT(sum(v1.setVal(0.0, 5, 100)) == 5);
        // set from 5 to 7 == 1
        CPPUNIT_ASSERT(sum(v1.setVal(1.0, 5, 7)) == 7);
    }
    
    void testBinaryOperations(){
        RVector v1(*v1_);
        RVector v2(*v1_);          
        double fac = 3.1415;
        v1.fill(fac); v2.fill(fac);
        CPPUNIT_ASSERT(RVector(v1 + v2) == RVector(v1_->size(), 2.0 * fac));
        CPPUNIT_ASSERT(RVector(v1 - v2) == RVector(v1_->size(), 0.0));
        CPPUNIT_ASSERT(RVector(v2 - v1) == RVector(v1_->size(), 0.0));
        CPPUNIT_ASSERT(RVector(v1 * v2) == RVector(v1_->size(), fac * fac));
        CPPUNIT_ASSERT(RVector(v1 / v2) == RVector(v1_->size(), 1.0));
        CPPUNIT_ASSERT(RVector(v2 / v1) == RVector(v1_->size(), 1.0));
        
        v2 *= 2.0;
        CPPUNIT_ASSERT(RVector(v1 - v2) == RVector(v1_->size(), -fac));
        CPPUNIT_ASSERT(RVector(v2 - v1) == RVector(v1_->size(),  fac));
        CPPUNIT_ASSERT(RVector(v1 - v2) != RVector(v1_->size(),  fac));
        CPPUNIT_ASSERT(RVector(v1 / v2) == RVector(v1_->size(),  0.5));
        CPPUNIT_ASSERT(RVector(v2 / v1) == RVector(v1_->size(),  2.0));
    }
    
    void testFunctions(){
        typedef Vector < double > Vec;
         
        Vec v(*v1_);
        v.fill(x__ + 1.0);
        CPPUNIT_ASSERT(sum(v) == (v[ 0 ] + v[ v1_->size() -1 ]) * 
                                                (::floor(v1_->size() / 2.0) + 
                                                    (v1_->size() / 2.0 - ::floor(v1_->size() / 2.0))));
        CPPUNIT_ASSERT(GIMLI::sum(v - v) == 0.0);
        CPPUNIT_ASSERT(min(v) == 1.0);
        CPPUNIT_ASSERT(max(v) == v1_->size());
        RVector vs(10);
        vs[ 1 ] = 1; vs[ 2 ] = 1; vs[ 3 ] = 2; vs[ 4 ] = 2; vs[ 7 ] = 1;
        CPPUNIT_ASSERT(unique(vs).size() == 6); // 0 1 2 0 1 0
        CPPUNIT_ASSERT(unique(sort(vs)).size() == 3); // 0 1 2
        
        vs.fill(x__ + 1.0);
        vs *= 1.111;
        vs.round(0.01);
        CPPUNIT_ASSERT(::fabs(vs[0] - 1.11) < TOLERANCE);
		std::cout << vs[4] << std::endl;
        CPPUNIT_ASSERT(::fabs(vs[4] - 5.56) < TOLERANCE);
        CPPUNIT_ASSERT(::fabs(vs[6] - 7.78) < TOLERANCE);
        CPPUNIT_ASSERT(::fabs(vs[9] - 11.11) < TOLERANCE);
        
        vs *=0;
        CPPUNIT_ASSERT(nonZero(vs) == false);
        CPPUNIT_ASSERT(zero(vs) == true);
        
        
    }
        
    void testStdVectorTemplates(){
        std::vector < int > iVec(10, 0);
        iVec[ 1 ] = 1; iVec[ 2 ] = 1; iVec[ 3 ] = 2; iVec[ 4 ] = 2; iVec[ 7 ] = 1;
        CPPUNIT_ASSERT(unique(iVec).size() == 6); // 0 1 2 0 1 0
        CPPUNIT_ASSERT(unique(sort(iVec)).size() == 3); // 0 1 2
        
        std::vector < double > dVec(10, 1.0);
        CPPUNIT_ASSERT(GIMLI::max(dVec) == 1.0); 
    }
    
    void testRVector3(){
        RVector3 p1(1.0, 1.0, 1.0);
        RVector3 p2(-1.0, -1.0, -1.0);
        CPPUNIT_ASSERT(p1 == p1);
        CPPUNIT_ASSERT(p1 != p2);
        CPPUNIT_ASSERT(p1 == -p2);
        CPPUNIT_ASSERT(p1 != p2);
        CPPUNIT_ASSERT(p1 != p1 * -1.0);
        CPPUNIT_ASSERT(p1 == (p2 * -1.0));
        CPPUNIT_ASSERT(p1 == p1 / 1.0);
        CPPUNIT_ASSERT(p1 * -1.0 == p1 / -1.0);
        CPPUNIT_ASSERT(p1 + p2 == p2 + p1);
        p2 *= -1.0;
        CPPUNIT_ASSERT(p1 == p2);
        p2 = p1; p2 += p1;  CPPUNIT_ASSERT(p1 + p1 == p2);
        p2 = p1; p2 += 1.0;  CPPUNIT_ASSERT(p1 + p1 == p2);
        p2 = p1; p2.translate(p1);  CPPUNIT_ASSERT(p1 + p1 == p2);
        p2 = p1; p2.scale(p1);  CPPUNIT_ASSERT(p1 == p2);
    }
        
    void testCVector(){
        size_t dim = 10;
//         RVector freq(dim);
//         Complex i_unit(0.0 , 1.0);
//         CVector k1(dim) ;
//         CVector g(dim, Complex(1.0, 0.0)); 
//         CVector z = toComplex(freq) / (k1 * g) * i_unit;
        
        CVector c1(dim, Complex(1.0, 0.0)); 
        RVector r1(dim, 1.0);
        CVector c2;
        
        c2 = c1 * r1; 
        CPPUNIT_ASSERT(c1 == c2);
        c2 = r1 * c1; 
        CPPUNIT_ASSERT(c1 == c2);
        c2 = toComplex(real(c1), imag(c1)); 
        CPPUNIT_ASSERT(c1 == c2);
        CPPUNIT_ASSERT(real(c1 + conj(c1)) == RVector(2.0 * real(c1)));
        CPPUNIT_ASSERT(imag(c1 - conj(c1)) == RVector(2.0 * imag(c1)));
    
        r1 *= 0.0;
        CPPUNIT_ASSERT(imag(c1 * conj(c1)) == r1);
        CPPUNIT_ASSERT(RVector(angle(c1) + angle(conj(c1))) == r1);

        CPPUNIT_ASSERT(abs(c1) == RVector(sqrt(real(c1) * real(c1)
                                                  + imag(c1) * imag(c1))));
        
        CPPUNIT_ASSERT(abs(c1) == RVector(exp(real(log(c1)))));
        
        //real, imag, angle, conj, abs
    }
    
    template < class ValueType > void testMatrix_(){
        typedef Matrix < ValueType > Mat;
        typedef Vector < ValueType > Vec;
        
        Mat A(5, 5);
        try{ A.row(11); CPPUNIT_ASSERT(0); } catch(...){}
        try{ A.row(-1); CPPUNIT_ASSERT(0); } catch(...){}
        try{ CPPUNIT_ASSERT(0); } catch(...){}
        CPPUNIT_ASSERT(A.rows() == 5);
        CPPUNIT_ASSERT(A.cols() == 5);
                        
        A.resize(3, 2); CPPUNIT_ASSERT(A.cols() == 2);
        A.resize(8, 9); CPPUNIT_ASSERT(A.cols() == 9);
        
        A[ 0 ][ 0 ] = 1.0;
        A[ 1 ] = A[ 0 ];

        CPPUNIT_ASSERT(A[ 0 ] == A[ 1 ]);
        CPPUNIT_ASSERT(A.row(2) != A[ 1 ]);
        CPPUNIT_ASSERT(A[ 0 ][ 0 ]  == 1.0);
        CPPUNIT_ASSERT(A[ 1 ][ 0 ]  == 1.0);
        
        CPPUNIT_ASSERT(fliplr(fliplr(A)) == A);
        
        A.push_back(A[ 0 ]); CPPUNIT_ASSERT(A.rows() == 9);
        CPPUNIT_ASSERT(A[ A.rows()-1 ] == A.back());
        
        Mat B(A);
        CPPUNIT_ASSERT(B.rows() == 9);
        CPPUNIT_ASSERT(B.cols() == 9);
        CPPUNIT_ASSERT(B[ 0 ][ 0 ]  == 1.0);
        
        B[ 1 ][ 1 ] = 1.0;
        A = B;
        CPPUNIT_ASSERT(A == B);
        CPPUNIT_ASSERT(A[ 1 ][ 1 ]  == 1.0);
        
        A.clear(); 
        CPPUNIT_ASSERT(A.rows() == 0);
        CPPUNIT_ASSERT(A.cols() == 0);
        
        A.resize(3, 4);
        A *= 0.0;
        A += 1.0;
        CPPUNIT_ASSERT(sum(A[ 0 ]) == A.cols());
        
        Vec a(A.rows()); a.fill(x__ * 2.5);
        Vec b(A.cols()); b.fill((x__ + 1.0) * 4.2);
        
        B = A;
        Mat C(B);
        scaleMatrix(B, a, b);
        rank1Update(C, a, b);
        CPPUNIT_ASSERT((C-B) == A);
        CPPUNIT_ASSERT(::fabs(sum(A * b) - sum(b) * A.rows()) < TOLERANCE);
        
        A *= 1.11;
        A.round(0.1);
        CPPUNIT_ASSERT(::fabs(A[0][0] -1.1) < TOLERANCE);
        
    }
    
    void testFind(){
        typedef Vector < double > Vec;
        Vec x(10); x.fill(x__ + 1.0);
        Vec y(10); y.fill(x__ + 2.0);
        
        x[5] = 10.0;
        
        BVector b(y > x);
        CPPUNIT_ASSERT((y > x).size() == 10);
        CPPUNIT_ASSERT(find(y > x).size() == 9);
        CPPUNIT_ASSERT(find(x > y).size() == 1);
        CPPUNIT_ASSERT(find(x > y)[0] == 5);
        CPPUNIT_ASSERT(x(find(x > y))[0] == 10.0);
        
        x[ 2 ] = 5.0;
        x[ 3 ] = 5.0;
	// x = [1, 2, 5, 5, 5, 6, 7, 8, 9, 10]
        CPPUNIT_ASSERT((x == 5).size() == 10);
        CPPUNIT_ASSERT(find(x == 5).size() == 3);
        CPPUNIT_ASSERT(find(~(x == 5)).size() == 7);
        CPPUNIT_ASSERT(find(x == 5)[ 0 ] == 2);
        CPPUNIT_ASSERT(find(x <= 5).size() == 5);
        CPPUNIT_ASSERT(find(x > 5).size() == 5);
        CPPUNIT_ASSERT(find(x < 5).size() == 2);
        CPPUNIT_ASSERT(find((x > 5) & (x < 5)).size() == 0);
        CPPUNIT_ASSERT(find((x > 5) | (x < 5)).size() == 7);
        CPPUNIT_ASSERT(find(x+x == 4.).size() == 1);
        CPPUNIT_ASSERT(find(x+x >= 4.).size() == 9);
        CPPUNIT_ASSERT(find(x+x > 4.).size() == 8);
        CPPUNIT_ASSERT(find(abs(x+x) > 4.).size() == 8);
        CPPUNIT_ASSERT(find(~(abs(x+x) > 4.)).size() == 2);
        x[0] = ::log(0);
        CPPUNIT_ASSERT(find(isInf(x)).size() == 1);
        CPPUNIT_ASSERT(find(isInf(x+x)).size() == 1);
        x[1] += ::sqrt(-1) ;
        CPPUNIT_ASSERT(find(isNaN(x)).size() == 1);
        CPPUNIT_ASSERT(find(isInfNaN(x)).size() == 2);
        //CPPUNIT_ASSERT(find(isnan(x+x)).size() == 0);
    }
    
    void testMatrix(){
        testMatrix_< double >();
//        testMatrix_< float >();
    }
    
    void testIO(){
        RVector v(100);
        randn(v);
        save(v, "test.vec");
        RVector t1("test.vec");
        CPPUNIT_ASSERT(v == t1);
        
        save(t1, "test", Binary);
        RVector t2("test", Binary);
        CPPUNIT_ASSERT(v == t2);
    }
    void setUp(){    
        v1_ = new RVector(10);
        v2_ = new RVector(*v1_);
    }

    void tearDown(){
        delete v1_;
        delete v2_;
    }

private:
    RVector * v1_, * v2_;
    
};

CPPUNIT_TEST_SUITE_REGISTRATION(VectorTest);