#include <iostream>
#include <Eigen/Dense>

#define USE_EIGEN3 1

#include "gimli.h"
#include "matrix.h"

using namespace GIMLI;

void toEigen(const RMatrix & A_, Eigen::MatrixXd & A){
    A.resize(A_.rows(), A_.cols());
    for (Index i=0; i < A.rows(); i ++){

        A(i, Eigen::all) = Eigen::Map <const Eigen::VectorXd>(&A_[i][0], A.cols());
    }    
}

void testMatMult(){
    Index m = 2;
    Index n = 3;
    Index k = 4;

    double *_A = new double[m * k];
    double *_B = new double[k * n];
    
    for (Index i = 0; i < m*k; i ++ ){_A[i] = i+1;}
    for (Index i = 0; i < k*n; i ++ ){_B[i] = i+1;}

    GIMLI::RMatrix A_(m, k, _A);
    GIMLI::RMatrix B_(k, n, _B);
    
    GIMLI::RMatrix C_;
    GIMLI::matMult_RM(A_, B_, C_, 1.0, 0.0);
    print(A_);
    print(B_);
    print(C_);
  
    SmallMatrix A;
    toEigenMatrix(A_, A);
        
    SmallMatrix B;
    toEigenMatrix(B_, B);
  
    std::cout << "mat:\n" << A << std::endl;
    std::cout << "mat:\n" << B << std::endl;

    std::cout << "mat:A*B\n" << A*B << std::endl;
  
    SmallMatrix C1;
    GIMLI::matMult(A, B, C1, 1.0, 0.0);
    std::cout << "C1\n" << C1 << std::endl;

    SmallMatrix C2;
    GIMLI::matMult(A, B.transpose(), C2, 1.0, 0.0);
    std::cout << "C2\n" << C2 << std::endl;

    
    GIMLI::RMatrix AT(k, m, _A);
    for (Index i = 0; i < k; i ++ ){
        for (Index j = 0; j < m; j ++ ){
            AT[i][j] = A_[j][i];
        }
    }
    GIMLI::matTransMult_RM(AT, B_, C_, 1.0, 0.0);
    print(C_);
    
    SmallMatrix C3;
    GIMLI::matTransMult(A.transpose(), B, C3, 1.0, 0.0);
    std::cout << "C3\n" << C3 << std::endl;

    SmallMatrix C4;
    GIMLI::matTransMult(A.transpose(), B.transpose(), C4, 1.0, 0.0);
    std::cout << "C4\n" << C4 << std::endl;
  
}

void testSetVal(){
    RMatrix A_(3, 3);
    RVector A0(2, 2);
    A_ += 1.0;

    A_(0).setVal(A0, 0, 2);
    print(A_);
    
    SmallMatrix A;
    toEigenMatrix(A_, A);
    A.array()+=1.0;
    
    SET_MAT_ROW_SLICE(A, 0, A0, 1.0, 0, 2);
    // print(Eigen::seqN(0, 2));

    // A(0, {0,1}) += Eigen::Map <const Eigen::VectorXd>(&A0[0], A0.size());
    // print(A);

    // A(0, Eigen::seqN(0, 2)) += Eigen::Map <const Eigen::VectorXd>(&A0[0], A0.size());
    // print(A);
    
    // A(0, Eigen::seq(0, 1)) += Eigen::Map <const Eigen::VectorXd>(&A0[0], A0.size());
    print(A);
    
    

}
    

int main(int argc, char ** argv){
    // testMatMult();
    testSetVal();
    return 0;
}
