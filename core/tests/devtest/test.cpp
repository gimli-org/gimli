#include <iostream>

#include "gimli.h"
#include "sparsemapmatrix.h"
#include "vector.h"

using namespace GIMLI;


int main(int argc, char ** argv){

    RVector a(2);
    print(a);

    RVector b{1.0, 2.0};
    
    print(b.size());
    print(b);


    RSparseMapMatrix A(3,3);
    for (Index i = 0; i < A.rows(); i ++ ){
        for (Index j = 0; j < A.cols(); j ++ ){
            A[i][j] = i*A.rows()+j;
        }
    }

    __MS(A.values())
    __MS(A.rowIDs())
    __MS(A.colIDs())
    A.reduce({0,1}, true);
    __MS(A.values())
    __MS(A.rowIDs())
    __MS(A.colIDs())
    //print(A)
    return 0;
}
