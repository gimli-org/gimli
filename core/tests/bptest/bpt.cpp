#include "boost/python.hpp"

// #include "boost/python/object.hpp"  //len function

// #include <iostream>
// int initNumpy(){
//     std::cout << "try to import sys ... " << std::endl;
//     int st;
//     PyObject *numpy = PyImport_ImportModule("sys");
//     if (numpy == NULL) {
//         std::cout << "can't load sys" << std::endl;
//         return -1;
//     }
//     std::cout << "done!" << std::endl;
//     // import_array();
//     // import_array1(0);
//     // import_array2("Cannot import numpy.core._multiarray_umath c-api for rvalue converters.", 0);
//     return 0;
// }
// const static int __numpy_initialized = initNumpy();

double t1(){ return 1.0;}
int t2(){ return 2;} 

namespace bp = boost::python;
BOOST_PYTHON_MODULE(_bpt_){
    bp::def("t1", t1, ""); // works
    bp::def("t2", t2, ""); // will segfault
}
