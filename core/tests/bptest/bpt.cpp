#include "boost/python.hpp"

double t1(){ return 1.0;}
int t2(){ return 2;} 

namespace bp = boost::python;
BOOST_PYTHON_MODULE(_bpt_){
    bp::def("t1", t1, ""); // works
    bp::def("t2", t2, ""); // will segfault
}
