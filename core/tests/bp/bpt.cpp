#include "boost/python.hpp"
#include <iostream>

void t(){
   std::cout << "test throw std::out_of_range" << std::endl;
   throw std::out_of_range ("std::out_of_range");
}

namespace bp = boost::python;
BOOST_PYTHON_MODULE(_bpt_){
    bp::def("t", t, "");
}
