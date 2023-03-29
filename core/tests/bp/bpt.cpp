#include "boost/python.hpp"
#include <iostream>

void t(){
   std::cout << "test throw std::out_of_range" << std::endl;
   throw std::out_of_range ("std::out_of_range");
}

class A{
public: 
    A(int dim=1){
        std::cout << "init A with: " << dim << std::endl;
    }
};

namespace bp = boost::python;

BOOST_PYTHON_MODULE(_bpt_){

    bp::def("t", t, "");

    typedef bp::class_< A > A_exposer_t;
    A_exposer_t A_exposer = A_exposer_t( "A", "class doctest", 
        bp::init< bp::optional< int > >((bp::arg("dim")=(int)(1)), "init doctest"));
}
