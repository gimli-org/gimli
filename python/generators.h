#ifndef PYGIMLI_GENERATORS__H
#define PYGIMLI_GENERATORS__H

#include "boost/python.hpp"
#include "vector.h"

namespace generators{
    
template< typename VectorIterator >

// don't use this until you know how to increase the reference counter for the underlying pythonobject for the vector
struct generator_maker_vector{
        
    typedef BOOST_DEDUCED_TYPENAME VectorIterator::value_type val_type;
    
    static void iter(const VectorIterator &){
        
    } //return_self call policies should be used
    
    static val_type next(VectorIterator & iter){
        if(!iter.hasMore()){
            boost::python::objects::stop_iteration_error();
            //will not come here
        }
        return iter.nextVal();
    }
    
    template< typename TNextCallPolicies, typename TPyClass>
    static void register_(TPyClass& py_cls ){
        typedef generator_maker_vector< VectorIterator > maker_type;
        
        py_cls.def("__iter__", &maker_type::iter, boost::python::return_self<>() );
        py_cls.def("__next__", &maker_type::next, TNextCallPolicies() );
        py_cls.def("next", &maker_type::next, TNextCallPolicies() );
    }
};
}

#endif