
#include "boost/python/object.hpp"  //len function
#include "boost/python/ssize_t.hpp" //ssize_t type definition
#include "boost/python/detail/none.hpp"
#include <boost/mpl/int.hpp>
#include <boost/mpl/next.hpp>
#include "tuples.hpp"

#include "gimli.h"
#include "pos.h"
#include "vector.h"

namespace bpl = boost::python;


// ** TODO check if we need a delete somewhere for all the new stuff

namespace r_values_impl{

struct PyTuple2RVector3{

    typedef boost::tuples::tuple< double, double > xy_type;
    typedef boost::tuples::tuple< double, double, double > xyz_type;
    typedef bpl::from_py_sequence< xy_type > xy_converter_type;
    typedef bpl::from_py_sequence< xyz_type > xyz_converter_type;
    
    typedef GIMLI::RVector3 xyz_t;

    static void * convertible(PyObject * obj){
        if (xy_converter_type::convertible(obj) || 
             xyz_converter_type::convertible(obj)){
            return obj;
        } else{
            return NULL;
        }
    }

    static void construct(PyObject* obj, bpl::converter::rvalue_from_python_stage1_data * data){

        typedef bpl::converter::rvalue_from_python_storage< xyz_t > xyz_storage_t;
        
        xyz_storage_t * the_storage = reinterpret_cast< xyz_storage_t * >(data);
        void * memory_chunk = the_storage->storage.bytes;

        double x(0.0), y(0.0), z(0.0);
        
        bpl::tuple py_tuple(bpl::handle<>(bpl::borrowed(obj)));
        
        if (3 == bpl::len(py_tuple)){
            boost::tuples::tie(x, y, z) = xyz_converter_type::to_c_tuple(obj);
        } else if (2 == bpl::len(py_tuple)){
            boost::tuples::tie(x, y) = xy_converter_type::to_c_tuple(obj);
        }        
                
        xyz_t * vec = new (memory_chunk) xyz_t(x, y, z);
        data->convertible = memory_chunk;
    }
};

struct PySequence2RVector{

    /*! Check if the object is convertible */
    static void * convertible(PyObject * obj){
                
        // is obj is a sequence
        if(!PySequence_Check(obj)){
            return NULL;
        }

        // has the obj a len method
        if(!PyObject_HasAttrString(obj, "__len__")){
            return NULL;
        }

        bpl::object py_sequence(bpl::handle<>(bpl::borrowed(obj)));
//         std::cout << "here am i 1 " << len(py_sequence) << std::endl;
        
        if (len(py_sequence) > 0) {
            
            bpl::object element = py_sequence[0];
            bpl::extract< double > type_checker(element);
            
            if(type_checker.check()){
                return obj;
            } else {
                std::cout << WHERE_AM_I << "element cannot converted to double" << std::endl;
            }
            
        } else {
            std::cout << WHERE_AM_I << " " << std::endl;
            return NULL;
        }
        // check if there is a valid converter
//         if(convertible_impl(py_sequence, boost::mpl::int_< 0 >(), length_type())){
//             return obj;
//         } else{
        return NULL;
        
    }

    /*! Convert obj into RVector */
    static void construct(PyObject* obj, bpl::converter::rvalue_from_python_stage1_data * data){
        
        bpl::object py_sequence(bpl::handle<>(bpl::borrowed(obj)));

        typedef bpl::converter::rvalue_from_python_storage< GIMLI::Vector< double > > storage_t;
         
        storage_t* the_storage = reinterpret_cast<storage_t*>(data);
        void* memory_chunk = the_storage->storage.bytes;
 
        GIMLI::Vector< double > * vec = new (memory_chunk) GIMLI::Vector< double >(len(py_sequence));
        data->convertible = memory_chunk;

        for (GIMLI::Index i = 0; i < vec->size(); i ++){
            vec->setVal(bpl::extract< double >(py_sequence[i]), i);
        }
    }
private:    
};

struct PySequence2StdVectorUL{

    /*! Check if the object is convertible */
    static void * convertible(PyObject * obj){
                
        // is obj is a sequence
        if(!PySequence_Check(obj)){
            return NULL;
        }

        // has the obj a len method
        if(!PyObject_HasAttrString(obj, "__len__")){
            return NULL;
        }

        bpl::object py_sequence(bpl::handle<>(bpl::borrowed(obj)));
        
        if (len(py_sequence) > 0) {
            
            bpl::object element = py_sequence[0];
            bpl::extract< GIMLI::Index > type_checker(element);
            
            if(type_checker.check()){
                return obj;
            } else {
                std::cout << WHERE_AM_I << "element cannot converted to GIMLI::Index: " << std::endl;
            }
            
        } else {
            std::cout << WHERE_AM_I << " " << std::endl;
            return NULL;
        }
        // check if there is a valid converter
//         if(convertible_impl(py_sequence, boost::mpl::int_< 0 >(), length_type())){
//             return obj;
//         } else{
        return NULL;
        
    }

    /*! Convert obj into RVector */
    static void construct(PyObject* obj, bpl::converter::rvalue_from_python_stage1_data * data){
        
        bpl::object py_sequence(bpl::handle<>(bpl::borrowed(obj)));

        typedef bpl::converter::rvalue_from_python_storage< std::vector < GIMLI::Index > > storage_t;
         
        storage_t* the_storage = reinterpret_cast<storage_t*>(data);
        void* memory_chunk = the_storage->storage.bytes;
 
        std::vector < GIMLI::Index > * vec = new (memory_chunk) std::vector < GIMLI::Index >(len(py_sequence));
        data->convertible = memory_chunk;

        for (GIMLI::Index i = 0; i < vec->size(); i ++){
            (*vec)[i] = bpl::extract< GIMLI::Index >(py_sequence[i]);
        }
    }
private:    
};

struct PySequence2StdVectorRVector3{

    /*! Check if the object is convertible */
    static void * convertible(PyObject * obj){
                
        // is obj is a sequence
        if(!PySequence_Check(obj)){
            return NULL;
        }

        // has the obj a len method
        if(!PyObject_HasAttrString(obj, "__len__")){
            return NULL;
        }

        bpl::object py_sequence(bpl::handle<>(bpl::borrowed(obj)));
        
        if (len(py_sequence) > 0) {
            
            bpl::object element = py_sequence[0];
            bpl::extract< GIMLI::Pos< double > > type_checker(element);
            
            if (type_checker.check()){
                return obj;
            } else {
                std::cout << WHERE_AM_I << "element cannot converted to GIMLI::Pos< double > " << std::endl;
            }
            
        } else {
            std::cout << WHERE_AM_I << " " << std::endl;
            return NULL;
        }
        // check if there is a valid converter
//         if(convertible_impl(py_sequence, boost::mpl::int_< 0 >(), length_type())){
//             return obj;
//         } else{
        return NULL;
        
    }

    /*! Convert obj into RVector */
    static void construct(PyObject* obj, bpl::converter::rvalue_from_python_stage1_data * data){
        
        bpl::object py_sequence(bpl::handle<>(bpl::borrowed(obj)));

        typedef bpl::converter::rvalue_from_python_storage< std::vector < GIMLI::Pos< double > > > storage_t;
         
        storage_t* the_storage = reinterpret_cast<storage_t*>(data);
        void* memory_chunk = the_storage->storage.bytes;
 
        std::vector < GIMLI::Pos < double > > *vec = new (memory_chunk) std::vector < GIMLI::Pos < double > >(len(py_sequence));
        data->convertible = memory_chunk;

        for (GIMLI::Index i = 0; i < vec->size(); i ++){
            (*vec)[i] = bpl::extract< GIMLI::Pos < double > >(py_sequence[i]);
        }
    }
private:    
};

} //r_values_impl

void register_pytuple_to_rvector3_conversion(){
    bpl::converter::registry::push_back(& r_values_impl::PyTuple2RVector3::convertible, 
                                        & r_values_impl::PyTuple2RVector3::construct, 
                                        bpl::type_id< GIMLI::Pos< double > >());
}

void register_pysequence_to_rvector_conversion(){
    bpl::converter::registry::push_back(& r_values_impl::PySequence2RVector::convertible, 
                                        & r_values_impl::PySequence2RVector::construct, 
                                        bpl::type_id< GIMLI::Vector< double > >());
}

void register_pysequence_to_StdVectorUL_conversion(){
    bpl::converter::registry::push_back(& r_values_impl::PySequence2StdVectorUL::convertible, 
                                        & r_values_impl::PySequence2StdVectorUL::construct, 
                                        bpl::type_id< std::vector< GIMLI::Index > >());
}

void register_pysequence_to_StdVectorRVector3_conversion(){
    bpl::converter::registry::push_back(& r_values_impl::PySequence2StdVectorRVector3::convertible, 
                                        & r_values_impl::PySequence2StdVectorRVector3::construct, 
                                        bpl::type_id< std::vector< GIMLI::Pos< double > > >());
}
