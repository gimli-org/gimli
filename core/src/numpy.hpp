#ifndef PYGIMLI__NUMPY__HPP
#define PYGIMLI__NUMPY__HPP

#include <numpy/arrayobject.h>

#include "gimli.h"


static int __numpy_initialized = 0;

// this fails for win since py39
// const static int __numpy_initialized = initNumpy();

inline int initNumpy(){
    // Needed or py* checks will segfault
    if (__numpy_initialized == 0){
        import_array2("Cannot import numpy.core.multiarray c-api for value converters.", 0);
        __numpy_initialized = 1;
    }
    return 0;
}

namespace GIMLI{


template < class ValueType > int npType(){
    __MS(sizeof(ValueType))
    __MS(typeid(ValueType).name())
    
    THROW_TO_IMPL
    return 0;
}

template < > inline int npType<bool>(){
    return NPY_BOOL;
}
template < > inline int npType<double>(){
    return NPY_DOUBLE;
}
template < > inline int npType<GIMLI::Complex>(){
    return NPY_COMPLEX128;
}
template < > inline int npType<GIMLI::Index>(){
    return NPY_UINT64;
}
template < > inline int npType<int>(){
    // int might by ambigous, check here .. MAC vs. WIN cs. LIN
    return NPY_INT32;
}

template < class ValueType, template < class > class Vec >
PyObject * toNumpyArray(Vec< ValueType > & vec){
    initNumpy();

    npy_intp length = (ssize_t)vec.size();
    
    //** this fails with typename==int!
    // PyObject * ret = PyArray_SimpleNewFromData(1, &length,
    //                                            npType<ValueType>(), 
    //                                            (void *)(&vec[0]));
    //** this fails with typename==int!

    PyObject * ret = PyArray_SimpleNew(1, &length, npType<ValueType>());
    std::memcpy(PyArray_DATA(reinterpret_cast<PyArrayObject*>(ret)),
                (void *)(&vec[0]), length * sizeof(ValueType));
    return ret;
}

// special case for std::vector < ValueType, ....> 
template < class ValueType, class Iter, template < class, class > class Vec >
PyObject * toNumpyArray(Vec< ValueType, Iter > & vec){
    initNumpy();
    
    npy_intp length = (ssize_t)vec.size();
 
    PyObject * ret = PyArray_SimpleNew(1, &length, npType<ValueType>());
    std::memcpy(PyArray_DATA(reinterpret_cast<PyArrayObject*>(ret)),
                (void *)(&vec[0]), length * sizeof(ValueType));
    return ret;
}



} // namespace GIMLI
#endif // PYGIMLI__NUMPY__HPP