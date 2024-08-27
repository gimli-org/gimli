#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import environment_for_pygimli_build


WRAPPER_DEFINITION_RVector3 = """
#include <numpy/arrayobject.h>

PyObject * RVector3_getArray(GIMLI::RVector3 & vec){
    import_array2("Cannot import numpy c-api from pygimli hand_make_wrapper2", NULL);
    npy_intp length = 3;
    PyObject * ret = PyArray_SimpleNew(1, &length, NPY_DOUBLE);
    std::memcpy(PyArray_DATA(reinterpret_cast<PyArrayObject*>(ret)),
                (void *)(&vec[0]), length * sizeof(double));

    // check if array is contiguous here
    // ** possible fixed due to memcpy here
    //PyArray_XINCREF(ret);
    //Py_INCREF(ret); // das scheint ignoriert zu werden darum muessen wir aussen noch kopieren
    //Py_DECREF(ret);
    return ret;
}

"""
WRAPPER_REGISTRATION_RVector3 = [
    """def("array",
       &RVector3_getArray,
       "PyGIMLI Helper Function: extract a numpy array object from a RVector3 ");""",
]

WRAPPER_DEFINITION_RVector = """
#include <numpy/arrayobject.h>

PyObject * RVector_getArray(GIMLI::RVector & vec){
    import_array2("Cannot import numpy c-api from pygimli hand_make_wrapper", NULL);
    npy_intp length = (ssize_t)vec.size();

    PyObject * ret = PyArray_SimpleNew(1, &length, NPY_DOUBLE);
    std::memcpy(PyArray_DATA(reinterpret_cast<PyArrayObject*>(ret)),
                (void *)(&vec[0]), length * sizeof(double));
        
    return ret;
}

"""
WRAPPER_REGISTRATION_RVector = [
    """def("array", &RVector_getArray,
       "PyGIMLI Helper Function: extract a numpy array object from a RVector ");""",
]

WRAPPER_DEFINITION_CVector = """
#include <numpy/arrayobject.h>

PyObject * CVector_getArray(GIMLI::CVector & vec){
    import_array2("Cannot import numpy c-api from pygimli hand_make_wrapper", NULL);
    npy_intp length = (ssize_t)vec.size();

    PyObject * ret = PyArray_SimpleNew(1, &length, NPY_COMPLEX128);
    std::memcpy(PyArray_DATA(reinterpret_cast<PyArrayObject*>(ret)),
                (void *)(&vec[0]), length * sizeof(GIMLI::Complex));
    
    return ret;
}

"""
WRAPPER_REGISTRATION_CVector = [
    """def("array", &CVector_getArray,
           "PyGIMLI Helper Function: extract a numpy array object from a CVector ")
    """,
]

WRAPPER_DEFINITION_BVector = """
#include <numpy/arrayobject.h>

PyObject * BVector_getArray(GIMLI::BVector & vec){
    import_array2("Cannot import numpy c-api from pygimli hand_make_wrapper", NULL);
    npy_intp length = (ssize_t)vec.size();

    PyObject * ret = PyArray_SimpleNew(1, &length, NPY_BOOL);
    std::memcpy(PyArray_DATA(reinterpret_cast<PyArrayObject*>(ret)),
                (void *)(&vec[0]), length * sizeof(bool));
        
    return ret;

    //PyObject * ret = PyArray_SimpleNew(1, &length, NPY_BOOL);
    //char * cout = (char *)PyArray_DATA(reinterpret_cast<PyArrayObject*>(ret));
    //for (ssize_t i=0; i<length; i++)  { cout[i] = vec[i]; }
    //return ret;
}

"""
WRAPPER_REGISTRATION_BVector = [
    """def("array", &BVector_getArray,
       "PyGIMLI Helper Function: extract a numpy array object from a BVector ");""",
]

WRAPPER_DEFINITION_IndexArray = """
#include <numpy/arrayobject.h>

PyObject * IndexArray_getArray(GIMLI::IndexArray & vec){
    import_array2("Cannot import numpy c-api from pygimli hand_make_wrapper", NULL);
    npy_intp length = (ssize_t)vec.size();

    PyObject * ret = PyArray_SimpleNew(1, &length, NPY_LONG);
    std::memcpy(PyArray_DATA(reinterpret_cast<PyArrayObject*>(ret)),
                (void *)(&vec[0]), length * sizeof(GIMLI::Index));
        
    return ret;
}

"""
WRAPPER_REGISTRATION_IndexArray = [
    """def("array", &IndexArray_getArray,
       "PyGIMLI Helper Function: extract a numpy array object from a IndexArray ");""",
]

WRAPPER_DEFINITION_R3Vector = """
#include <numpy/arrayobject.h>

PyObject * R3Vector_getArray(GIMLI::R3Vector & vec){
    import_array2("Cannot import numpy c-api from pygimli hand_make_wrapper2", NULL);
    npy_intp length = (ssize_t)vec.size();

#ifdef MS_WIN64
    //long long int dim2 [] = {length, 3};
    npy_intp dim2 [] = {length, 3};
#else
    npy_intp dim2 [] = {length, 3};
#endif

    PyObject * ret = PyArray_SimpleNew(2, dim2, NPY_DOUBLE);

    // check if array is contiguous here

    std::memcpy(PyArray_DATA(reinterpret_cast<PyArrayObject*>(ret)),
                (void *) &GIMLI::toArray(vec)[0],
                (length * 3) * sizeof(double));
    return ret;
}

"""
WRAPPER_REGISTRATION_R3Vector = [
    """def("array", &R3Vector_getArray,
       "PyGIMLI Helper Function: extract a numpy array object from a R3Vector ");""",
]


WRAPPER_DEFINITION_General = """
bool checkDataWrapper()
{
    std::cout << "checkDataWrapper() called\\n";
    return true;
}

GIMLI::RVector * General_createRVector(boost::python::list listin){
    std::cout << "HEREAM_I GIMLI::RVector * General_createRVector(boost::python::list listin){" << std::endl;
    GIMLI::RVector * ret = new GIMLI::RVector(0);
    return ret;
}

"""
WRAPPER_REGISTRATION_General = [
    """bp::def( "checkDataWrapper", &checkDataWrapper);""",
    """bp::def( "General_createRVector", &General_createRVector,
                "PyGIMLI Helper Function: create a GIMLI::RVector from given list. Check with custom_rvalue",
                bp::return_value_policy< bp::reference_existing_object, bp::default_call_policies >());""",
]

##################################################################


def iter_as_generator_vector(cls):
    # print("ITER:", cls.name)

    try:
        code = os.linesep.join(
            [
                "typedef %(cls)s iter_type;",
                "generators::generator_maker_vector< iter_type >::register_< %(call_policies)s >( %(exposer_name)s );",
            ]
        )
        cls.add_registration_code(
            code
            % {
                "cls": cls.decl_string,
                "call_policies": cls.mem_fun("nextVal").call_policies.create_type(),
                "exposer_name": cls.class_var_name,
            },
            works_on_instance=False,
        )
        cls.include_files.append("generators.h")
        # print("OK")
    except:
        raise
        # print("FAILED ")


##########################################################################
##########################################################################


def apply_reg(class_, code):
    for c in code:
        class_.add_registration_code(c)


def apply(mb):
    print("Register 'Vector<double>' handmade wrapper")
    rt = mb.class_("Vector<double>")
    rt.add_declaration_code(WRAPPER_DEFINITION_RVector)
    apply_reg(rt, WRAPPER_REGISTRATION_RVector)

    print("Register 'Vector<Complex>' handmade wrapper")
    rt = mb.class_("Vector< std::complex< double > >")
    rt.add_declaration_code(WRAPPER_DEFINITION_CVector)
    apply_reg(rt, WRAPPER_REGISTRATION_CVector)

    print("Register 'Vector<bool>' handmade wrapper")
    rt = mb.class_("Vector<bool>")
    rt.add_declaration_code(WRAPPER_DEFINITION_BVector)
    apply_reg(rt, WRAPPER_REGISTRATION_BVector)

    # rt = mb.class_('IndexArray')
    # rt.add_declaration_code(WRAPPER_DEFINITION_IndexArray)
    # apply_reg(rt, WRAPPER_REGISTRATION_IndexArray)

    # rt = mb.class_('Vector<GIMLI::Index>')
    # rt.add_declaration_code(WRAPPER_DEFINITION_IndexArray)
    # apply_reg(rt, WRAPPER_REGISTRATION_IndexArray)

    rt = mb.class_("Vector< GIMLI::Pos >")
    rt.add_declaration_code(WRAPPER_DEFINITION_R3Vector)
    apply_reg(rt, WRAPPER_REGISTRATION_R3Vector)

    try:
        rt = mb.class_("Pos")
        rt.add_declaration_code(WRAPPER_DEFINITION_RVector3)
        apply_reg(rt, WRAPPER_REGISTRATION_RVector3)
        # rt.add_registration_code ("""def(bp::init< PyObject * >((bp::arg("value"))))""")

        mb.add_declaration_code(WRAPPER_DEFINITION_General)
        apply_reg(mb, WRAPPER_REGISTRATION_General)
    except:
        pass

    # vec_iterators = mb.classes(lambda cls: cls.name.startswith('R'))

    vec_iterators = mb.classes(lambda cls: cls.name.startswith("VectorIterator"))
    for cls in vec_iterators:
        iter_as_generator_vector(cls)
