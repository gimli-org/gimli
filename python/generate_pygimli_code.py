#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import sys
import string

from optparse import OptionParser

from environment_for_pygimli_build import settings

optionParser = OptionParser("usage: %prog [options]")
optionParser.add_option("", "--extra-includes", dest="extraIncludes")
optionParser.add_option("", "--extra-path", dest="extraPath")
optionParser.add_option("", "--gccxml", dest="gccxml")

(options, args) = optionParser.parse_args()

if options.gccxml:
    settings.gccxml_path = options.gccxml

if options.extraPath:
    sys.path.append(options.extraPath)


import hand_made_wrappers

from pygccxml import parser
from pygccxml import declarations
from pygccxml.declarations import access_type_matcher_t
from pyplusplus import code_creators, module_builder, messages, decl_wrappers
from pyplusplus.module_builder import call_policies
from pyplusplus.decl_wrappers.doc_extractor import doc_extractor_i

import hashlib

MAIN_NAMESPACE = 'GIMLI'


def samefile(sourcefile, destfile):
    """
    """
    if not os.path.exists(destfile):
        return False
    if not os.path.exists(sourcefile):
        return False

    inhash = hashlib.md5()
    inhash.update(open(sourcefile).read().encode('utf-8'))

    outhash = hashlib.md5()
    outhash.update(open(destfile).read().encode('utf-8'))

    if inhash.digest() != outhash.digest():
        return False

    # probably don't need these as hash should handle it fine..
    if os.stat(sourcefile).st_mtime > os.stat(destfile).st_mtime:
        return False
    if os.stat(sourcefile).st_size != os.stat(destfile).st_size:
        return False

    return True


class decl_starts_with (object):

    def __init__(self, prefix):
        self.prefix = prefix

    def __call__(self, decl):
        return self.prefix in decl.name


def exclude(method, return_type='', name='', symbol=''):
    for funct in return_type:
        if len(funct):
            fun = method(return_type=funct, allow_empty=True)

            for f in fun:
                print("exclude return type", f)
                f.exclude()

    for funct in name:
        if len(funct):
            fun = method(name=funct, allow_empty=True)

            for f in fun:
                print("exclude name", f)
                f.exclude()

    for funct in symbol:
        if len(funct):
            fun = method(symbol=funct, allow_empty=True)

            for f in fun:
                print("exclude symbol", f)
                f.exclude()


def setMemberFunctionCallPolicieByReturn(mb, MemberRetRef, callPolicie):
    for ref in MemberRetRef:
        memFuns = mb.global_ns.member_functions(
            return_type=ref,
            allow_empty=True)
        print(ref, len(memFuns))

        for memFun in memFuns:
            if memFun.call_policies:
                continue
            else:
                memFun.call_policies = \
                    call_policies.return_value_policy(callPolicie)


class docExtractor(doc_extractor_i):

    def __init__(self):
        doc_extractor_i.__init__(self)
        pass

    # def __call__(self, decl):
        # print "__call__(self, decl):", decl
        # print decl.location.file_name
        # print decl.location.line
        # return ""Doku here""

    def escape_doc(self, doc):
        return '"'+doc+'"'

    def extract(self, decl):
        print(decl.location.file_name)
        print(decl.location.line)
        print("extract(self, decl):", decl)
        return "Doku coming soon"


def generate(defined_symbols, extraIncludes):
#    messages.disable(
# messages.W1005 # using a non public variable type for argucments or returns
# Warnings 1020 - 1031 are all about why Py++ generates wrapper for class X
# , messages.W1009 # check this
# , messages.W1014 # check this
#        , messages.W1020
#        , messages.W1021
#        , messages.W1022
#        , messages.W1023
#        , messages.W1024
#        , messages.W1025
#        , messages.W1026
#        , messages.W1027
#        , messages.W1028
#        , messages.W1029
#        , messages.W1030
#        , messages.W1031
# , messages.W1036 # check this
#)

    print("Install SRC:  ", os.path.abspath(__file__))
    print("Execute from: ", os.getcwd())

    sourcedir = os.path.dirname(os.path.abspath(__file__))
    sourceHeader = os.path.abspath(sourcedir + "/" + r"pygimli.h")
    gimliInclude = os.path.dirname(
        os.path.abspath(
            sourcedir +
            "/../src/" +
            r"gimli.h"))
    settings.includesPaths.append(gimliInclude)

    xml_cached_fc = parser.create_cached_source_fc(
        sourceHeader,
        settings.module_name +
        '.cache')
    #xml_cached_fc = parser.create_cached_source_fc(os.path.join(r"pygimli.h"), settings.module_name + '.cache')

    import platform

    defines = ['PYGIMLI_GCCXML', 'HAVE_BOOST_THREAD_HPP']

    if platform.architecture()[0] == '64bit' and platform.system() == 'Windows':

        if sys.platform == 'darwin':
            pass
        else:
            defines.append('_WIN64')
            print('Marking win64 for gccxml')

    for define in [settings.gimli_defines, defined_symbols]:
        if len(define) > 0:
            defines.append(define)

    try:
        if sys.platform == 'win32':
            # os.name == 'nt' (default on my mingw) results in wrong commandline
            # for gccxml
            os.name = 'mingw'
            gccxmlpath = settings.gccxml_path.replace('\\', '\\\\')
            gccxmlpath = settings.gccxml_path.replace('/', '\\')
            if not '.exe' in gccxmlpath:
                gccxmlpath += '\\gccxml.exe'

        else:
            gccxmlpath = settings.gccxml_path
    except Exception as e:
        print (str(e))
        raise Exception("Problems determine gccxml binary")

    settings.includesPaths.insert(0, os.path.abspath(extraIncludes))

    print("gccxml_path=%s" % gccxmlpath)
    print("working_directory=%s" % settings.gimli_path)
    print("include_paths=%s" % settings.includesPaths)
    print("define_symbols=%s" % defines)
    print("indexing_suite_version=2")

    mb = module_builder.module_builder_t([xml_cached_fc],
                                         gccxml_path=gccxmlpath,
                                         working_directory=settings.gimli_path,
                                         include_paths=settings.includesPaths,
                                         define_symbols=defines,
                                         indexing_suite_version=2)

    mb.classes().always_expose_using_scope = True
    mb.calldefs().create_with_signature = True

    hand_made_wrappers.apply(mb)

    global_ns = mb.global_ns
    global_ns.exclude()
    main_ns = global_ns.namespace(MAIN_NAMESPACE)
    main_ns.include()

    # START manual r-value converters
    rvalue_converters = [
        'register_pysequence_to_StdVectorUL_conversion',
        'register_pytuple_to_rvector3_conversion',
        'register_pysequence_to_rvector_conversion',
        'register_pysequence_to_StdVectorRVector3_conversion'
        ]

    for converter in rvalue_converters:
        mb.add_declaration_code('void %s();' % converter)
        mb.add_registration_code('%s();' % converter)

    # END manual r-value converters

    custom_rvalue_path = os.path.join(
        os.path.abspath(
            os.path.dirname(__file__)),
        'custom_rvalue.cpp')

    exclude(
        main_ns.variables,
        name=[
            'Triangle6_S1',
            'Triangle6_S2',
            'Triangle6_S3',
            'HexahedronFacesID',
            'Hexahedron20FacesID',
            'TetrahedronFacesID',
            'HexahedronSplit5TetID',
            'HexahedronSplit6TetID',
            'TriPrismFacesID',
            'TriPrimSplit3TetID',
            'NodeCoordinates',
            'EdgeCoordinates',
            'TriCoordinates',
            'QuadCoordinates',
            'TetCoordinates',
            'HexCoordinates',
            'PrismCoordinates',
            'PyramidCoordinates',
            'PyramidFacesID',
            'Tet10NodeSplit',
            'Tet10NodeSplitZienk',
            'Hex20NodeSplit',
            'Prism15NodeSplit',
            'Pyramid13NodeSplit'])

    for f in main_ns.declarations:
        if isinstance(f, decl_wrappers.calldef_wrapper.free_function_t):
            if (str(f.return_type).find('GIMLI::VectorExpr') != -1):
                f.exclude()

    exclude(
        main_ns.free_functions,
        return_type=[
            'float *',
            'float &',
            "::GIMLI::__VectorExpr< double, GIMLI::__VectorUnaryExprOp< double, GIMLI::VectorIterator< double >, GIMLI::ABS_ > >"],
        name=[
            'strReplaceBlankWithUnderscore',
            'toStr',
            'toInt',
            'toFloat',
            'toDouble',
            'str',
            'getRowSubstrings',
            'getNonEmptyRow',
            'getSubstrings',
            'abs',
            'type'])

    exclude(main_ns.free_operators,    name=[''],
            return_type=['::std::ostream &', '::std::istream &'])

    exclude(
        main_ns.classes, name=['ABS_', 'ACOT', 'ATAN', 'COS', 'COT', 'EXP',
                               'ABS_', 'LOG', 'LOG10', 'SIGN', 'SIN', 'SQRT', 'SQR', 'TAN', 'TANH',
                               'PLUS', 'MINUS', 'MULT', 'DIVID', 'BINASSIGN', 'cerrPtr',
                               'cerrPtrObject', 'coutPtr', 'coutPtrObject', 'deletePtr', 'edge_',
                               'distancePair_', 'IPCMessage', 'PythonGILSave', '__VectorExpr',
                               '__VectorUnaryExprOp', '__VectorBinaryExprOp', '__VectorValExprOp',
                               '__ValVectorExprOp',
                               'GIMLI::__ValVectorExprOp< double, GIMLI::VectorIterator< double >, GIMLI::MULT > >',
                               '::GIMLI::__VectorValExprOp< double, GIMLI::__VectorIterator< double >, GIMLI::MULT >',
                               '::GIMLI::__VectorBinaryExprOp< double, GIMLI::__VectorIterator< double >, GIMLI::__VectorIterator< double >, GIMLI::MULT>',
                               '::GIMLI::__VectorExpr< double, GIMLI::__VectorBinaryExprOp< double, GIMLI::__VectorIterator< double >, GIMLI::__VectorIterator< double >, GIMLI::MULT > >',
                               '::GIMLI::__VectorExpr< double, GIMLI::__VectorValExprOp< double, GIMLI::__VectorIterator< double >, GIMLI::MULT > >',
                               '::GIMLI::__VectorExpr< double, GIMLI::__VectorUnaryExprOp< double, GIMLI::__VectorIterator< double >, GIMLI::LOG10 > >',
                               '::GIMLI::__VectorExpr< double, GIMLI::__VectorUnaryExprOp< double, GIMLI::__VectorIterator< double >, GIMLI::LOG > >',
                               '::GIMLI::__VectorUnaryExprOp< double, GIMLI::__VectorIterator< double >, GIMLI::LOG10 >',
                               '::GIMLI::__VectorUnaryExprOp< double, GIMLI::__VectorIterator< double >, GIMLI::LOG >',
                               '::GIMLI::__VectorExpr<double, GIMLI::__VectorValExprOp<double, GIMLI::__VectorExpr<double, GIMLI::__ValVectorExprOp<double, GIMLI::__VectorIterator<double >, GIMLI::MULT > >, GIMLI::MULT > >',
                               '::GIMLI::__VectorValExprOp<double, GIMLI::__VectorExpr<double, GIMLI::__ValVectorExprOp<double, GIMLI::__VectorIterator<double >, GIMLI::MULT > >, GIMLI::MULT >.pypp.hpp',
                               'GIMLI::Expr<GIMLI::ExprIdentity>'])

    exclude(
        main_ns.member_functions,
        name=[
            'begin',
            'end',
            'val'],
        return_type=[''])

    exclude(main_ns.member_operators,  symbol=[''])

    mb.calldefs(access_type_matcher_t('protected')).exclude()
    mb.calldefs(access_type_matcher_t('private')).exclude()

    # setMemberFunctionCallPolicieByReturn(mb, [ '::GIMLI::Node &'
                                                #, '::GIMLI::Cell &'
                                                #, '::GIMLI::Boundary &'
                                                #, '::GIMLI::Shape &'
                                                #, '::GIMLI::Node *'
                                                #, '::GIMLI::Cell *'
                                                #, '::GIMLI::Boundary *'
                                                #, '::GIMLI::Shape *'
                                        #]
                                        #, call_policies.reference_existing_object)

    setMemberFunctionCallPolicieByReturn(
        mb,
        ['::std::string *', 'float *', 'double *', 'int *', 'long *',
         'long int *', 'long long int *', 'unsigned long long int *',
         '::GIMLI::Index *', '::GIMLI::SIndex *', 'bool *'],
        call_policies.return_pointee_value)

    setMemberFunctionCallPolicieByReturn(mb, ['::std::string &',
                                              'float &',
                                              'double &',
                                              'int &',
                                              'long &',
                                              'long int &',
                                              'long long int &',
                                              'unsigned long long int &',
                                              '::GIMLI::Index &',
                                              '::GIMLI::SIndex &',
                                              'bool &'
                                              ], call_policies.return_by_value)

    # setMemberFunctionCallPolicieByReturn(mb, ['::GIMLI::VectorIterator<double> &']
                                        #, call_policies.copy_const_reference)
    # setMemberFunctionCallPolicieByReturn(mb, [
                                                #,  'double &' ]
                                                #, call_policies.reference_existing_object)

    # call_policies.return_value_policy(call_policies.reference_existing_object)
    # call_policies.return_value_policy(call_policies.copy_non_const_reference)
    # call_policies.return_value_policy(call_policies.copy_const_reference)

     # addAutoConversions(mb)

   # excludeMemberByReturn(main_ns, ['::DCFEMLib::SparseMatrix<double> &'])
    #main_ns.classes(decl_starts_with(['STLMatrix']), allow_empty=True).exclude()
    #fun = mb.global_ns.member_functions('begin', allow_empty=True)
    # for f in fun:
        # f.exclude()

    # excludeFreeFunctionsByName(main_ns, ['strReplaceBlankWithUnderscore'
                               #'toStr', 'toInt', 'toFloat', 'toDouble',
                               #'getRowSubstrings', 'getNonEmptyRow', 'getSubstrings' ])

    #excludeFreeFunctionsByReturn(main_ns, [ 'float *', 'float &' ])
    #fun = ns.free_operators(return_type=funct, allow_empty=True)

    #excludeMemberOperators(main_ns, ['++', '--', '*'])

    # exclude all that does not match any predefined callpolicie

    excludeRest = True

    if excludeRest:
        mem_funs = mb.calldefs()

        for mem_fun in mem_funs:
            if mem_fun.call_policies:
                continue
            if not mem_fun.call_policies and \
                    (declarations.is_reference(mem_fun.return_type) or declarations.is_pointer(mem_fun.return_type)):
                # print mem_fun
                # mem_fun.exclude()
                mem_fun.call_policies = call_policies.return_value_policy(
                    call_policies.reference_existing_object)
                # mem_fun.call_policies = \
                #    call_policies.return_value_policy(call_policies.return_pointee_value)
                # mem_fun.call_policies = \
                #    call_policies.return_value_policy(call_policies.return_opaque_pointer)
                # mem_fun.call_policies = \
                 #   call_policies.return_value_policy(call_policies.copy_non_const_reference)

    # Now it is the time to give a name to our module
    from doxygen import doxygen_doc_extractor
    extractor = doxygen_doc_extractor()

    mb.build_code_creator(settings.module_name, doc_extractor=extractor)

    # It is common requirement in software world - each file should have license
    #mb.code_creator.license = '//Boost Software License(http://boost.org/more/license_info.html)'

    # I don't want absolute includes within code
    mb.code_creator.user_defined_directories.append(os.path.abspath('.'))

    # And finally we can write code to the disk
    def ignore(val):
        pass
    mb.split_module('./generated', on_unused_file_found=ignore)

    additional_files = [
        os.path.join(
            os.path.abspath(
                os.path.dirname(__file__)), 'custom_rvalue.cpp'), os.path.join(
            os.path.abspath(
                os.path.dirname(__file__)), 'generators.h'), os.path.join(
            os.path.abspath(
                os.path.dirname(__file__)), 'tuples.hpp')]

    for sourcefile in additional_files:
        p, filename = os.path.split(sourcefile)
        destfile = os.path.join('./generated', filename)

        if not samefile(sourcefile, destfile):
            shutil.copy(sourcefile, './generated')
            print("Updated ", filename, "as it was missing or out of date")


if __name__ == '__main__':

    defined_symbols = ''

    #from optparse import OptionParser

    #optionParser = OptionParser("usage: %prog [options]")
    #optionParser.add_option("", "--extra-includes", dest="extraIncludes")
    #optionParser.add_option("", "--extra-path", dest="extraPath")
    #optionParser.add_option("", "--gccxml", dest="gccxml")

    #(options, args) = optionParser.parse_args()
    # if options.gccxml:
        #settings.gccxml_path = options.gccxml

    # if options.extraPath:
        # sys.path.append(options.extraPath)

    generate(defined_symbols, options.extraIncludes)
