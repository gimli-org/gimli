#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os

this_module_dir_path = os.path.abspath ( os.path.dirname( sys.modules[__name__].__file__) )
project_root = os.path.abspath( os.path.join( this_module_dir_path, '..','..', '..' ) )
print project_root
complete_path = lambda *args: os.path.join( project_root, *args )

import distutils.sysconfig
print distutils.sysconfig.get_python_version() # Get the major Python version without patchlevel
print distutils.sysconfig.get_python_lib(standard_lib=True)

class settings:
    module_name = '_pygimli_'
    gimli_path = '../src'
    #gccxml_path                 = complete_path( 'gccxml_bin', 'v09', sys.platform, 'bin' )
    #pygccxml_path               = complete_path( 'pygccxml_dev' )
    #pyplusplus_path             = complete_path( 'pyplusplus_dev' )
    gimli_defines               = ''

    @staticmethod
    def setup_environment():
        sys.path.append( settings.pygccxml_path )
        sys.path.append( settings.pyplusplus_path )

if sys.platform == 'linux2':
    settings.gccxml_path        = '../../../gccxml-bin/bin/gccxml'
    settings.pygccxml_path      = '../../../pygccxml'
    settings.pyplusplus_path    = '../../../pyplusplus'
    
elif sys.platform == 'win32':
    settings.gimli_defines      = 'MINGW'
    settings.python_libs_path   = distutils.sysconfig.get_python_lib(standard_lib=True)
    settings.pygccxml_path      = '../../../pygccxml'
    settings.pyplusplus_path    = '../../../pyplusplus'
    
    if os.path.exists('../../../gccxml-bin/bin'):
        settings.gccxml_path = os.path.abspath('../../../gccxml-bin/bin')
    elif os.path.exists('../../gccxml-bin/bin'):
        settings.gccxml_path = os.path.abspath('../../gccxml-bin/bin')

else:
    raise RuntimeError( 'There is no configuration for "%s" platform.' % sys.platform )

setup_environment = settings.setup_environment()
