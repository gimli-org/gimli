#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

this_module_dir_path = os.path.abspath(os.path.dirname(sys.modules[__name__].__file__))
project_root = os.path.abspath(os.path.join(this_module_dir_path, "..", "..", ".."))

complete_path = lambda *args: os.path.join(project_root, *args)

import distutils.sysconfig

# Get the major Python version without patchlevel
# print("Python version:", distutils.sysconfig.get_python_version())
# print("Python lib:", distutils.sysconfig.get_python_lib(standard_lib=True))


class settings(object):
    module_name = "_pygimli_"
    gimli_path = "../src"
    gimli_defines = ""
    includesPaths = []

    @staticmethod
    def setup_environment():
        settings.pygccxml_path = "../../../pygccxml"
        settings.pyplusplus_path = "../../../pyplusplus"

        if os.path.exists("../../../gccxml-bin/bin"):
            settings.gccxml_path = os.path.abspath("../../../gccxml-bin/bin")
        elif os.path.exists("../../gccxml-bin/bin"):
            settings.gccxml_path = os.path.abspath("../../gccxml-bin/bin")

        sys.path.append(settings.pygccxml_path)
        sys.path.append(settings.pyplusplus_path)


if sys.platform == "linux2" or sys.platform == "linux":
    pass

elif sys.platform == "darwin":
    pass

elif sys.platform == "win32":
    settings.gimli_defines = "MINGW"
    settings.python_libs_path = distutils.sysconfig.get_python_lib(standard_lib=True)

    if os.path.exists("../../../boost/include"):
        settings.includesPaths.append(os.path.abspath("../../../boost/include"))

else:
    raise RuntimeError('There is no configuration for "%s" platform.' % sys.platform)

settings.setup_environment()
