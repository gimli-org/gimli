# -*- coding: utf-8 -*-
"""
Setup script for pygimli and python apps.

Usage:

    python setup.py install --user # for local installation
    sudo python setup.py install # for system-wide installation
    python setup.py develop # for developers (creates symlinks to source)

"""

import os
from setuptools import setup, find_packages
#import versioneer

apps = ["apps/" + app for app in os.listdir('apps') if not "template" in app]

try:
    with open(os.path.join("README.rst")) as f:
        long_description = f.read()
except:
    long_description = "Geophysical Inversion and Modelling Library Core"

setup(name="pygimli",
    ############################################################################
    # Do not edit next two lines, do "git tag v1.x.x; git push --tags" instead.
    version="1.2.0",
    # cmdclass=versioneer.get_cmdclass(),
    ############################################################################
    description="Geophysical Inversion and Modelling Library Core",
    long_description=long_description,
    author="Carsten Rücker, Thomas Günther, Florian Wagner",
    author_email="mail@pygimli.org",
    license="Apache 2.0",
    url="https://www.pygimli.org",
    packages=find_packages(),
    package_data={'': ['*.so', '*.dll', '*.pyd']},
    scripts=apps
    )
