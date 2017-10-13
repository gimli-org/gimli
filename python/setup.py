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

apps = ["apps/" + app for app in os.listdir('apps') if not "template" in app]

try:
    # fails with py3 due to ascii encoding problem.
    with open(os.path.join("../README.rst")) as f:
        long_description = f.read()
except:
    long_description = "Geophysical Inversion and Modelling Library"

setup(name="pygimli",
      version="1.0.2",
      description="Geophysical Inversion and Modelling Library",
      long_description=long_description,
      author="Carsten Rücker, Thomas Günther, Florian Wagner",
      author_email="mail@pygimli.org",
      license="Apache",
      url="https://www.pygimli.org",
      packages=find_packages(),
      package_data={'': ['*.so', '*.dll', '*.pyd']},
      scripts=apps)
