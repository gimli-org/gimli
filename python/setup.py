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

with open(os.path.join("../README.rst"), encoding='utf-8') as f:
    long_description = f.read()

setup(name="pygimli",
      version="1.0rc",
      description="Geophysical Inversion and Modelling Library",
      long_description=long_description,
      author="Carsten Rücker, Thomas Günther, Florian Wagner",
      author_email="mail@pygimli.org",
      license="GPL",
      url="http://www.pygimli.org",
      packages=find_packages(),
      package_data={'': ['*.so']},
      scripts=apps)
