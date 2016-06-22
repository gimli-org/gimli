<!---
Readme for Github repository only. (Get's selected before *.rst file)
-->

![GIMLi](https://raw.githubusercontent.com/gimli-org/gimli/master/doc/_static/gimli.png)

[![Build Status](http://www.pygimli.org/build_status.svg)](http://www.pygimli.org/build.html)
[![Code Health](https://landscape.io/github/gimli-org/gimli/master/landscape.svg)](https://landscape.io/github/gimli-org/gimli/master)
[![Issue Stats](http://issuestats.com/github/gimli-org/gimli/badge/issue?style=flat)](http://issuestats.com/github/gimli-org/gimli)
[![Code Issues](https://www.quantifiedcode.com/api/v1/project/d0d835a5d75e4334a1c58389cafccaa0/badge.svg)](https://www.quantifiedcode.com/app/project/d0d835a5d75e4334a1c58389cafccaa0)


GIMLi is an open-source multi-method library for solving inverse
and forward modelling tasks.

##### What GIMLi is good for?:

- creating inversion applications (C++) and scripts (Python) for existing modules
- add your own forward calculations and build a robust inversion quickly
- combining different geophysical methods in various ways
- doing modelling of different PDEs

##### What GIMLi is **NOT** good for?:

- for people that expect a ready-made GUI for interpreting their data

##### Install via curl
```bash
curl -Ls install.pygimli.org | bash
```

##### For Anaconda users (currently Linux only)
```bash
# Installation
conda install -c gimli pygimli
# Update
conda update -c gimli -f pygimli
```

##### Usage
```python
import pygimli as pg
print(pg.__version__)
```

Check www.pygimli.org for additional information.
