<!---
Readme for Github repository only. (Get's selected before *.rst file)
-->

<a href="http://www.pygimli.org">
  <img src="http://www.pygimli.org/_static/gimli_logo.svg" width="50%">
</a>

[![Build Status](http://www.pygimli.org/build_status.svg)](http://www.pygimli.org/build.html)
[![Code Health](https://landscape.io/github/gimli-org/gimli/master/landscape.svg)](https://landscape.io/github/gimli-org/gimli/master)
[![Issue Stats](http://issuestats.com/github/gimli-org/gimli/badge/issue?style=flat)](http://issuestats.com/github/gimli-org/gimli)
[![Anaconda-Server Badge](https://anaconda.org/gimli/pygimli/badges/license.svg)](https://anaconda.org/gimli/pygimli)

pyGIMLi is an open-source library for modelling and inversion and in geophysics. The object-oriented library provides management for structured and unstructured meshes in 2D and 3D, finite-element and finite-volume solvers, various geophysical forward operators, as well as Gauss-Newton based frameworks for constrained, joint and fully-coupled inversions with flexible regularization.

What pyGIMLi is suited for?

-   analyze, visualize and invert geophysical data in a reproducible manner
-   forward modelling of (geo)physical problems on complex 2D and 3D geometries
-   inversion with flexible controls on a-priori information and regularization
-   combination of different methods in constrained, joint and fully-coupled inversions
-   teaching applied geophysics (e.g. in combination with [Jupyter notebooks])

What pyGIMLi is **NOT** suited for?

-   for people that expect a ready-made GUI for interpreting their data

[jupyter notebooks]: http://jupyter-notebook.readthedocs.io/en/latest/notebook.html#notebook-documents

##### Install via curl

```bash
curl -Ls install.pygimli.org | bash
```

##### For Anaconda users (currently Linux only)

[![Anaconda-Server Badge](https://anaconda.org/gimli/pygimli/badges/installer/conda.svg)](https://conda.anaconda.org/gimli)
[![Anaconda-Server Badge](https://anaconda.org/gimli/pygimli/badges/downloads.svg)](https://anaconda.org/gimli/pygimli)

```bash
# Add gimli and conda-forge channels (only once)
conda config --add channels gimli --add channels conda-forge

# Install pygimli
conda install -f pygimli

# Test pygimli
python -c "import pygimli; pygimli.test()"

# Update pygimli
conda update -f pygimli
```

##### Import convention

```python
import pygimli as pg
print(pg.__version__)
```

Check www.pygimli.org for additional information, detailed installation
instructions and many examples.

##### License

pyGIMLi is distributed under the terms of the **Apache 2.0** license. Details on
the license agreement can be found [here].

[here]: https://www.pygimli.org/license.html
