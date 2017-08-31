<!---
Readme for Github repository only. (Get's selected before *.rst file)
-->

<a href="https://www.pygimli.org">
  <img src="https://www.pygimli.org/_static/gimli_logo.svg" width="50%">
</a>

[![Build Status](https://www.pygimli.org/build_status.svg)](https://www.pygimli.org/build.html)
[![Code Health](https://landscape.io/github/gimli-org/gimli/master/landscape.svg)](https://landscape.io/github/gimli-org/gimli/master)
[![Issue Stats](http://issuestats.com/github/gimli-org/gimli/badge/issue?style=flat)](http://issuestats.com/github/gimli-org/gimli)
[![license](https://img.shields.io/github/license/gimli-org/gimli.svg?style=flat-square)](https://pygimli.org/license.html)
[![Anaconda-Server Badge](https://anaconda.org/gimli/pygimli/badges/license.svg)](https://anaconda.org/gimli/pygimli)

pyGIMLi is an open-source library for modelling and inversion and in geophysics. The object-oriented library provides management for structured and unstructured meshes in 2D and 3D, finite-element and finite-volume solvers, various geophysical forward operators, as well as Gauss-Newton based frameworks for constrained, joint and fully-coupled inversions with flexible regularization.

What is pyGIMLi suited for?

-   analyze, visualize and invert geophysical data in a reproducible manner
-   forward modelling of (geo)physical problems on complex 2D and 3D geometries
-   inversion with flexible controls on a-priori information and regularization
-   combination of different methods in constrained, joint and fully-coupled inversions
-   teaching applied geophysics (e.g. in combination with [Jupyter notebooks])

What is pyGIMLi **NOT** suited for?

-   for people that expect a ready-made GUI for interpreting their data

[jupyter notebooks]: http://jupyter-notebook.readthedocs.io/en/latest/notebook.html#notebook-documents

##### Binaries (Windows)

See binaries on <https://www.pygimli.org/installation.html#win>

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

##### Install via curl

```bash
curl -Ls install.pygimli.org | bash
```

##### Import convention

```python
import pygimli as pg
print(pg.__version__)
```

Check www.pygimli.org for additional information, detailed installation
instructions and many examples.

#### Citing pyGIMLi

If you use pyGIMLi for your work, please cite as:

> Rücker, C., Günther, T., Wagner, F.M., 2017. pyGIMLi: An open-source library for modelling and inversion in geophysics, Computers and Geosciences, 109, 106-123, doi: 10.1016/j.cageo.2017.07.011 ([Download PDF]).

[download pdf]: http://www.sciencedirect.com/science/article/pii/S0098300417300584/pdfft?md5=44253eaacd5490e3fb32210671672496&pid=1-s2.0-S0098300417300584-main.pdf

BibTeX code:

```sourceCode
@article{Ruecker2017,
  title = "{pyGIMLi}: An open-source library for modelling and inversion in geophysics",
  journal = "Computers and Geosciences",
  volume = "109",
  number = "",
  pages = "106--123",
  year = "2017",
  issn = "0098-3004",
  doi = "10.1016/j.cageo.2017.07.011",
  url = "http://www.sciencedirect.com/science/article/pii/S0098300417300584",
  author = "Carsten R\"ucker and Thomas G\"unther and Florian M. Wagner"
}
```

##### License

pyGIMLi is distributed under the terms of the **Apache 2.0** license. Details on
the license agreement can be found [here].

[here]: https://www.pygimli.org/license.html
