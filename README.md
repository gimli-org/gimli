<!---
Readme for Github repository only. (Get's selected before *.rst file)
-->

<a href="http://www.pygimli.org">
  <img src="http://www.pygimli.org/_static/gimli_logo.svg" width="50%">
</a>

[![Build Status](http://www.pygimli.org/build_status.svg)](http://www.pygimli.org/build.html)
[![Code Health](https://landscape.io/github/gimli-org/gimli/master/landscape.svg)](https://landscape.io/github/gimli-org/gimli/master)
[![Issue Stats](http://issuestats.com/github/gimli-org/gimli/badge/issue?style=flat)](http://issuestats.com/github/gimli-org/gimli)
[![license](https://img.shields.io/github/license/gimli-org/gimli.svg?style=flat-square)](https://pygimli.org/license.html)
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

#### Citing pyGIMLi

If you use pyGIMLi for your work, please cite as:

> Rücker, C., Günther, T., Wagner, F.M., pyGIMLi: An open-source library for modelling and inversion in geophysics, Computers and Geosciences (2017), doi: 10.1016/j.cageo.2017.07.011 ([Download PDF]).

BibTeX code:

```sourceCode
@article{Ruecker2017,
  title = "pyGIMLi: An open-source library for modelling and inversion in geophysics",
  journal = "Computers and Geosciences",
  volume = "",
  number = "",
  pages = "",
  year = "2017",
  note = "",
  issn = "0098-3004",
  doi = "10.1016/j.cageo.2017.07.011",
  author = "Carsten R\"ucker and Thomas G\"unther and Florian M. Wagner"
}
```

[download pdf]: https://www.pygimli.org/paper/Ruecker2017_CG_pyGIMLi.pdf

##### License

pyGIMLi is distributed under the terms of the **Apache 2.0** license. Details on
the license agreement can be found [here].

[here]: https://www.pygimli.org/license.html
