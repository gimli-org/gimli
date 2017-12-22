# Changelog

## v1.0.5 (21/12/2017)

#### building and distribution

- [**building and distribution**] Make setup.py and conda's meta.yaml aware of pg.__version__ [#110](https://github.com/gimli-org/gimli/issues/110)

---

## v1.0.4 (17/11/2017)

#### bug

- [**bug**][**building and distribution**] make -DPYVERSION more flexible [#21](https://github.com/gimli-org/gimli/issues/21)

#### building and distribution

- [**building and distribution**] Binary Installers for Windows [#73](https://github.com/gimli-org/gimli/issues/73)
- [**building and distribution**] support for g++ 7.x.x [#64](https://github.com/gimli-org/gimli/issues/64)
- [**building and distribution**] Support python 3.6 [#59](https://github.com/gimli-org/gimli/issues/59)
- [**building and distribution**] os - depending python/boost_python setting in cmake [#28](https://github.com/gimli-org/gimli/issues/28)
- [**building and distribution**] Make binary installers available [#2](https://github.com/gimli-org/gimli/issues/2)

#### closed

- [**closed**] petro.resistivityArchie() returning empty vector [#109](https://github.com/gimli-org/gimli/issues/109)
- [**closed**] ERT.Simulate -- AttributeError 'NoneType' object has no attribute 'N' [#108](https://github.com/gimli-org/gimli/issues/108)
- [**closed**] petro.resistivityArchie interpolation returns 'corrupt mesh' error [#107](https://github.com/gimli-org/gimli/issues/107)
- [**closed**] Argument Error in pg.RMatrix in winpython [#106](https://github.com/gimli-org/gimli/issues/106)
- [**closed**] Assigning regions to large numbers of cells [#105](https://github.com/gimli-org/gimli/issues/105)
- [**closed**] List of keywords arguments (**kwargs) for pg.show [#104](https://github.com/gimli-org/gimli/issues/104)
- [**closed**] "gimli" CodeBlocks project + "bert" [#103](https://github.com/gimli-org/gimli/issues/103)
- [**closed**] A simple DC resistivity modelling case [#102](https://github.com/gimli-org/gimli/issues/102)
- [**closed**] Expanded mesh interpolation documentation? [#101](https://github.com/gimli-org/gimli/issues/101)
- [**closed**] pg.mplviewer.drawModel to view data [#100](https://github.com/gimli-org/gimli/issues/100)
- [**closed**] Using a pointcloud with pg.meshtools.createMesh  [#98](https://github.com/gimli-org/gimli/issues/98)
- [**closed**] Correct usage of the resample function from harmfit [#97](https://github.com/gimli-org/gimli/issues/97)
- [**closed**] Modifying grid.py (appendTriangleBoundary) to different y-max [#96](https://github.com/gimli-org/gimli/issues/96)
- [**closed**] Argument error in cg17/example-2_modelling.py [#95](https://github.com/gimli-org/gimli/issues/95)
- [**closed**] pg.createGrid fails with negative number range [#90](https://github.com/gimli-org/gimli/issues/90)
- [**closed**] Importing forward modeling results and performing ERT simulations [#84](https://github.com/gimli-org/gimli/issues/84)
- [**closed**] conda installation failing [#83](https://github.com/gimli-org/gimli/issues/83)
- [**closed**] ERT with forward modeling solution sourced from outside pygimli? [#82](https://github.com/gimli-org/gimli/issues/82)
- [**closed**] Install pygimli on ubuntu encounter the problem [#78](https://github.com/gimli-org/gimli/issues/78)
- [**closed**] How to install pyGIMLi for windows user ?  [#76](https://github.com/gimli-org/gimli/issues/76)
- [**closed**] drawPLC [#71](https://github.com/gimli-org/gimli/issues/71)
- [**closed**] Pygimli Windows 7 Installation [#70](https://github.com/gimli-org/gimli/issues/70)
- [**closed**] slownessWyllie import error [#69](https://github.com/gimli-org/gimli/issues/69)
- [**closed**] Two tests are failing with py 2.7 only [#68](https://github.com/gimli-org/gimli/issues/68)
- [**closed**] Compilation problems on Windows 7 with MSys64  [#51](https://github.com/gimli-org/gimli/issues/51)
- [**closed**] Base gallery build on version-controlled .matplotlibrc [#49](https://github.com/gimli-org/gimli/issues/49)

#### documentation

- [**documentation**] New structure for examples on website [#63](https://github.com/gimli-org/gimli/issues/63)

#### enhancement

- [**enhancement**] automatic rvalue conversion from array(dtype=int|long) to pg.RVector on demand [#99](https://github.com/gimli-org/gimli/issues/99)
- [**enhancement**][**feature request**] HDF format support [#65](https://github.com/gimli-org/gimli/issues/65)
- [**enhancement**][**feature request**] custom rvalue converter from/to RMatrix [#17](https://github.com/gimli-org/gimli/issues/17)
- [**enhancement**][**feature request**] Allow user-defined fillValue in pg.interpolate [#16](https://github.com/gimli-org/gimli/issues/16)

#### feature request

- [**feature request**] Easy python access to c++ inversion framework [#54](https://github.com/gimli-org/gimli/issues/54)
- [**feature request**] Importer: tetgen-1.4* and tetgen 1.5* [#7](https://github.com/gimli-org/gimli/issues/7)

---

## v0.91 (01/06/2017)

#### API change

- [**API change**] Introduction of 2 new submodules [#33](https://github.com/gimli-org/gimli/issues/33)
- [**API change**] I want to remove cell.attribute() [#27](https://github.com/gimli-org/gimli/issues/27)

#### bug

- [**bug**] Cannot show rotated PLCs crossing y=0 due to bug in Triangle(-Wrapper) [#47](https://github.com/gimli-org/gimli/issues/47)
- [**bug**][**building and distribution**] triangle build on MacOSX [#24](https://github.com/gimli-org/gimli/issues/24)

#### building and distribution

- [**building and distribution**] boost 1.60 compatibility issues [#31](https://github.com/gimli-org/gimli/issues/31)
- [**building and distribution**] make curl installers work on CI platforms [#26](https://github.com/gimli-org/gimli/issues/26)
- [**building and distribution**] Gimli install [#23](https://github.com/gimli-org/gimli/issues/23)
- [**building and distribution**] Remove old und unused binary apps [#10](https://github.com/gimli-org/gimli/issues/10)
- [**building and distribution**] Reduce dependencies [#3](https://github.com/gimli-org/gimli/issues/3)

#### closed

- [**closed**] missing argument in polytools [#56](https://github.com/gimli-org/gimli/issues/56)
- [**closed**] Ubuntu: cmake issue - "cannot find -lpthreads" [#53](https://github.com/gimli-org/gimli/issues/53)
- [**closed**] Use castxml binaries to avoid clang and llvm related issues [#52](https://github.com/gimli-org/gimli/issues/52)
- [**closed**] Ship self-testing binaries [#50](https://github.com/gimli-org/gimli/issues/50)
- [**closed**] Problem with dependencies to libpython3.4m.so, libpython3.4m.so.1.0 and libboost_python-py34.so while making pygimli under python3.5 [#45](https://github.com/gimli-org/gimli/issues/45)
- [**closed**] Issue with gimli installation [#44](https://github.com/gimli-org/gimli/issues/44)
- [**closed**] Strange aspect of pyGIMLi logo [#42](https://github.com/gimli-org/gimli/issues/42)
- [**closed**] Some examples are broken due to API changes [#41](https://github.com/gimli-org/gimli/issues/41)
- [**closed**] castXML cannot be build with llvm-3.8 [#40](https://github.com/gimli-org/gimli/issues/40)
- [**closed**] Runtime error  [#39](https://github.com/gimli-org/gimli/issues/39)
- [**closed**] Problem building pygimli with python3.5 [#38](https://github.com/gimli-org/gimli/issues/38)
- [**closed**] Make convertMesh gmsh compatible (also for anaconda users) [#37](https://github.com/gimli-org/gimli/issues/37)
- [**closed**] Problem with pygimli installation [#36](https://github.com/gimli-org/gimli/issues/36)
- [**closed**] Problem with castxml [#35](https://github.com/gimli-org/gimli/issues/35)
- [**closed**] Remove cppunit as a hard dependency [#34](https://github.com/gimli-org/gimli/issues/34)
- [**closed**] ImportError: No module named _pygimli_ [#32](https://github.com/gimli-org/gimli/issues/32)
- [**closed**] building castxml on ubuntu 14.04 [#30](https://github.com/gimli-org/gimli/issues/30)
- [**closed**] ImportError: dynamic module does not define init function (PyInit___pygimli_) [#29](https://github.com/gimli-org/gimli/issues/29)

#### documentation

- [**documentation**] Documentation issues [#12](https://github.com/gimli-org/gimli/issues/12)

#### enhancement

- [**enhancement**][**testing**] pg.test() writes files in current working directory [#43](https://github.com/gimli-org/gimli/issues/43)

---

## v0.9 (06/10/2015)

#### API change

- [**API change**][**bug**][**building and distribution**] make pygimli fails [#18](https://github.com/gimli-org/gimli/issues/18)

#### bug

- [**bug**][**documentation**] Documentation is not building with sphinx 1.3.1 [#11](https://github.com/gimli-org/gimli/issues/11)

#### building and distribution

- [**building and distribution**][**enhancement**] Avoid static link to libgimli.so in _pygimli_.so [#13](https://github.com/gimli-org/gimli/issues/13)
- [**building and distribution**] Force dependency update via cmake [#9](https://github.com/gimli-org/gimli/issues/9)
- [**building and distribution**] make: *** [pygimli] Error 2 [#8](https://github.com/gimli-org/gimli/issues/8)
- [**building and distribution**] Provide support for newer compilers via CastXML [#6](https://github.com/gimli-org/gimli/issues/6)

#### closed

- [**closed**] bin directory empty after make [#20](https://github.com/gimli-org/gimli/issues/20)

#### enhancement

- [**enhancement**][**feature request**] Array conversion methods for more/all vector/matrix formats [#5](https://github.com/gimli-org/gimli/issues/5)

#### feature request

- [**feature request**] Export triangle poly [#4](https://github.com/gimli-org/gimli/issues/4)
- [**feature request**] importer for plc (triangle) [#1](https://github.com/gimli-org/gimli/issues/1)

#### testing

- [**testing**] pg.test() should consider tests in python/tests [#19](https://github.com/gimli-org/gimli/issues/19)
