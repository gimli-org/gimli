# Changelog

## v1.0.6 (03/04/2018)
#### new functionality
- `pg.show(geom, markers=True)` to inspect region and boundary markers of PLCs (also works for meshes)
- `mt.interpolateAlongCurve`, add simple extrapolation for to little curve bounds.
- `pg.warn(*args), pg.info(*args), pg.debug(*args)` etc. easier logging interface #14
- `Pos=Vector3` alias
- `mt.createGrid`  proxy function
- static `ndim` property to `pg.Vector`, `pg.Matrix` etc.
- convenience calls: `pg.x`, `pg.y`, and `pg.z` to return appropriate arrays for (data, mesh, R3Vector, ndarray, ..)
- core: improve mesh entities selection via (BVector), e.g., `mesh.cells(mesh.cellMarkers() == 1)` 
is the same as `mesh.cells(pg.find(mesh.cellMarkers() == 1))`
- mplviewer new flag: showBoundary
- `pg.getConfigPath()` returns path to user config
- add more flexible way to define callables for boundary conditions

#### other changes
- `pg.test()` compatible with numpy 1.14 (and 1.13)
- documentation compatible with Sphinx > 1.6 (new `-i` flag in sphinx-autogen)
- increased visibility of hyperlinks in documentation
- removed io/ path #124

#### backward compatibility

- [**backward compatibility**] Backward compatibility of pg.interpolate [#131](https://github.com/gimli-org/gimli/issues/131)

#### bug

- [**bug**] Plotting issue for multiple subplots in Jupyter Notebook [#125](https://github.com/gimli-org/gimli/issues/125)
- fix implementation of Neumann boundary conditions.

#### enhancement

- [**enhancement**] pg.meshtools.cellDataToNodeData does not take scalar fields [#120](https://github.com/gimli-org/gimli/issues/120)
- CellBrowser: Highlight style, click now toggle, data and mesh can be changed

#### closed issues

- [**closed**] Drawing a colorbar with pg.mplviewer.drawModel [#130](https://github.com/gimli-org/gimli/issues/130)
- [**closed**] 'MemoryError' for pb.invert [#129](https://github.com/gimli-org/gimli/issues/129)
- [**closed**] FDEM dataset [#122](https://github.com/gimli-org/gimli/issues/122)
- [**closed**] petro.resistivityArchie interpolation returns 'corrupt mesh' error [#107](https://github.com/gimli-org/gimli/issues/107)
- [**closed**] API doc not building with sphinx 1.6.3 [#74](https://github.com/gimli-org/gimli/issues/74)
---

## v1.0.5 (22/12/2017)
Install with Ana- or Miniconda:

``` bash
conda install pygimli=v1.0.5
```
#### Fixes:
- automatic RValue conversion for scalar numpy data types to pg.core long 
- Crank-Nicolson solver is now working again with default theta setting

#### building and distribution

- [**building and distribution**] Make setup.py and conda's meta.yaml aware of pg.__version__ [#110](https://github.com/gimli-org/gimli/issues/110)

---

## v1.0.4 - Start of consequent versioning (22/12/2017)
Install with Ana- or Miniconda:

``` bash
conda install pygimli=v1.0.4
```