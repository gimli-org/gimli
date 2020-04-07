import pygimli as pg

from .utils import pgMesh2pvMesh
pv = pg.optImport('pyvista', requiredFor="proper visualization in 3D")


def drawMesh(ax, mesh, notebook=False, **kwargs):
    """

    Parameters
    ----------
    ax: pyvista.Plotter [optional]
        The plotter to draw everything. If none is given, one will be created.
    mesh: pg.Mesh
        The mesh to show.
    notebook: bool [False]
        Sets the plotter up for jupyter notebook/lab.
    cmap: str ['viridis']
        The colormap string.

    Returns
    -------
    ax: pyvista.Plotter [optional]
        The plotter
    """
    # sort out a few kwargs to not confuse the plotter initialization
    show_edges = kwargs.pop('show_edges', True)
    opacity = kwargs.pop('alpha', 1)
    cmap = kwargs.pop('cmap', None)
    color = kwargs.pop('color', 'k')
    style = kwargs.pop('style', 'wireframe')
    returnActor = kwargs.pop('returnActor', False)

    if ax is None:
        ax = pv.Plotter(notebook=notebook, **kwargs)

    ax.show_bounds(all_edges=True, minor_ticks=True)
    ax.add_axes()

    if isinstance(mesh, pg.Mesh):
        mesh = pgMesh2pvMesh(mesh)

    _actor = ax.add_mesh(
        mesh,  # type: pv.UnstructuredGrid
        cmap=cmap,
        color=color,
        style=style,
        show_edges=show_edges,
        opacity=opacity,
        )

    if returnActor:
        return ax, _actor
    else:
        return ax


def drawModel(ax=None, mesh=None, data=None, **kwargs):
    """
    Draw the mesh with given data.

    Parameters
    ----------
    ax: pv.Plotter() [None]
        Pyvista's basic Plotter to add the mesh to.
    mesh: pg.Mesh
        The Mesh to plot.
    data: iterable
        Data that should be displayed with the mesh.
    """
    if all(v is None for v in [ax, mesh, data]):
        pg.critical("At least mesh or data should not be None")
        return None

    mesh = pgMesh2pvMesh(mesh, data, kwargs.pop('label', 'data'))
    if data is not None:
        kwargs['style'] = 'surface'
        kwargs['color'] = None

    if 'cmap' not in kwargs:
        kwargs['cmap'] = 'viridis'
    return drawMesh(ax, mesh, **kwargs)


def drawSensors(ax, sensors, diam=0.01, color='grey', **kwargs):
    """
    Draw the sensor positions to given mesh or the the one in given plotter.

    Parameters
    ----------
    ax: pyvista.Plotter
        The plotter to draw everything. If none is given, one will be created.
    sensors: iterable
        Array-like object containing tuple-like (x, y, z) positions.
    diam: float [0.01]
        Radius of sphere markers.
    color: str ['grey']
        Color of sphere markers.

    Returns
    -------
    ax: pyvista.Plotter
        The plotter containing the mesh and drawn electrode positions.
    """
    for pos in sensors:
        s = pv.Sphere(radius=diam/2, center=pos)
        ax.add_mesh(s, color=color, **kwargs)

    return ax


def drawSlice(ax, mesh, normal=[1, 0, 0], **kwargs):
    """

    Note
    ----
    Possible kwargs are:
    normal: tuple(float), str
    origin: tuple(float)
    generate_triangles: bool, optional
    contour: bool, optional

    They can be found at https://docs.pyvista.org/core/filters.html?highlight=slice_orthogonal#1pyvista.CompositeFilters.slice
    """
    label = kwargs.pop('label', 'data')
    data = kwargs.pop('data', None)
    mesh = pgMesh2pvMesh(mesh, data, label)

    try:
        single_slice = mesh.slice(normal, **kwargs)

    except AssertionError as e:
        # 'contour' kwarg only works with point data and breaks execution
        pg.error(e)
    else:
        # REVIEW: bounds and axes might be confused with the outline..?!
        outline = mesh.outline()
        ax.add_mesh(outline, color="k")
        ax.add_mesh(single_slice)

    return ax


def drawStreamLines(ax, mesh, data, label=None, radius=0.01, **kwargs):
    """
    Draw streamlines of given data.

    PyVista streamline needs a vector field of gradient data per cell.

    Parameters
    ----------
    ax: pv.Plotter() [None]
        The plotter that should be used for visualization.
    mesh: pv.UnstructuredGrid|pg.Mesh [None]
        Structure to plot the streamlines in to.
        If its a pv grid a check is performed if the data set is already contained.
    data: iterable [None]
        Values used for streamlining.
    label: str
        Label for the data set. Will be searched for within the data.
    radius: float [0.01]
        Radius for the streamline tubes.

    Examples
    --------
    >>> import pyvista as pv
    >>> import pygimli as pg
    >>> from pygimli.viewer.pv import drawStreamLines
    >>>
    >>> mesh = pg.createGrid(40,20,20)
    >>> data = pg.x(mesh.positions()) * pg.y(mesh.positions())
    >>>
    >>> ax, _ = pg.show(mesh, notebook=True, hold=True, alpha=0.1)
    >>> drawStreamLines(ax, mesh, data, radius=.1, source_radius=20, n_points=500)
    >>> _ = ax.show()

    Note
    ----
    All kwargs will be forwarded to pyvistas streamline filter:
    https://docs.pyvista.org/core/filters.html?highlight=streamlines#pyvista.DataSetFilters.streamlines
    """

    if isinstance(mesh, pg.Mesh):
        # create gradient of cell data
        grad = pg.solver.grad(mesh, data)
        # ensure that its point/node data in the mesh
        if len(grad) == mesh.cellCount():
            grad = pg.meshtools.cellDataToNodeData(mesh, grad)
        # add data to the mesh and convert to pyvista grid
        mesh = pgMesh2pvMesh(mesh, grad.T, label)

    elif isinstance(mesh, pv.UnstructuredGrid):
        if label not in mesh.point_arrays:  # conversion needed
            mesh.cell_data_to_point_data()

    if label is None:
        label = list(mesh.point_arrays.keys())[0]

    kwargs['vectors'] = label

    streamlines = mesh.streamlines(**kwargs)

    ax.add_mesh(streamlines.tube(radius=radius))
