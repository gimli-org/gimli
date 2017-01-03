# coding=utf-8
"""
Tools to assess the quality of unstructured meshes.
---------------------------------------------------

Quality measures are based on the following review article:

Field, D. A. (2000), Qualitative measures for initial meshes. Int. J. Numer.
Meth. Engng., 47: 887–906.
"""

import pygimli as pg
import matplotlib.pyplot as plt
import numpy as np

# Helper functions
def boundaryLengths(cell):
    """Return boundary lengths of a given cell."""
    return np.array([cell.boundary(i).size() for i in range(cell.boundaryCount())])

def unitVector(vector):
    """Return the unit vector of the vector."""
    return vector / np.linalg.norm(vector)

def angleBetween(v1, v2):
    """Return the angle between vectors v1 and v2.

    Examples
    --------
    >>> angleBetween((1, 0, 0), (1, 0, 0))
    0.0
    >>> angleBetween((1, 0, 0), (-1, 0, 0))
    180.0
    """
    v1_u = unitVector(v1)
    v2_u = unitVector(v2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    return np.degrees(angle)

def cellAngles(cell):
    """Return angles of a triangular cell.

    Examples
    --------
    >>> mesh = pg.Mesh()
    >>> for pos in (0.,0.), (0.,1.), (1.,0.):
    ...     n = mesh.createNode(pos[0], pos[1], 0.0)
    >>> cell = mesh.createCell((0, 1, 2))
    >>> print(cellAngles(cell)[0])
    90.0
    """
    if cell.nodeCount() > 3:
        raise Exception("Cell %d is not a triangular cell." % cell.id())

    a = cell.node(0).pos()
    b = cell.node(1).pos()
    c = cell.node(2).pos()

    ab = b - a
    ac = c - a
    bc = c - b

    alpha = angleBetween(ab, ac)
    beta = angleBetween(-ab, bc)
    gamma = angleBetween(-ac, -bc)
    assert np.allclose(gamma, 180.0 - alpha - beta)

    return alpha, beta, gamma

# Quality measures
def eta(cell):
    r"""Return default triangle quality (eta) of a given cell.

    The quality measure relates the area of the triangle (a)
    to its edge lengths (l1, l2, l3).

    .. math::

        \eta = \frac{4\sqrt{3}a}{l_1^2 + l_2^2 + l_3^2}
    """
    return 4 * np.sqrt(3) * cell.size() / np.sum(boundaryLengths(cell)**2)

def minimumAngle(cell):
    """Return the normalized minimum angle of a given cell."""
    return np.min(cellAngles(cell))/60.

def nsr(cell):
    r"""Return the normalized shape ratio (NSR) for a given cell.

    Also referred to as the radius ratio, as it is described by
    the ratio between the inradius (r) and the circumradius (R).

    .. math::

        \rho = \frac{2r}{R}
    """
    a, b, c = boundaryLengths(cell)
    r = 2 * cell.size() / (a + b + c) # inradius
    R = 0.25 * a * b * c / cell.size() # circumradius

    return 2 * r / R

# Main function
def quality(mesh, measure="eta", show=False):
    """Return the quality of a given triangular mesh.

    Parameters
    ----------
    mesh : mesh object
        Mesh for which the quality is calculated.
    measure : quality measure, str
        Can be either "eta", "nsr", or "minimumAngle".
    show : boolean
        Show mesh quality and corresponding histogram.

    Examples
    --------
    >>> # no need to import matplotlib
    >>> import pygimli as pg
    >>> from pygimli.meshtools import polytools as plc
    >>> from pygimli.meshtools import quality
    >>> # Create Mesh
    >>> world = plc.createWorld(start=[-10, 0], end=[10, -10],
    ...                         marker=1, worldMarker=False)
    >>> c1 = plc.createCircle(pos=[0.0, -5.0], radius=3.0, area=.3)
    >>> mesh = pg.meshtools.createMesh([world, c1], quality=21.3)
    >>> # Compute and show quality
    >>> q = quality(mesh, measure="nsr")
    >>> ax, _ = pg.show(mesh, q, cmap="RdYlGn", grid=True, cMin=0.5, cMax=1.0,
    ...                 label="Normalized shape ratio")
    """

    # Implemented quality measures (Triangular meshses only.)
    measures = {
        "eta": eta,
        "nsr": nsr,
        "minimumAngle": minimumAngle
    }

    m = measures[measure]
    qualities = [m(cell) for cell in mesh.cells()]

    if show:
        fig, axes = plt.subplots(1,2)
        axes[1].hist(qualities, color="grey")
        pg.show(mesh, qualities, ax=axes[0], cMin=0.5, cMax=1, hold=True,
                logScale=False, label="Mesh quality", cmap="RdYlGn", grid=True)
        s = "Min: %.2f, Mean: %.2f, Max: %.2f" % (np.min(qualities),
                                                  np.mean(qualities),
                                                  np.max(qualities))
        axes[1].set_title(s)
        axes[1].set_xlabel("Mesh quality")
        axes[1].set_ylabel("Frequency")

        # Figure resizing according to mesh dimesions
        x = mesh.xmax() - mesh.xmin()
        y = mesh.ymax() - mesh.ymin()
        width, height = fig.get_size_inches()
        fig.set_figwidth(width * 1.1 * (x/y))
    else:
        return qualities
