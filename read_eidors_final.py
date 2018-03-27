# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 15:00:47 2018

@author: srao
"""

import numpy as np
import scipy.io as spio

import pygimli as pg


# import matlab structure ----------------------------------------------------


def loadmat(filename):
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return check_keys(data)


def check_keys(dict):
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = todict(dict[key])
    return dict


def todict(matobj):
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = todict(elem)
        else:
            dict[strg] = elem
    return dict


def get_nested(data, *args):
    if args and data:
        element = args[0]
        if element:
            value = data.get(element)
            return value if len(args) == 1 else get_nested(value, *args[1:])


#-----------------------------------------------------------------------------


def read_eidors(filename, matlab_varname):
    """Reads finite element model in EIDORS format and returns pygimli mesh.

    Parameters
    ----------
    filename : str
        name of the .mat file containing the EIDORS model
    matlab_varname : str
        variable name of .mat file in MATLAB workspace
    """

    print('Reading %s... ' % filename)
    matlab_eidors = loadmat(filename)
    python_eidors = get_nested(matlab_eidors, matlab_varname)
    # if input eidors data is forward model instead of an image
    if 'nodes' in python_eidors.keys():
        nodes = get_nested(python_eidors, "nodes")
        elems = get_nested(python_eidors, "elems")
        boundary_numbers = get_nested(python_eidors, "boundary_numbers")

    # it is an image with elem_data and fwd_model
    else:
        nodes = get_nested(python_eidors, "fwd_model", "nodes")
        elems = get_nested(python_eidors, "fwd_model", "elems")
        boundary_numbers = get_nested(python_eidors, "fwd_model",
                                      "boundary_numbers")

    dim_nodes = np.size(nodes, 1)
    dim_elems = np.size(elems, 1)

    print('found %s  %s-dimensional nodes... ' % (np.size(nodes, 0), dim_nodes))
    if dim_elems == 3:
        print('found %s triangles... ' % np.size(elems, 0))
    elif dim_elems == 4:
        print('found %s tetrahedrons... ' % np.size(elems, 0))

    if 'elem_data' in python_eidors.keys():
        elem_data = get_nested(python_eidors, "elem_data")
        print('found %s element data... ' % len(elem_data))
    else:
        print("found no element data... ")

    region_markers = np.unique(boundary_numbers)
    no_of_regions = len(region_markers)

    if boundary_numbers is not None:
        print('Found %s Unique Regions with Region Markers' % no_of_regions)
        print('%s' % region_markers)

    if boundary_numbers is None:
        print('Found no Unique Region with Region Markers')

    # create nodes from eidors model
    mesh = pg.Mesh()

    # if two dimensional eidors model
    if (dim_nodes == 2 and dim_elems == 3):
        print('converting to %d-D pygimli mesh... ' % dim_nodes)
        for i in range(len(nodes)):
            mesh.createNode(pg.RVector3(nodes[i, 0], nodes[i, 1], 0.))
        for i in range(len(elems)):
            mesh.createTriangle(
                mesh.node(int(elems[i, 0]) - 1),
                mesh.node(int(elems[i, 1]) - 1),
                mesh.node(int(elems[i, 2]) - 1), 1)

    # for non-planar 2D models
    if (dim_nodes == 3 and dim_elems == 3):
        print('converting to pygimli mesh...')
        print('found 3d nodes with 2d elements')
        for i in range(len(nodes)):
            mesh.createNode(pg.RVector3(nodes[i, 0], nodes[i, 1], nodes[i, 2]))
        for i in range(len(elems)):
            mesh.createTriangle(
                mesh.node(int(elems[i, 0]) - 1),
                mesh.node(int(elems[i, 1]) - 1),
                mesh.node(int(elems[i, 2]) - 1), 1)

    # if three dimensional eidors model
    if (dim_nodes == 3 and dim_elems == 4):
        print('converting to %d-D pygimli mesh... ' % dim_nodes)
        for i in range(len(nodes)):
            mesh.createNode(pg.RVector3(nodes[i, 0], nodes[i, 1], nodes[i, 2]))
        for i in range(len(elems)):
            mesh.createTetrahedron(
                mesh.node(int(elems[i, 0]) - 1),
                mesh.node(int(elems[i, 1]) - 1),
                mesh.node(int(elems[i, 2]) - 1),
                mesh.node(int(elems[i, 3]) - 1), 1)

    if boundary_numbers is not None:
        mesh.setCellMarkers(boundary_numbers.astype(int))

    return mesh


#-------- test the function --------------------------------------------------

filename = 'A.mat'
matlab_varname = "A"
mesh = read_eidors(filename, matlab_varname)

pg.show(mesh, mesh.cellMarkers())
