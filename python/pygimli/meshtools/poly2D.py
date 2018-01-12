# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 15:48:00 2015.

@author: Marcus

No official maintenance by the GIMLi team.

"""

import numpy as np
# do not use this here import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.meshtools import createMesh, createParaMeshPLC

from pygimli.mplviewer.meshview import drawMesh
from pygimli.io import opt_import
from pygimli.meshtools import writePLC


class Poly2D(object):
    """
    Utility class to create polygon model that can be meshed and used
    for synthetic modelling.

    ca: This class can currently read a given geometry from
    xml file and create a mesh.
    """

    def __init__(self, polyfile):
        """Read a file that specifies topography and subsurface regions.

        The polygons should be specified by an xml file
        """
        self.sensorsPos = None
        self.poly = None
        self.paramaxcellsize = 0
        self.load(polyfile)

    def __str__(self):
        pass

    def __repr__(self):
        pass

    def load(self, polyfile):
        """Read polygon info from XML file (see example.xml)."""
        ET = opt_import("xml.etree.cElementTree", "read in XML files")
        self.doc = ET.parse(polyfile)
        self.parse()

    def parse(self):
        """Parses all polygon regions."""
        root = self.doc.getroot()

        for region in root.findall("line"):
            if region.get('name') == 'surface':
                self.parseSurface(region)
            else:
                self.parseRegion(region)

    def extractPoints(self, points):
        """
        Converts the text from the list of XML 'Elements' to a list of
        RVector3's that can be used for creating the polygon.

        Will create points that have x and y coordinates with z=0 always.
        """

        positions = np.empty(len(points), dtype=pg.RVector3)

        idx = 0
        for p in points:
            x = float(p.findtext('x'))
            y = float(p.findtext('y'))
            z = float(p.findtext('z'))

            if z != 0.0:
                positions[idx] = pg.RVector3(x, z, 0.0)
            else:
                positions[idx] = pg.RVector3(x, y, 0.0)
            idx += 1

        return positions

    def extractMarker(self, region):
        """Extract positions and marker id from a region and return them."""

        marker_id = int(region.get('id'))
        pos = self.extractPoints(region.findall("point[@type='marker']"))[0]

        return pos, marker_id

    def parseSurface(self, surf):
        """Parses the surface and creates the parametric domain.

        The parameter domain will have marker '2' and the outer region will
        have marker '1'.

        TODO: Add override capability via **kwargs
        """

        self.sensorsPos = self.extractPoints(
            surf.findall("point[@type='node']"))
        paradx = float(surf.get('paradx', 0.5))
        paradepth = float(surf.get('paradepth', 0.0))
        paraboundary = float(surf.get('paraboundary', 2.0))
        boundary = float(surf.get('boundary', -1.0))
        self.paramaxcellsize = float(surf.get('maxcellsize', 0.0))

        self.poly = createParaMeshPLC(self.sensorsPos, paradx, paradepth,
                                      paraboundary, self.paramaxcellsize,
                                      boundary)

        marker_pos, marker_id = self.extractMarker(surf)
        self.poly.addRegionMarker(marker_pos, marker_id)

    def parseRegion(self, region):
        """Parse a region.

        Takes care of closing regions if needed and sets
        a marker according to what is specified in the XML.
        """

        name = region.get('name')
        region_type = region.get('type')  # closed loop or not
        boundarycond = self.getBcMarkerType(region.get('bc'))
        points = region.findall("point[@type='node']")
        node_ids = np.empty(len(points), dtype=int)

        print(region)
        print('name: {}'.format(name))
        print("num points: {}".format(len(points)))

        point_pos = self.extractPoints(points)

        idx = 0
        for p in point_pos:
            node = self.poly.createNode(p)

            print('x: {}, y: {}, z: {}'.format(p.x(), p.y(), p.z()))

            node_ids[idx] = node.id()
            idx += 1

        for n in range(len(points)-1):
            self.poly.createEdge(self.poly.node(int(node_ids[n])),
                                 self.poly.node(int(node_ids[n+1])),
                                 int(boundarycond))

        if region_type == 'closed':
            self.poly.createEdge(self.poly.node(int(node_ids[-1])),
                                 self.poly.node(int(node_ids[0])),
                                 int(boundarycond))

            marker_pos = pg.RVector3()
            for i in node_ids:
                marker_pos += self.poly.node(i).pos()
            marker_pos /= float(len(node_ids))

        if region_type == 'layer':
            marker_pos, _ = self.extractMarker(region)
            # mk = region.findall("point[@type='marker']")
            # marker_pos = self._extract_points(mk)[0]

            # marker = int(region.find("point[@type='marker']").findtext('id'))
        marker = int(region.get('id'))
        print(marker)
        if marker > 0:
            self.poly.addRegionMarker(marker_pos, marker, self.paramaxcellsize)

    def getBcMarkerType(self, bc_marker):
        """Returns a pygimli marker."""

        if bc_marker == "hom_neumann":
            pg_marker = pg.MARKER_BOUND_HOMOGEN_NEUMANN
        elif bc_marker == "mixed":
            pg_marker = pg.MARKER_BOUND_MIXED
        elif bc_marker == "hom_dirichlet":
            pg_marker = pg.MARKER_BOUND_HOMOGEN_DIRICHLET
        elif bc_marker == "zero":
            pg_marker = 0
        else:
            raise NotImplementedError(
                "Boundary marker '{}' not implemented!".format(bc_marker))
        return pg_marker

    def createMesh(self, only_pd_mesh=False, verbose=False, **kwargs):
        """Generate a mesh from the polygon.

        Use the supplied quality if not already specified in the XML.

        TODO: Allow overrides through **kwargs
        """

        if not hasattr(self, 'poly'):
            raise AttributeError("No polygon created yet!")

        quality = kwargs.pop("quality", 34.2)
        smooth = kwargs.pop("smooth", (1, 10))
        surf = self.doc.getroot().find("line[@name='surface']")
        qual = float(surf.get("meshquality", quality))
        area = kwargs.pop('area', 0.0)

        m_with_bg = createMesh(self.poly, qual, smooth=smooth,
                               node_move=False, area=area, verbose=verbose)

        if only_pd_mesh:
            pd = pg.Mesh(2)
            marker_start = 2
            marker_end = -1
            pd.createMeshByMarker(m_with_bg, marker_start, marker_end)
            return pd

        return m_with_bg

    def show(self, ax=None, **kwargs):
        """Diplays the polygonal model using pg.show()."""

        import matplotlib.pyplot as plt
        if ax is None:
            figkeys = ('nrows', 'ncols', 'sharex', 'sharey', 'figsize')
            figargs = dict((k, kwargs.pop(k)) for k in figkeys if k in kwargs)
            _, ax = plt.subplots(**figargs)

        drawMesh(ax, self.poly, **kwargs)

    def export(self, filename):
        """Export Triangle poly file."""
        writePLC(self.poly, filename)

    def drawLines(self, ax=None, color='white'):
        """Draw edges of a mesh as lines (e.g. onto an inversion result)."""
        for b in self.poly.boundaries():
            ax.plot([b.node(0).x(), b.node(1).x()],
                    [b.node(0).y(), b.node(1).y()], color=color)

if __name__ == '__main__':
    plc = Poly2D('example.xml')  # Poly2D handles loading/creating a 2D-PLCs
    plc.show()
    mesh = plc.createMesh()
    pg.show(mesh, mesh.cellMarkers())
