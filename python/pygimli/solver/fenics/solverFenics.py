#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    http://fenicsproject.org/
"""

import pygimli as pg

dolfin = None
dlf = None
try:
    import dolfin
    dlf = dolfin
    from dolfin_utils.meshconvert import xml_writer
except:
    BaseException("Cannot find fenics packages")


def convertDolfinMesh(mesh):
    """
        Convert :term:`Dolfin` mesh instance into a gimliapi:`GIMLI::Mesh`
        
        The resulting mesh does not separate between different boundary domains. 
        Instead each outer boundary have the marker 1.
        
        TODO 
        * 1D Mesh
        * 3D Mesh
                
    Parameters
    ----------
    mesh: Dolfin::Mesh
    
    Returns
    -------
    out: gimliapi:`GIMLI::Mesh`
    
    """
    out = pg.Mesh(2)
    for v in mesh.coordinates(): out.createNode(v)
    for c in mesh.cells(): out.createCell(pg.IndexArray(c))
    
    out.createNeighbourInfos()
    for b in out.boundaries():
        if b.leftCell() is None or b.rightCell() is None:
            b.setMarker(1)
            
    return out

def exportDolfinXML(mesh, ofile):
    """
        Write gimliapi:`GIMLI::Mesh` mesh into dolfin xml format
        
          TODO 
        * 1D Mesh
        * 3D Mesh
        * quad meshs
        * mixed meshs
        * tests
        
    """
    # Write everything out
    xml_writer.write_header_mesh(ofile, "triangle", 2)
    xml_writer.write_header_vertices(ofile, mesh.nodeCount())
    
    
    for n in mesh.nodes():
        xml_writer.write_vertex(ofile, n.id(), n.pos()[0], n.pos()[1], 0.0)
    
    xml_writer.write_footer_vertices(ofile)
    xml_writer.write_header_cells(ofile, mesh.cellCount())
    
    for c in mesh.cells():
        if c.nodeCount() == 4:
            raise BaseException('pls support quads')
        xml_writer.write_cell_triangle(ofile, c.id(), 
                                       c.node(0).id(),
                                       c.node(1).id(),
                                       c.node(2).id()
                                       )
    xml_writer.write_footer_cells(ofile)
    
    #edge_markers_local = 1
    #if len(edge_markers_local) > 0:
        #xml_writer.write_header_domains(ofile)
        #xml_writer.write_header_meshvaluecollection(ofile, \
                            #"edge markers", 1, len(edge_markers_local), "uint")
        #for tri, local_edge, marker in edge_markers_local:
             #xml_writer.write_entity_meshvaluecollection(ofile, \
                                            #1, tri+tri_off, marker, local_edge)
        #xml_writer.write_footer_meshvaluecollection(ofile)
        #xml_writer.write_footer_domains(ofile)
    
    xml_writer.write_footer_mesh(ofile)
    
    for i in range(i):
        afilename = ofilename.replace(".xml", ".attr"+str(i)+".xml")
        afile = open(afilename, "w")
        xml_writer.write_header_meshfunction2(afile)
        xml_writer.write_header_meshvaluecollection(afile, \
                             "triangle attribs "+str(i), 2, mesh.cellCount, "double")
        for c in mesh.cells():
             xml_writer.write_entity_meshvaluecollection(afile, \
                                            2, c.id(), c.attribute(), 0)
        xml_writer.write_footer_meshvaluecollection(afile)
        xml_writer.write_footer_meshfunction(afile)
        print("triangle attributes from .ele file written to "+afilename)
        afile.close()

    # Close files
    ofile.close()

if __name__ == '__main__':

    pass