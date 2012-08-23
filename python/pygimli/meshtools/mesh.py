# -*- coding: utf-8 -*-

import pygimli as g
import pygimli.utils

def readHydrus2dMesh(filename='MESHTRIA.TXT'):
    '''
        import Hydrus2D mesh (MESHTRIA.TXT)
        mesh=readHydrus2dMesh(filename) 
    '''
    fid=open(filename)
    line=fid.readline().split()
    nnodes=int(line[1])
    ncells=int(line[3])
    mesh=g.Mesh()
    for i in range(nnodes):
        line=fid.readline().split()
        mesh.createNode(g.RVector3(float(line[1])/100.,float(line[2])/100.,0.))
    
    for i in range(3): 
        line=fid.readline()
        
    for i in range(ncells):
        line=fid.readline().split()
        if len(line)==4:
            mesh.createTriangle(mesh.node(int(line[1])-1),mesh.node(int(line[2])-1),mesh.node(int(line[3])-1),1)
        elif len(line)==5:
            mesh.createTetrahedron(mesh.node(int(line[1])-1),mesh.node(int(line[2])-1),
                                   mesh.node(int(line[3])-1),mesh.node(int(line[4])-1),1)
    
    fid.close()
    return mesh
    

def readHydrus3dMesh(filename='MESHTRIA.TXT'):
    '''
        import regular Hydrus3D mesh (MESHTRIA.TXT)
        mesh=readHydrus3dMesh(filename) 
    '''
    f=open(filename,'r')
    for i in range(6):
        line1=f.readline()
        
    nnodes=int(line1.split()[0])
    ncells=int(line1.split()[1])
    print nnodes, ncells
    line1=f.readline()
    nodes=[]
    dx=0.01
    mesh=g.Mesh()
    for ni in range(nnodes):
        pos=f.readline().split()
        p=g.RVector3(float(pos[1])*dx,float(pos[2])*dx,float(pos[3])*dx*(-1.))
        n=mesh.createNode(p)
        nodes.append(n)

    line1=f.readline()
    line1=f.readline()
    cells=[]
    for ci in range(ncells):
        pos=f.readline().split()
        i,j,k,l=int(pos[1]),int(pos[2]),int(pos[3]),int(pos[4]),
        c=mesh.createTetrahedron(nodes[i-1],nodes[j-1],nodes[k-1],nodes[l-1])
        cells.append(c)

    f.close()
    return mesh

def rot2DGridToWorld( mesh, start, end ):
    '''
        rotate a given 2D grid in 
    '''
    mesh.rotate( g.degToRad( g.RVector3( -90.0, 0.0, 0.0 ) ) )

    src = g.RVector3( 0.0, 0.0, 0.0 ).norm( g.RVector3( 0.0,  0.0, -10.0 ),
                                            g.RVector3( 10.0, 0.0, -10.0 ) )
    dest = start.norm( start - g.RVector3( 0.0, 0.0, 10.0 ), end )
                        
    q = g.getRotation( src, dest )
    rot = g.RMatrix( 4, 4 )
    q.rotMatrix( rot )
    mesh.transform( rot )
    mesh.translate( start )
#def rot2DGridToWorld( ... )


def merge2Meshes( m1, m2 ):
    '''
        Merge two meshes into one a new mesh
        return the new mesh.
    '''
    
    mesh = g.Mesh( m1 )

    for c in m2.cells():
        mesh.copyCell( c )

    for key in mesh.exportDataMap().keys():
        d = mesh.exportDataMap()[ key ]
        print d
        d.resize( mesh.cellCount() )
        d.setVal( m1.exportDataMap()[ key ], 0, m1.cellCount() )
        d.setVal( m2.exportDataMap()[ key ], m1.cellCount(), m1.cellCount() + m2.cellCount() )
        mesh.addExportData( key, d )
        
    return mesh;
    
def mergeMeshes( meshlist ):
    '''
        Merge several meshes into one a new mesh, meshlist have to be a list of at least 2 meshes
        return the new mesh.
    '''
    if type( meshlist ) is not list:
        raise "argument meshlist is no list"
    
    if len( meshlist ) < 2:
        raise "to few meshes in meshlist"
    
    mesh = meshlist[ 0 ]
    
    for m in range( 1, len( meshlist ) ):
        mesh = merge2Meshes( mesh, meshlist[ m ] )
        
    return mesh
    
        
def createParaMesh2dGrid( sensors, paraDX = 1, paraDZ = 1, paraDepth = 0, nLayers = 11, boundary = -1, paraBoundary = 2, verbose = False,  *args, **kwargs ):
    '''
        Create gimli parameter mesh for a given list of sensor positions
        sensors .. list of RVector3
        paraDX - Horizontal distance between sensors, relative regarding sensor distance. Value must be greater than 0 otherwise 1 is assumed
        paraDZ - Vertical distance to the first depth layer, relative regarding sensor distance.
        Value must be greater than 0 otherwise 1 is assumed
        paraDepth - Maximum depth for parametric domain, 0[default] means 0.4 * max sensor range
        nLayers - Number of depth layers
        boundary - triangle boundary to be appended for domain prolongation  values lover 0 force boundary to be 4 times para domain width
        paraBoundary - offset for parameter domain in absolute sensor distance 2[default]
        verbose .. be verbose
    '''
    mesh = g.Mesh( 2 )
    
    # maybee separete x y z and sort
    sensorX = g.x( sensors )
    eSpacing = abs( sensorX[ 1 ] - sensorX[ 0 ] )

    xmin = min( sensorX ) - paraBoundary * eSpacing
    xmax = max( sensorX ) + paraBoundary * eSpacing

    if paraDX == 0: paraDX = 1.
    if paraDZ == 0: paraDZ = 1.
    
    dx = eSpacing * paraDX
    dz = eSpacing * paraDZ
    
    if paraDepth == 0:
        paraDepth = 0.4 * (xmax-xmin)
        
    x = g.utils.grange( xmin, xmax, dx = dx )
    
    y = -g.increasingRange( dz, paraDepth, nLayers )

    mesh.createGrid( x, y )
        
    map( lambda cell: cell.setMarker( 2 ), mesh.cells() )

    paraXLimits = [ xmin, xmax ]
    paraYLimits = [ min( y ), max( y ) ]
    
    if boundary < 0:
        boundary = abs( ( paraXLimits[ 1 ] - paraXLimits[ 0 ] ) * 4.0 )
        
    mesh = g.meshtools.appendTriangleBoundary( mesh, xbound = boundary, ybound = boundary, marker = 1, 
                                                *args, **kwargs )
  
    return mesh
# def createParaMesh2dGrid( ... )
    