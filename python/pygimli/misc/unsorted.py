# -*- coding: utf-8 -*-

import pygimli as g

def rot2DGridToWorld( mesh, start, end ):
    print mesh, start, end
    mesh.rotate( g.degToRad( g.RVector3( -90.0, 0.0, 0.0 ) ) )

    src = g.RVector3( 0.0, 0.0, 0.0 ).norm( g.RVector3( 0.0,  0.0, -10.0 ), g.RVector3( 10.0, 0.0, -10.0 ) )
    dest = start.norm( start - g.RVector3( 0.0, 0.0, 10.0 ), end )

    q = g.getRotation( src, dest )
    rot = g.RMatrix( 4, 4 )
    q.rotMatrix( rot )
    mesh.transform( rot )
    mesh.translate( start )


def streamline( mesh, field, start, dLength, maxSteps = 1000, verbose = False, koords=[0,2] ):
    xd = []
    yd = []
    counter = 0
    
    # search downward
    pos = g.RVector3( start )
    c = mesh.findCell( pos );
    lastU = 1e99;
    
    while c is not None and len( xd ) < maxSteps:
        d = c.grad( pos, field )
        u = c.interpolate( pos, field )
        #print "cell:", c.id(), u
        # always go u down
        if u > lastU:
            #print u, lastU
            break;
            #pass
        
        pos -= d/d.length() * dLength * max( 1.0, ( (start-pos).length() ) )
        xd.append( pos[ koords[ 0 ] ] )
        yd.append( pos[ koords[ 1 ] ] )
        c = mesh.findCell( pos, False );
        
        lastU = u
        if verbose:
            print pos, u

    xu = []
    yu = []
    #return xu, yu
    # search upward
    pos = g.RVector3( start )
    c = mesh.findCell( pos );

    lastu=-1e99
    while c is not None and len( xu ) < maxSteps:
        d = c.grad( pos, field )
        u = c.interpolate( pos, field )
        
        # always go u up
        if u < lastU:
            break;
    
        pos += d/d.length() * dLength * max( 1.0, ( (start-pos).length() ) )
        xu.append( pos[ koords[ 0 ] ] )
        yu.append( pos[ koords[ 1 ] ] )
        c = mesh.findCell( pos, False );
        lastU = u
            
    xu.reverse()
    yu.reverse()
    x = xu + xd
    y = yu + yd
    return x,y

def boundaryPlaneIntersectionLines( boundaries, plane ):
    '''
        Create Lines from boundaries that intersect a plane
    '''
    lines = []

    for b in boundaries:
        ps = []
        for i, n in enumerate( b.shape().nodes() ):
            line = g.Line( n.pos(), b.shape().node( (i+1)%b.shape().nodeCount() ).pos() )
            p = plane.intersect( line, 1e-8, True )
            if p.valid():
                ps.append( p )
        
        if len( ps ) == 2:
            lines.append( zip( [ ps[0].x(), ps[1].x() ],
                               [ ps[0].z(), ps[1].z() ] ) )
    return lines
# def intersectionLines



def number_of_processors():
    import sys, os
    ''' number of virtual processors on the computer '''
    # Windows
    if os.name == 'nt':
        return int(os.getenv('NUMBER_OF_PROCESSORS'))
    # Linux
    elif sys.platform == 'linux2':
        retv = 0
        with open('/proc/cpuinfo','rt') as cpuinfo:
            for line in cpuinfo:
                if line[:9] == 'processor': retv += 1
        return retv

    # Please add similar hacks for MacOSX, Solaris, Irix,
    # FreeBSD, HPUX, etc.
    else:
        raise RuntimeError, 'unknown platform'
    
    
def assembleDC( mesh, source = g.RVector3( 0.0, 0.0, 0.0 ) ):
    '''
        assemble stiffness matrix for 3d dc forward problem using fem
    '''
    S = g.DSparseMatrix()
    S.buildSparsityPattern( mesh )
#    se = g.DElementMatrix()

#    for c in mesh.cells():
#        se.ux2uy2uz2( c )
#        S += se
    g.dcfemDomainAssembleStiffnessMatrix( S, mesh, 0.0, False )
    g.dcfemBoundaryAssembleStiffnessMatrix( S, mesh, source, 0.0 )
    return S
#def assembleDC
    
def assembleCEM( S, mesh, marker, zi, nodeID = -1, verbose = False ):
    '''
        add dc-cem to stiffness system, return new Matrix and sum of electrodes surface
    '''

    sumArea = 0

    if nodeID == -1:
        for b in mesh.findBoundaryByMarker( marker ):
            sumArea += b.shape().domainSize()
        print "addCEM: ", marker, sumArea, 'm^2', zi/sumArea, 'Ohm\n', 
    else:
        sumArea = 1
        print "addCEM: node"

    mapS = g.DSparseMapMatrix( S );
    oldSize = S.size()

    se = g.DElementMatrix()
    
    mapS.setRows( oldSize + 1 );
    mapS.setCols( oldSize + 1 );
        
    if nodeID == -1:
        for b in mesh.findBoundaryByMarker( marker ):
            se.u( b )
            se /= -zi
            mapS.addToCol( oldSize, se );
            mapS.addToRow( oldSize, se );
            se.u2(b)
            se /= zi
            mapS += se

        mapS.setVal( oldSize,  oldSize, sumArea /zi )
    else:
        mapS.addVal( nodeID, nodeID,  1.0 )
        mapS.addVal( oldSize, nodeID,  - 1.0 )
        mapS.addVal( nodeID,  oldSize, - 1.0 )
        mapS.addVal( oldSize,  oldSize, 1.0 )

    return g.DSparseMatrix( mapS ), sumArea
#def assembleCEM