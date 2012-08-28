# -*- coding: utf-8 -*-

import pygimli as g

def assembleStiffnessMatrixHomogenDirichletBC( S, u0, rhs = None ):
    '''
        Create homogeneous dirichlet boundary condition and apply them to stiffness matrix S
    '''
    
    # start interpreting homogeneous Dirichlet boundary condition
    u0_Ids = []
    
    if len( u0 ):
       
        if ( type( u0 ) == g._pygimli_.stdVectorBounds ):
            #assuming 
            for b in u0:
                for n in b.nodes():
                    u0_Ids.append( n.id() )
        else:
            raise Exception( "cannot interpret u0 sequence: " + u0 )
            
    for i in u0_Ids:
        S.cleanRow( i )
        S.cleanCol( i )
        S.setVal( i, i, 1.0 )
        
        if rhs is not None:
            rhs[ i ] = 0.
# def assembleStiffnessMatrixHomogenDirichletBC( ... )
    
def assembleStiffnessMatrixDirichletBC( S, boundaries, rhs ):
    u_Ids = []
    
    if len( boundaries ):
       
        if ( type( boundaries ) == g._pygimli_.stdVectorBounds ):
            #assuming 
            for b in boundaries:
                for n in b.nodes():
                    u_Ids.append( n.id() )
        else:
            raise Exception( "cannot interpret u0 sequence: " + boundaries )
            
    ud = g.RVector( S.rows(), 0.0 )
    
    for i in u_Ids:
    
        #S.cleanRow( i )
        #S.cleanCol( i )
        #S.setVal( i, i, 1.0 )
        
        ud[ i ] = 1.

    print ud
    rhs -= S * ud
        
    #for i in u_Ids:
    
        #S.cleanRow( i )
        #S.cleanCol( i )
        #S.setVal( i, i, 1.0 )
        
        #rhs[ i ] = 0.
        
    print "rhs-3", rhs
    
    
def solvePoisson( mesh, a = None, f = 1, u0 = None, Btest = None, verbose = False ):
    '''
    '''
    if verbose:
        print( mesh )

    swatch = g.Stopwatch( True )
    
    # define an empty stiffness matrix
    S = g.DSparseMatrix()
    
    # create matrix structure regarding the mesh
    S.buildSparsityPattern( mesh )

    # define a local stiffness matrix 
    Se = g.DElementMatrix()

    # check for material parameter
    if a is None:
        a = g.RVector( mesh.cellCount(), 1.0 )
    
    if len( a ) != mesh.cellCount():
        raise Exception( "Material array 'a' has the wrong size: " + len( a ) + " != " +  mesh.cellCount() )
    
    # assemble the stiffness matrix
    rhs = g.RVector( mesh.nodeCount(), 0 )
    
    for c in mesh.cells():
        Se.ux2uy2uz2( c )
        Se *= a[ c.id() ] 
        S += Se
    
        Se.u( c )
        for i, idx in enumerate( Se.idx() ):
            rhs[ idx ] += Se.row( 0 )[ i ] * f
    
    
    if Btest:
        ud = assembleStiffnessMatrixDirichletBC( S, Btest, rhs )    
    
    assembleStiffnessMatrixHomogenDirichletBC( S, u0, rhs )
           
    u = g.RVector( rhs.size(), 0.0 )
    
    if verbose:
        print( "asssemblation takes: ", swatch.duration( True ) )

    print S, S.rows(), S.cols()
    print rhs, u
        
    solver = g.LinSolver( S, verbose )
    solver.solve( rhs, u )
    
    
    if verbose:
        print( "lin solving takes: ", swatch.duration( True ) )
    
    print min(u), max(u)
    
    return u
