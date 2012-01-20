#!/usr/bin/env python

import pygimli as g

class TestMatrix( g.MatrixBase ):
    def __init__( self ):
        g.MatrixBase.__init__( self )
        
    def rows( self ): return 1
    
    def cols( self ): return 1
    
    def mult( self, b ):
        ret = g.RVector( self.rows() )
        print "TestMatrix::mult"
        return ret
    
    def transMult( self, b ):
        ret = g.RVector( self.cols() )
        print "TestMatrix::transMult"
        return ret
        
    def save( self, name ):
        print "TestMatrix::save", name
        
class TestModelling( g.ModellingBase ):

    def __init__( self ):
        g.ModellingBase.__init__( self, True )
        self.regionManager().setParameterCount( 1 )
        self.mat = TestMatrix( )    
        
    def response( self, model ):
        print "TestModelling::response"
        res = g.RVector( 1, 1.0)
        return res
    
    def pJacobian( self ):
        print "TestModelling::jacobian()"
        return self.mat
    
    def createJacobian( self, model ):
        print "TestModelling::createJacobian"
            

    
F = TestModelling( )

dat = g.RVector( 1, 1)
err = g.RVector( 1, 0.00001)

inv = g.RInversion( dat, F, True, True )
inv.setError( err )

inv.run()
