#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
    Was macht das ding

'''

import pygimli as g

import pylab as P
from pygimli import FDEM1dModelling, RVector, asvector, RTrans, RTransLog, RTransLogLU, RInversion

def importMaxminData( filename, verbose = False ):
    """ pure import function reading in positions, data, frequencies and geometry """
    delim = None
    fid = open(filename)
    coilspacing = 0.
    freq = []
    for i, aline in enumerate( fid ):
        if aline.split()[0][0].isdigit(): #number found
            break
        elif aline.find('COIL') > 0:     #[:6] == '/ COIL':
            coilspacing = float( aline.split()[-2] )
        elif aline.find('FREQ') > 0:   #[:6] == '/ FREQ':
            freq = P.array( [float(aa) for aa in aline[aline.find(':')+1:].replace(',',' ').split() if aa[0].isdigit()] )
    
    fid.close()
    
    if verbose: print "CS=", coilspacing, "F=", freq
    if aline.find(',')>0: delim=','
    
    nf = len( freq )
    if verbose: print "delim=", delim, "nf=", nf
    A = P.loadtxt( filename, skiprows=i, delimiter=delim ).T
    x, IP, OP = A[0], A[2:nf*2+2:2].T, A[3:nf*2+2:2].T

    return x, freq, coilspacing, IP, OP

class FDEMData():
    def __init__( self, filename = None):
        """ initialize data class and load data """
        self.x, self.f, self.cs, self.IP, self.OP = None,None,None,None,None
        
        if filename:
            self.x, self.f, self.cs, self.IP, self.OP = importMaxminData( filename )
            
        self.height = 1.0
        self.activeFreq = ( self.f > 0.0 )
        
    def showInfos( self ):
        if isinstance( self.x, float ):
            print "Soundings with", len(self.f), "frequencies"
        else:
            print len(self.x), "soundings with each", len(self.f), "frequencies"

    def deactivate( self, fr ):
        """ deactivate a single frequency """
        fi = P.find( P.absolute( self.f / fr - 1.) < 0.1 )
        self.activeFreq[ fi ] = False
        
    def freq( self ):
        """ return active frequencies """
        return self.f[ self.activeFreq ]
    
    def FOP( self, nlay = 2 ):
        """ retrieve forward modelling operator """
        return FDEM1dModelling( nlay, asvector( self.freq() ), self.cs, -self.height )
    
    def selectData( self, xpos=0 ):
        """ retrieve inphase and outphase vector from index or near given position """
        if isinstance( xpos, int ) and ( xpos < len( self.x ) ) and ( xpos >= 0 ): # index
            n = xpos
        else:
            n = P.argmin( P.absolute( self.x - xpos ) )
        
        return self.IP[ n, self.activeFreq ], self.OP[ n, self.activeFreq ]

    def datavec( self, xpos=0 ):
        """ extract data vector (stacking inphase and outphase """
        ip, op = self.selectData( xpos )
        return asvector( P.hstack( ( ip, op ) ) )
    
    def invBlock( self, xpos=0, nlay=2, noise=1.0, stmod=10., lam=100., lBound=1., uBound=0., verbose=False ):
        """ yield gimli inversion instance for block inversion """
        """ inv(xpos,nlay) where nlay can be a FOP or a number of layers """
        self.transThk = RTransLog()
        self.transRes = RTransLogLU( lBound, uBound )
        self.transData = RTrans()
        # EM forward operator
        if isinstance( nlay, FDEM1dModelling ):
            self.fop = nlay
        else:
            self.fop = self.FOP( nlay )
        
        data = self.datavec( xpos )
        
        self.fop.region(0).setTransModel( self.transThk )
        self.fop.region(1).setTransModel( self.transRes )
        
        if isinstance( noise, float ):
            noiseVec = RVector( len(data), noise)
        else:
            noiseVec = asvector( noise )
        
        # independent EM inversion
        self.inv = RInversion( data, self.fop, self.transData, verbose )
        if isinstance( stmod, float): # real model given
            model = RVector( nlay * 2 - 1, stmod )
            model[0] = 2.
        else:
            if len( stmod ) == nlay*2-1:
                model = asvector( stmod )
            else:
                model = RVector( nlay*2-1, 30. )
        
        self.inv.setAbsoluteError( noiseVec )
        self.inv.setLambda( lam )
        self.inv.setMarquardtScheme( 0.8 )
        self.inv.setModel( model )
        self.inv.setReferenceModel( model )
        return self.inv
    
    def plotData( self, xpos=0, response = None, marker='bo-', rmarker='rx-', clf=True ):
        """ plot data as curves at given position """
        ip, op = self.selectData( xpos )
        fr = self.freq()
        if clf: P.clf()
        P.subplot(121)
        P.semilogy( ip, fr, marker, label='obs' )
        P.axis('tight')
        P.grid(True)
        P.xlabel('inphase [%]')
        P.ylabel('f [Hz]')
        if response is not None:
            rip = P.asarray( response )[:len(ip)]
            P.semilogy( rip, fr, rmarker, label='syn' )
        
        P.legend( loc='best' )
        
        P.subplot(122)
        P.semilogy( op, fr, marker, label='obs' )
        if response is not None:
            rop = P.asarray( response )[len(ip):]
            P.semilogy( rop, fr, rmarker, label='syn' )
        
        P.axis('tight')
        P.grid(True)
        P.xlabel('outphase [%]')
        P.ylabel('f [Hz]')
        P.legend( loc='best' )
        return
    
    def plotAllData( self, allF = True, orientation='vertical' ):
        """ plot data along a profile as image plots for IP and OP """
        nt = range( 0, len( self.x ), 5 )
        freq = self.freq()
        nf = len( freq )
        P.clf()
        ax1 = P.subplot(211)
        P.imshow( self.IP[:,self.activeFreq].T, interpolation='nearest' )
        P.imshow( self.IP[:,self.activeFreq].T, interpolation='nearest' )
        P.ylim(P.ylim()[::-1])
        ax1.set_xticks(nt)
        ax1.set_xticklabels(["%g" % xi for xi in self.x[nt]])
        ax1.set_yticks(range(0,nf+1,2))
        ax1.set_yticklabels(["%g" % freq[i] for i in range(0,nf,2)])
        P.colorbar(orientation=orientation,aspect=30)
        P.xlabel('x [m]')
        P.ylabel('f [Hz]')
        P.title('inphase percent')
        ax2 = P.subplot(212)
        P.imshow( self.OP[:,self.activeFreq].T, interpolation='nearest' )
        P.ylim(P.ylim()[::-1])
        ax2.set_xticks(nt)
        ax2.set_xticklabels(["%g" % xi for xi in self.x[nt]])
        ax2.set_yticks(range(0,nf+1,2))
        ax2.set_yticklabels(["%g" % freq[i] for i in range(0,nf,2)])
        P.colorbar(orientation=orientation,aspect=30)
        P.xlabel('x [m]')
        P.ylabel('f [Hz]')
        P.title('outphase percent')
        return

#    def FOP2d( nlay ):
#        """ 2d forward modelling operator """
#        return FDEM1dModelling( nlay, asvector( self.freq() ), self.cs, -self.height )
#    

if __name__ == "__main__":
    import sys
    from optparse import OptionParser

    parser = OptionParser( "usage: %prog [options] fdem", version="%prog: " + g.versionStr()  )
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true"
                            , help="be verbose", default=False )
    
    (options, args) = parser.parse_args()

    if options.verbose:
        __verbose__ = True 
        
    if len( args ) == 0:
        parser.print_help()
        print "Please add a mesh or model name."
        sys.exit( 2 )
    else:
        if len( args ) == 0:
            parser.print_help()
            print "Please add a mesh or model name."
            sys.exit( 2 )
        else:
            fdem = FDEMData( args[ 0 ] )


