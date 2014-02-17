#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Was macht das ding."""

import numpy as np
import pygimli as pg

from pygimli import FDEM1dModelling, RVector, asvector, RTrans, RTransLog, RTransLogLU, RInversion

def importEmsysAsciiData( filename, verbose = False ):
    """pure import function reading in positions, data, frequencies, error and
    geometry."""

    xx, sep, f, pf, ip, op, hmod, q = np.loadtxt( filename, skiprows=1, usecols=(1,4,6,8,9,12,15,16 ), unpack=True )

    err = q / pf  * 100. # percentage of primary field

    if len( np.unique( sep ) ) > 1:
        print("Warning! Several coil spacings present in file!")

    coilspacing = np.median( sep )

    f = np.round_( f )
    freq, mf, nf = np.unique( f, True, True )
    x, mx, nx = np.unique( xx, True, True )
    IP = np.ones( ( len(x), len(freq) ) ) * np.nan
    OP = np.ones( ( len(x), len(freq) ) ) * np.nan
    ERR = np.ones( ( len(x), len(freq) ) ) * np.nan

    for i in range( len( f ) ):
        #print i, nx[i], nf[i]
        IP[ nx[i],nf[i] ] = ip[i]
        OP[ nx[i],nf[i] ] = op[i]
        ERR[nx[i],nf[i] ] = err[i]

    return x, freq, coilspacing, IP, OP, ERR

def importMaxminData( filename, verbose = False ):
    """pure import function reading in positions, data, frequencies and
    geometry."""
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
            freq = np.array( [float(aa) for aa in aline[aline.find(':')+1:].replace(',',' ').split() if aa[0].isdigit()] )

    fid.close()

    if verbose: print("CS=", coilspacing, "F=", freq)
    if aline.find(',')>0: delim=','

    nf = len( freq )
    if verbose: print("delim=", delim, "nf=", nf)
    A = np.loadtxt( filename, skiprows=i, delimiter=delim ).T
    x, IP, OP = A[0], A[2:nf*2+2:2].T, A[3:nf*2+2:2].T

    return x, freq, coilspacing, IP, OP

def xfplot( ax, DATA, x, freq, everyx=5, orientation='vertical', aspect=40 ):
    """plots a matrix according to x and frequencies."""
    nt = list(range( 0, len( x ), everyx))
    plt.imshow( DATA.T, interpolation='nearest' )
    plt.ylim(plt.ylim()[::-1])
    ax.set_xticks(nt)
    ax.set_xticklabels(["%g" % xi for xi in x[nt]])
    ax.set_yticks(list(range(0,len(freq)+1,2)))
    ax.set_yticklabels(["%g" % freq[i] for i in range(0,len(freq),2)])
    plt.colorbar(orientation=orientation,aspect=aspect)
    plt.xlabel('x [m]')
    plt.ylabel('f [Hz]')


class FDEMData():
    def __init__( self, filename = None):
        """initialize data class and load data."""
        self.x, self.f, self.cs, self.IP, self.OP, self.ERR = None,None,None,None,None,None
        self.height = 1.0

        if filename:
            # check if filename extension is TXT or CSV
            if filename.lower().rfind( '.txt' ) > 0 or filename.lower().rfind( '.csv' ) > 0:
                self.x, self.f, self.cs, self.IP, self.OP, self.ERR = importEmsysAsciiData( filename )
            else:
                self.x, self.f, self.cs, self.IP, self.OP = importMaxminData( filename )

        self.activeFreq = ( self.f > 0.0 )

    def __repr__( self ):
        if isinstance( self.x, float ):
            return "<FDEMdata: sounding with %d frequencies, coilspacing is %.1f>" % (len(self.f), self.cs)
        else:
            return  "<FDEMdata: %d soundings with each %d frequencies, coilspacing is %.1f>" % (len(self.x), len(self.f), self.cs)

    def showInfos( self ): # only for old scripts using it
        print(__repr__( self ))

    def deactivate( self, fr ):
        """deactivate a single frequency."""
        fi = np.find( np.absolute( self.f / fr - 1.) < 0.1 )
        self.activeFreq[ fi ] = False

    def freq( self ):
        """return active frequencies."""
        return self.f[ self.activeFreq ]

    def FOP( self, nlay = 2 ):
        """retrieve forward modelling operator using a block discretization."""
        return FDEM1dModelling( nlay, asvector( self.freq() ), self.cs, -self.height )

    def FOPsmooth( self, zvec ):
        """retrieve forward modelling operator using fixed layers (smooth
        inversion)"""
        return FDEM1dRhoModelling( asvector( zvec ), asvector( self.freq() ), self.cs, -self.height )

    def selectData( self, xpos=0 ):
        """retrieve inphase and outphase vector from index or near given
        position."""
        if isinstance( xpos, int ) and ( xpos < len( self.x ) ) and ( xpos >= 0 ): # index
            n = xpos
        else:
            n = np.argmin( np.absolute( self.x - xpos ) )

        if self.ERR is not None:
            return self.IP[ n, self.activeFreq ], self.OP[ n, self.activeFreq ], self.ERR[ n, self.activeFreq ]
        else:
            return self.IP[ n, self.activeFreq ], self.OP[ n, self.activeFreq ], None

    def datavec( self, xpos=0 ):
        """extract data vector (stacking inphase and outphase."""
        ip, op, err = self.selectData( xpos )
        return asvector( np.hstack( ( ip, op ) ) )

    def errorvec( self, xpos=0, minvalue=0.0 ):
        """extract error vector."""
        ip, op, err = self.selectData( xpos )
        return asvector( np.tile( np.maximum( err * 0.7071, minvalue ) ) )

    def invBlock( self, xpos=0, nlay=2, noise=1.0, stmod=10., lam=100., lBound=1., uBound=0., verbose=False ):
        """yield gimli inversion instance for block inversion."""
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

    def plotData( self, xpos=0, response = None, error=None, ax=None, marker='bo-', rmarker='rx-', clf=True, addlabel='', nv=2 ):
        """plot data as curves at given position."""
        ip, op = self.selectData( xpos )
        fr = self.freq()
        if ax is None:
            if clf: plt.clf()
            plt.subplot(1,nv,nv-1)
        else:
            plt.sca( ax[0] )

        markersize = 4
        if error is not None:
            markersize = 2

        plt.semilogy( ip, fr, marker, label='obs'+addlabel, markersize=markersize )
        if error is not None and len(error) == len( ip ):
            plt.errorbar( ip, fr, xerr=error )

        plt.axis('tight')
        if error is not None:
            plt.ylim((min(fr)*.98,max(fr)*1.02))


        plt.grid(True)
        plt.xlabel('inphase [%]')
        plt.ylabel('f [Hz]')
        if response is not None:
            rip = np.asarray( response )[:len(ip)]
            plt.semilogy( rip, fr, rmarker, label='syn'+addlabel )

        plt.legend( loc='best' )

        if ax is None:
            plt.subplot(1,nv,nv)
        else:
            plt.sca( ax[1] )

        plt.semilogy( op, fr, marker, label='obs'+addlabel, markersize=markersize )
        if error is not None and len(error) == len( ip ):
            plt.errorbar( op, fr, xerr=error )

        if response is not None:
            rop = np.asarray( response )[len(ip):]
            plt.semilogy( rop, fr, rmarker, label='syn'+addlabel )

        plt.axis('tight')
        if error is not None:
            plt.ylim((min(fr)*.98,max(fr)*1.02))

        plt.grid(True)
        plt.xlabel('outphase [%]')
        plt.ylabel('f [Hz]')
        plt.legend( loc='best' )
        plt.subplot( 1, nv, 1 )
        return


    def plotDataOld( self, xpos=0, response = None, marker='bo-', rmarker='rx-', clf=True ):
        """plot data as curves at given position."""
        ip, op = self.selectData( xpos )
        fr = self.freq()
        if clf: plt.clf()
        plt.subplot(121)
        plt.semilogy( ip, fr, marker, label='obs' )
        plt.axis('tight')
        plt.grid(True)
        plt.xlabel('inphase [%]')
        plt.ylabel('f [Hz]')
        if response is not None:
            rip = np.asarray( response )[:len(ip)]
            plt.semilogy( rip, fr, rmarker, label='syn' )

        plt.legend( loc='best' )

        plt.subplot(122)
        plt.semilogy( op, fr, marker, label='obs' )
        if response is not None:
            rop = np.asarray( response )[len(ip):]
            plt.semilogy( rop, fr, rmarker, label='syn' )

        plt.axis('tight')
        plt.grid(True)
        plt.xlabel('outphase [%]')
        plt.ylabel('f [Hz]')
        plt.legend( loc='best' )
        plt.show()
        return

    def showModelAndData( self, model, xpos=0, response=None ):
        plt.clf()
        model = np.asarray( model )
        nlay = ( len( model ) + 1 ) / 2
        thk = model[:nlay-1]
        res = model[nlay-1:2*nlay-1]
        ax1 = plt.subplot(131)
        draw1dmodel( res, thk )
        ax2 = plt.subplot(132)
        ax3 = plt.subplot(133)
        self.plotData( xpos, response, (ax2, ax3), clf=False )

    def plotAllData( self, allF = True, orientation='vertical', outname=None ):
        """plot data along a profile as image plots for IP and OP."""

        freq = self.freq()
        nf = len( freq )
        np = 2
        if self.ERR is not None:
            np = 3

        plt.clf()
        ax1 = plt.subplot(np,1,1)
        xfplot( ax1, self.IP[:,self.activeFreq], self.x, freq, orientation=orientation )
        plt.title('inphase percent')
        ax2 = plt.subplot(np,1,2)
        xfplot( ax2, self.OP[:,self.activeFreq], self.x, freq, orientation=orientation )
        plt.title('outphase percent')
        if self.ERR is not None:
            ax3 = plt.subplot(np,1,3)
            xfplot( ax3, self.ERR[:,self.activeFreq], self.x, freq, orientation=orientation )
            plt.title('error percent')

        if outname is not None:
            plt.savefig( outname )

        plt.show()
        return

#    def FOP2d( nlay ):
#        """ 2d forward modelling operator """
#        return FDEM1dModelling( nlay, asvector( self.freq() ), self.cs, -self.height )
#

if __name__ == "__main__":
    import sys
    from optparse import OptionParser

    parser = OptionParser( "usage: %prog [options] fdem", version="%prog: " + pg.__version__  )
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true"
                            , help="be verbose", default=False )

    (options, args) = parser.parse_args()

    if options.verbose:
        __verbose__ = True

    if len( args ) == 0:
        parser.print_help()
        print("Please add a mesh or model name.")
        sys.exit( 2 )
    else:
        if len( args ) == 0:
            parser.print_help()
            print("Please add a mesh or model name.")
            sys.exit( 2 )
        else:
            fdem = FDEMData( args[ 0 ] )


