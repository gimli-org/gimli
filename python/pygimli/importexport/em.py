import numpy as np
import matplotlib.pyplot as plt
from pygimli import FDEM1dModelling, RVector, asvector, RTrans, RTransLog, RTransLogLU, RInversion
from pygimli.utils import draw1dmodel

def readusffile( filename, DATA = [] ):
    ''' read data from single USF (universal sounding file) file
        DATA = readusffile( filename )
        DATA = readusffile( filename, DATA ) will append to DATA '''
    
    columns = []
    nr = 0
    sounding = {}
    sounding['FILENAME'] = filename
    isdata = False
    fid = open( filename )
    for line in fid:
        zeile = line.replace('\n','').replace(',','') # commas are useless here
        if zeile: # anything at all
            if zeile[0] == '/': # comment-like
                if zeile[1:4] == 'END': # end of a sounding
                    if isdata: # already read some data
                        sounding[ 'data' ] = columns
                        for i, cn in enumerate( sounding['column_names'] ):
                            sounding[cn] = columns[:,i]
                        
                        DATA.append( sounding )
                        sounding = {}
                    
                    isdata = not isdata # turn off data mode
                elif zeile.find(':') > 0: # key-value pair
                    key, value = zeile[1:].split(':')
                    try:
                        val = float( value )
                        sounding[key] = val
                    except:    
                        sounding[key] = value
    
            else:
                if isdata:
                    values = zeile.split()
                    try:
                        for i, v  in enumerate( values ):
                            columns[ nr, i ] = float( v )
                        
                        nr += 1
                    except:
                        sounding['column_names'] = values
                        columns = np.zeros( ( int(sounding['POINTS']), len( values ) ) )
                        nr = 0
    
    fid.close()
    return DATA

def readusffiles( filenames, DATA = [] ):
    ''' read all soundings data from a list of usf files
        DATA = readusffiles( filenames ) '''
    DATA = []
    for onefile in filenames:
        DATA = readusffile( onefile, DATA )
    
    return DATA
    
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
        elif aline.upper().find('COIL') > 0:     #[:6] == '/ COIL':
            coilspacing = float( aline.split()[-2] )
        elif aline.upper().find('FREQ') > 0:   #[:6] == '/ FREQ':
            freq = np.array( [float(aa) for aa in aline[aline.find(':')+1:].replace(',',' ').split() if aa[0].isdigit()] )
    
    fid.close()
    
    if verbose: print("CS=", coilspacing, "F=", freq)
    if aline.find(',')>0: delim=','
    
    nf = len( freq )
    if verbose: print("delim=", delim, "nf=", nf)
    A = np.loadtxt( filename, skiprows=i, delimiter=delim ).T
    x, IP, OP = A[0], A[2:nf*2+2:2].T, A[3:nf*2+2:2].T

    return x, freq, coilspacing, IP, OP

class FDEMData():
    def __init__( self, filename, height=1.0, verbose=False ):
        """initialize data class and load data."""
        self.x, self.f, self.cs, self.IP, self.OP = importMaxminData( filename, verbose )
        self.height = height
        self.activeFreq = ( self.f > 0.0 )
        
    def showInfos( self ):
        if isinstance( self.x, float ):
            print("Soundings with", len(self.f), "frequencies")
        else:
            print(len(self.x), "soundings with each", len(self.f), "frequencies")

    def deactivate( self, fr ):
        """deactivate a single frequency."""
        fi = np.find( np.absolute( self.f / fr - 1.) < 0.1 )
        self.activeFreq[ fi ] = False
        
    def freq( self ):
        """return active frequencies."""
        return self.f[ self.activeFreq ]
    
    def FOP( self, nlay = 2 ):
        """retrieve forward modelling operator."""
        return FDEM1dModelling( nlay, asvector( self.freq() ), self.cs, -self.height )
    
    def selectData( self, xpos=0 ):
        """retrieve inphase and outphase vector from index or near given
        position."""
        if isinstance( xpos, int ) and ( xpos < len( self.x ) ) and ( xpos >= 0 ): # index
            n = xpos
        else:
            n = Pnpargmin( np.absolute( self.x - xpos ) )
        
        return self.IP[ n, self.activeFreq ], self.OP[ n, self.activeFreq ]

    def datavec( self, xpos=0 ):
        """extract data vector (stacking inphase and outphase."""
        ip, op = self.selectData( xpos )
        return asvector( np.hstack( ( ip, op ) ) )
    
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
    
    def plotData( self, xpos=0, response = None, ax=None, marker='bo-', rmarker='rx-', clf=True, addlabel='', nv=2 ):
        """plot data as curves at given position."""
        ip, op = self.selectData( xpos )
        fr = self.freq()
        if ax is None:
            if clf: plt.clf()
            plt.subplot(1,nv,nv-1)
        else:
            plt.sca( ax[0] )
            
        plt.semilogy( ip, fr, marker, label='obs'+addlabel )
        plt.axis('tight')
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
        
        plt.semilogy( op, fr, marker, label='obs'+addlabel )
        if response is not None:
            rop = np.asarray( response )[len(ip):]
            plt.semilogy( rop, fr, rmarker, label='syn'+addlabel )
        
        plt.axis('tight')
        plt.grid(True)
        plt.xlabel('outphase [%]')
        plt.ylabel('f [Hz]')
        plt.legend( loc='best' )
        plt.subplot( 1, nv, 1 )
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
        
    def plotAllData( self, allF = True, orientation='vertical' ):
        """plot data along a profile as image plots for IP and OP."""
        nt = list(range( 0, len( self.x ), 5))
        freq = self.freq()
        nf = len( freq )
        plt.clf()
        ax1 = plt.subplot(211)
        plt.imshow( self.IP[:,self.activeFreq].T, interpolation='nearest' )
        plt.imshow( self.IP[:,self.activeFreq].T, interpolation='nearest' )
        ax1.set_xticks(nt)
        ax1.set_xticklabels(["%g" % xi for xi in self.x[nt]])
        ax1.set_yticks(list(range(0,nf+1,2)))
        ax1.set_yticklabels(["%g" % freq[i] for i in range(0,nf,2)])
        plt.ylim((-0.5,nf-0.5))
        plt.colorbar(orientation=orientation,aspect=30,shrink=0.8)
        plt.xlabel('x [m]')
        plt.ylabel('f [Hz]')
        plt.title('inphase percent')
        ax2 = plt.subplot(212)
        plt.imshow( self.OP[:,self.activeFreq].T, interpolation='nearest' )
        ax2.set_xticks(nt)
        ax2.set_xticklabels(["%g" % xi for xi in self.x[nt]])
        ax2.set_yticks(list(range(0,nf+1,2)))
        ax2.set_yticklabels(["%g" % freq[i] for i in range(0,nf,2)])
        plt.ylim((-0.5,nf-0.5))
        plt.colorbar(orientation=orientation,aspect=30,shrink=0.8)
        plt.xlabel('x [m]')
        plt.ylabel('f [Hz]')
        plt.title('outphase percent')
        return

#    def FOP2d( nlay ):
#        """ 2d forward modelling operator """
#        return FDEM1dModelling( nlay, asvector( self.freq() ), self.cs, -self.height )
#    
        
# Data file example        
#/ MAXMIN ELECTROMAGNETIC SURVEY
#/ FILENAME: Muell150.xyz
#/ PROJECT NUMBER:         1
#/ OPERATOR NUMBER:        3
#/ MAXMIN EQUIPMENT: MMI-9 S/N 3395
#/ SLOPE METHOD: NO SLOPES
#/ COIL SEPARATION:    150.0 METRES 
#/ STATION SPACING:     50.0 METRES 
#/ MODE: MAX1 (Horizontal Coplanar)
#/ FREQUENCIES ON SPREE150.dat FILE: 110, 220, 440, 880, 1760, 3520, 7040, 14080 Hz
#LINE      1.0 Surveying:  N    ELEV        110             220             440             880            1760            3520            7040           14080        110C   220C   440C   880C  1760C  3520C  7040C 14080C    BFC  ERROR
#	0	0	6.61	7.97	8.07	12.42	14.14	19.5	27.66	28.03	45.82	23.67	63.08	11.45	68.98	-8.62	58.82	-20.77
#	20	0	5.04	5.56	7.11	10.31	13.22	16.28	25.06	21.91	37.18	14.17	57.3	4.67	52.07	-17.81	42.18	-31.07

if __name__ == '__main__':
    print("print do some tests here")
