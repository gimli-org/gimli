import numpy as np
from pygimli import FDEM1dModelling, RVector, asvector, RTrans, RTransLog, RTransLogLU, RInversion
from pygimli.utils import draw1dmodel
from math import sqrt, pi

def rhoafromU( UbyI, t, Tx, Rx=None ):
    ''' apparent resistivity curve from classical TEM '''
    if Rx is None:
        Rx = Tx # assume single/coincident loop
    
    mu0 = 4e-7 * pi
#    rhoa = ( Rx * Tx * mu0 / 20. / UbyI )**(2./3.) * t**(-5./3.) * mu0 / pi
    rhoa = ( Rx * Tx * mu0 / 20. / UbyI )**(2./3.) * t**(-5./3.) * mu0 / pi
    return rhoa

def rhoafromB( B, t, Tx, I=1 ):
    ''' apparent resistivity from B-field TEM '''
    mu0 = 4e-7 * pi
    rhoa = ( I * Tx * mu0 / 30. / B )**(2./3.) * mu0 / pi / t
    return rhoa

def rhoa( snd, cal=260e-9, corrramp=True ):
    ''' compute apparent resistivity from sounding '''
    Tx = np.prod( [float(a) for a in snd['LOOP_SIZE'].split()] )
    if 'COIL_SIZE' in snd:
        Rx = snd['COIL_SIZE']
    else:
        Rx = Tx
    
    v = snd['VOLTAGE']
    istart, istop = 0, len(v) # default: take all
    mav = np.find( v == max(v) )
    if len(mav) > 1: #several equal big ones: start after
        istart = max(mav)+1

    if min(v) < 0.0: # negative values: stop at first
        istop = min( np.find( v < 0.0 ) )
    
    v = v[istart:istop]
    dv = snd['ST_DEV'][istart:istop] #/ snd['CURRENT']
    t = snd['TIME'][istart:istop]
    if corrramp and 'RAMP_TIME' in snd: 
        t = t - snd['RAMP_TIME']
    
    if Rx == 1: # apparently B-field
        rhoa = rhoafromB( v * cal, t, Tx )
    else:
        rhoa = rhoafromU( v, t, Tx, Rx )

    rhoaerr = dv / v * (2./3.)
    return rhoa, t, rhoaerr

def readusffile( filename ):
    ''' read data from single USF (universal sounding file) file
        DATA = readusffile( filename )
        DATA = readusffile( filename, DATA ) will append to DATA '''
    
    DATA = []
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
    ALLDATA = []
    for onefile in filenames:
        ALLDATA.extend( readusffile( onefile ) )
    
    return ALLDATA

def readSiroTEMData( fname ):
    ''' read TEM data from siroTEM instrument dump 
        DATA = readSiroTEMData( filename )
        .. list of soundings with USF and siro-specific keys 
    '''
    Time_ST = np.array( [487.,887.,1287.,1687.,2087.,2687.,3487.,4287.,5087.,5887.,7087.,8687.,10287.,11887.,13487.,
                    15887.,19087.,22287.,25487.,28687.,33487.,39887.,46287.,52687.,59087.,68687.,81487.,94287.,
                    107090.,119890.,139090.,164690.,190290.,215890.,241490.,279890.,331090.,382290.,433490.,
                    484690.,561490.,663890.,766290.,868690.,971090.,1124700.,1329500.,1534300.,1739100.,1943900.] )
    Time_ET = np.array( [0.05, 0.1, 0.15, 0.25, 0.325, 0.425, 0.525, 0.625, 0.725, 0.875, 1.075, 1.275, 1.475, 1.675,
                        1.975, 2.375, 2.775, 3.175, 3.575, 4.175, 4.975, 5.775, 6.575, 7.375, 8.575, 10.175, 11.775,
                        13.375, 14.975, 17.375, 20.575, 23.775, 26.975, 30.175, 34.975, 41.375, 47.775, 54.175, 60.574,
                        70.175, 82.975, 95.775,108.575,121.375,140.575,166.175,191.775,217.375,242.975,281.375,332.575] )
    
    fid = open( fname )
    # read in file header until : sign
    line = 'a'
    while len( line ) > 0 and line[0] != ':':
        line = fid.readline()

    DATA = []
    line = fid.readline()
    while line[0] != ';':
        header = line[1:-6].split(',')

        snd = {} # dictionary, uppercase corresponds to USF format keys
        snd['INSRTUMENT'] = 'siroTEM'
        snd['dtype'] = int(header[3])
        dstring = header[1]
        snd['DATE'] = int( '20' + dstring[6:8] + dstring[3:4] + dstring[0:1] )
        snd['win0'], snd['win1'], ngain, snd['conf'], snd['nch'], snd['SOUNDING_NUMBER']  = \
            [int(h) for h in header[5:11]]
        
        snd['GAIN_FACTOR'] = [ 0.1, 1.0, 10.0, 100.0 ][ ngain ] # predefined gain factors
        snd['STACK_SIZE'] = int( header[14] )
        snd['ttype'] = int( header[20] ) # 1-composite,2-earlytime,3-standard,4-highresolution
        snd['CURRENT'], snd['RAMP_TIME'] = float( header[17] ), float( header[18] )*1e-6
        snd['TIME_DELAY'] = float( header[19] )
        snd['LOOP_SIZE'], snd['COIL_SIZE'] = float( header[21]), float( header[22])

        dummy = fid.readline()
        data = []
        line = fid.readline()[:-1] # trim CR+LF newline
        while len(line) >0:
            while line[-1] == '/':
                line = line[:-1] + fid.readline()[:-1].replace('\t','')
                aline = line

            nums = [float(el[-7:-2])*10**(float(el[-2:])) for el in line[1:-5].split(',')[1:]]
            data.append( np.array( nums ) )
            line = fid.readline()[:-1] # trim newline

        snd['VOLTAGE'] = data[0]
        if snd['ttype'] == 2: # early time
            snd['TIME'] = Time_ET[ snd['win0']-1:snd['win1'] ] * 1e-3
        if snd['ttype'] == 3: # standard time
            snd['TIME'] = Time_ST[ snd['win0']-1:snd['win1'] ] * 1e-6
        
        snd['ST_DEV'] = data[1]
        if snd['dtype'] >0: #normal measurement
            DATA.append( snd )

        line = fid.readline()

    fid.close() 
    return DATA

class TDEMData():
    """ TEM class holding data etc. """
    def __init__( self, filename, height=1.0, verbose=False ):
        """ initialize data class and load data """

    def __init__( self, filename = None):
        """ initialize data class and load data """
        self.x, self.f, self.cs, self.IP, self.OP, self.ERR = None,None,None,None,None,None
        self.height = 1.0
        
        if filename:
            # check if filename extension is usf
            if filename.lower().rfind( '.usf' ) > 0:
                if filename.find( '*' ) > 0:
                    self.DATA = readusffiles( filename )
                else:
                    self.DATA = readusffile( filename )
            elif filename.lower().rfind( '.txt' ) > 0:
                self.DATA = readSiroTEMData( filename )
            
        self.activeFreq = ( self.f > 0.0 )

    def __repr__( self ):
        if isinstance( self.x, float ):
            return "<FDEMdata: sounding with %d frequencies, coilspacing is %.1f>" % (len(self.f), self.cs)
        else:
            return  "<FDEMdata: %d soundings with each %d frequencies, coilspacing is %.1f>" % (len(self.x), len(self.f), self.cs)
        
    def showInfos( self ): # only for old scripts using it
        print(__repr__( self ))


    
if __name__ == '__main__':
    print("print do some tests here")
