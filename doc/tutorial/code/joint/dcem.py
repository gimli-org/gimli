#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as g
import pylab as P
from pygimli.utils.base import draw1dmodel

class DCEM1dModelling( g.ModellingBase ):
    """
        modelling jointing DC and EM 1Dforward operators
    """
    def __init__( self, nlay, ab2, mn2, freq, coilspacing, verbose = False  ):
        """
            constructor, number of layers, AB/2 & MN/2, frequencies & coil spacing */
        """
        g.ModellingBase.__init__( self, verbose )
        self.nlay_ = nlay
        self.fDC_ = g.DC1dModelling( nlay, ab2, mn2, verbose )
        self.fEM_ = g.FDEM1dModelling( nlay, freq, coilspacing, verbose )
        self.mesh_ = g.createMesh1DBlock( nlay )
        self.setMesh( self.mesh_ )

    def response( self, model ):
        return g.cat( self.fDC_( model ), self.fEM_( model ) )


# options to play with
noiseEM = 1. # absolute (per cent of primary signal)
noiseDC = 3. # in per cent
lamEM, lamDC, lamDCEM = 300., 500., 500.
verbose = False

# synthetic model 
nlay = 3
model = g.RVector( nlay * 2 - 1 )
for i in range( nlay - 1 ): model[ i ] = 15.
for i in range( nlay ): model[ i + nlay -1 ] = 200.
model[ nlay ] = 10.
model[ nlay + 1 ] = 50.
thk = model( 0, nlay -1 )
res = model( nlay - 1, nlay * 2 - 1 )

# EM forward operator and synthetic data
coilspacing = 50.
nf = 10
freq = g.RVector( nf, 110. )
for i in range( nf - 1 ): freq[ i +1 ] = freq[ i ] * 2.
fEM = g.FDEM1dModelling( nlay, freq, coilspacing )
dataEM = fEM(model)
for i in range(len(dataEM)): dataEM[i] += P.randn(1)[0] * noiseEM

# model transformations
transRhoa = g.RTransLog()
transThk = g.RTransLog()
transRes = g.RTransLogLU(1., 1000.)
transEM = g.RTrans()
fEM.region(0).setTransModel( transThk )
fEM.region(1).setTransModel( transRes )

# independent EM inversion
invEM = g.RInversion( dataEM, fEM, transEM, verbose )
modelEM = g.RVector( nlay * 2 - 1, 50. )
invEM.setModel( modelEM )
invEM.setAbsoluteError( noiseEM )
invEM.setLambda( lamEM )
invEM.setMarquardtScheme( 0.9 )
modelEM = invEM.run()
respEM = invEM.response()

# DC forward operator and synthetic data
ab2 = g.RVector( 20, 3. )
na = len( ab2 )
mn2 = g.RVector( na, 1.0 )
for i in range( na - 1 ): ab2[ i +1 ] = ab2[ i ] * 1.3
fDC = g.DC1dModelling( nlay, ab2, mn2 )
dataDC = fDC(model)
for i in range(len(dataDC)): dataDC[i] *= 1. + P.randn(1)[0] * noiseDC / 100.
fDC.region(0).setTransModel( transThk )
fDC.region(1).setTransModel( transRes )

# independent DC inversion
invDC = g.RInversion( dataDC, fDC, transRhoa, verbose )
modelDC = g.RVector( nlay * 2 - 1, 20. )
invDC.setModel( modelDC )
invDC.setRelativeError( noiseDC / 100. )
invDC.setLambda( lamDC )
invDC.setMarquardtScheme( 0.9 )
modelDC = invDC.run()
respDC = invDC.response()

# joint forward operator
fDCEM = DCEM1dModelling( nlay, ab2, mn2, freq, coilspacing )
fDCEM.region(0).setTransModel( transThk )
fDCEM.region(1).setTransModel( transRes )

# joint inversion
transData = g.RTransCumulative()
transData.push_back( transRhoa, na )
transData.push_back( transEM, nf*2 )
invDCEM = g.RInversion( g.cat( dataDC, dataEM ), fDCEM, transData, verbose )
modelDCEM = g.RVector( nlay * 2 - 1, 20. )
invDCEM.setModel( modelDCEM )
err = g.cat( dataDC * noiseDC / 100., g.RVector( len(dataEM), noiseEM) )
invDCEM.setAbsoluteError( err )
invDCEM.setLambda( lamDCEM )
invDCEM.setMarquardtScheme( 0.9 )
modelDCEM = invDCEM.run()
respDCEM = invDCEM.response()

# comparison
invDC.echoStatus()
invEM.echoStatus()
invDCEM.echoStatus()
[invDC.chi2(), invEM.chi2(), invDCEM.chi2()]

# plot results
fig=P.figure(1)
fig.clf()
ax1=fig.add_subplot(131)
draw1dmodel( res, thk )
draw1dmodel( modelEM( nlay-1, nlay*2-1 ), modelEM(0,nlay-1) )
draw1dmodel( modelDC( nlay-1, nlay*2-1 ), modelDC(0,nlay-1) )
draw1dmodel( modelDCEM( nlay-1, nlay*2-1 ), modelDCEM(0,nlay-1) )
P.legend(('Syn','EM','DC','JI'))
P.xlim((10.,1000.))
P.ylim((40.,0.))
P.grid(which='both')
ax2=fig.add_subplot(132)
P.semilogy( dataEM( 0, nf ), freq, 'bx', dataEM(nf,nf*2), freq, 'bo')
P.semilogy( respEM( 0, nf ), freq, 'g-', respEM(nf,nf*2), freq, 'g-')
P.semilogy( respDCEM( na, na+nf ), freq, 'c-', respDCEM(na+nf, na+nf*2),freq,'c-')
P.ylim((min(freq),max(freq)))
P.grid(which='both')
P.legend(("syn","","EM","","DCEM",""),loc="lower left")
ax3=fig.add_subplot(133)
P.loglog( dataDC, ab2, 'bx-', respDC, ab2, 'r-', respDCEM(0,na), ab2, 'c-' )
P.axis('tight')
P.ylim( ( max(ab2), min(ab2) ) )
P.grid(which='both')
P.xlabel(r'\rho_a in \Omegam')
P.ylabel("AB/2 in m");
P.legend( ("syn", "DC","EMDC"), loc="lower left" )
P.show()