#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as g
import pylab as P
from pygimli.utils.base import draw1dmodel

nlay = 4
lam = 200.
errPerc = 3.

abmnr=g.RMatrix()
g.loadMatrixCol(abmnr,"sond1-100.ves")
ab2  = abmnr[0]
mn2  = abmnr[1]
rhoa = abmnr[2]

transRho  = g.RTransLogLU(1.,1000.)
transThk  = g.RTransLog()
transRhoa = g.RTransLog()

f=g.DC1dModelling( nlay, ab2, mn2 )
f.region( 0 ).setTransModel( transThk )
f.region( 1 ).setTransModel( transRho )

paraDepth = max( ab2 ) / 3
f.region( 0 ).setStartValue( max(ab2) / 3. / nlay / 2. )
f.region( 1 ).setStartValue( P.median( rhoa ) )

model = f.createStartVector()
model[ nlay ] *= 1.5

inv = g.RInversion( rhoa, f, True )
inv.setModel( model )
inv.setTransData( transRhoa )
inv.setRelativeError( errPerc / 100.0 )
inv.setLambda( lam )
inv.setMarquardtScheme( 0.9 )
model = inv.run()

fig = P.figure(1)
fig.clf()
ax1=fig.add_subplot(121)
P.loglog( rhoa, ab2, 'rx-', inv.response(), ab2, 'b-' )
P.axis('tight')
P.ylim((max(ab2),min(ab2)))
P.grid(which='both')
P.xlabel(r"\rho_a in \Omegam")
P.ylabel("AB/2 in m");
P.legend(("measured", "fitted"),loc="upper left")
ax2=fig.add_subplot(122)
res,thk = model(nlay-1,nlay*2-1),model(0,nlay-1)
draw1dmodel(res,thk,r'\rho in \Omega m')
draw1dmodel([100.,500.,20.,1000.],[0.5,3.5,6.])
P.grid(which='both')
P.show()