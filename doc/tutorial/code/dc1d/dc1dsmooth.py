#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as g
import pylab as P
from pygimli.utils.base import draw1dmodel

nlay = 30
lam = 20.
errPerc = 3.
isBlocky = False

abmnr=g.RMatrix()
g.loadMatrixCol(abmnr,"sond1-100.ves")
ab2  = abmnr[0]
mn2  = abmnr[1]
rhoa = abmnr[2]

maxDep = max( ab2 ) / 2.
print "Maximum depth estimated to ", maxDep
thk = g.RVector( nlay - 1, maxDep / ( nlay - 1 ) )
thk *= ( maxDep / sum( thk ) / 3 )

transRho  = g.RTransLog()
transRhoa = g.RTransLog()

f=g.DC1dRhoModelling( thk, ab2, mn2 )

inv = g.RInversion( rhoa, f, True )
model = g.RVector( nlay, P.median( rhoa ) )
inv.setModel( model )
inv.setTransData( transRhoa )
inv.setTransModel( transRho )
inv.setRelativeError( errPerc / 100.0 )
inv.setLambda( lam )
model = inv.run()

model2 = g.RVector( nlay, P.median( rhoa ) )
inv.setModel( model2 )
inv.setBlockyModel( True )
model2 = inv.run()

fig = P.figure(1)
fig.clf()
ax1=fig.add_subplot(121)
P.loglog( rhoa, ab2, 'rx-', inv.response(), ab2, 'b-' )
P.axis('tight')
P.ylim((max(ab2),min(ab2)))
P.grid(which='both')
P.xlabel(r"\rho_a in \Omegam");
P.ylabel("AB/2 in m");
P.legend(("measured", "fitted"),loc="upper left");
ax2=fig.add_subplot(122)
draw1dmodel(model,thk,r'\rho in \Omega m')
draw1dmodel(model2,thk,r'\rho in \Omega m')
draw1dmodel([100, 500, 20, 1000],[0.5,3.5,6])
P.show()