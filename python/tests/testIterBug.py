#import sys, os
# sys.path.append(os.path.abspath('../pygimli'))
#import _pygimli_ as pg

import pygimli as pg

a = pg.RVector(10, 1)

# das geht schief wegen fehlendem referenzcounter. der Iter nutzt das
# temporäre Object a(0,9) das nicht weiter gezählt wird
print(a(0, 9).beginPyIter()[0], "!=", a[0])

# das geht weil wir einen eigenen iter definieren der die referenz hällt
for ai in a(0, 9):
    print(ai)
