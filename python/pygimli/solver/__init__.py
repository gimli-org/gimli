# -*- coding: utf-8 -*-

from .solver import *
from .green import *

# unsorted stuff need

def divergence(mesh, F=lambda r:r, order=1):
    """"""
    div = 0
    directionCheck = False
    
    if mesh.cellCount() > 0:
        directionCheck = True
        
    for b in mesh.boundaries():
               
        if directionCheck:
            if b.leftCell() is None and b.rightCell() is None:
                print(b.id(), b.leftCell(), b.rightCell())  
                sw = g.Stopwatch(True)
                mesh.createNeighbourInfos()
                print("NeighbourInfos()", sw.duration(True))
                ##return gauss(grid, F)
              
            if not b.leftCell() is None and not b.rightCell() is None: continue
               
        tmpdiv = 0
        shape = b.shape()
        
        if order == 1:
            tmpdiv = shape.norm().dot(F(shape.center())) * shape.domainSize()
        else:
            weights = g.IntegrationRules.instance().weights(shape, order)
            abscissa = g.IntegrationRules.instance().abscissa(shape, order)
                    
            for i, p in enumerate(abscissa): 
                rPos = shape.xyz(p)
                tmpdiv += shape.norm().dot(F(rPos)) * weights[i] * shape.domainSize()
       
        
        if directionCheck and b.leftCell() is None:
            tmpdiv *= -1
            #raise Exception("invalid mesh: left is None .. every boundary need leftCell")
                    
        div += tmpdiv
        
    return div
# def divergence