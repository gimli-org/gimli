# -*- coding: utf-8 -*-

import pygimli as pg
#import pygimli.utils

import numpy as np
#from numpy import arange, ndarray, array, ma
#import matplotlib as mpl

#from colorbar import *

#def annotateSeparationAxis(ax, shemeID, grid = False):
    #'''
        #Draw y-axes tick labels corresponding to the separation
    #'''
    #prefix = DataShemeManager().sheme(shemeID).prefix
    
    #def sepName(l):
        #suffix = ""
        
        #if l == 0:
            #return ''
        #elif l > 0:
            #suffix = "'"
                
        #if grid:
            #ax.plot(ax.get_xlim(), [l,l], color = 'black', linewidth = 1, linestyle='dotted')
        
        #return prefix + ' $' + str(abs(int(l))) + suffix +'$'

    #ax.yaxis.set_ticklabels(map(lambda l: sepName(l), ax.yaxis.get_ticklocs()))
    
## END def annotateSeparations(...)

def drawSensorAsMarker(ax, data):
    '''
        Draw Sensor marker, these marker are pickable
    '''
    elecsX = []
    elecsY = []
    
    for i in range(len(data.sensorPositions())):
        elecsX.append(data.sensorPositions()[i][ 0 ])
        elecsY.append(data.sensorPositions()[i][ 1 ])
    
    electrodeMarker, =  ax.plot(elecsX, elecsY, 'x', color = 'black', picker = 5.)
    
    ax.set_xlim([ data.sensorPositions()[0][0]-1., data.sensorPositions()[ data.sensorCount() -1][0] + 1. ])
    #print electrodeMarker
    return electrodeMarker
# END def drawElectrodesAsMarker(...)


def drawTravelTimeData(a, data):
    '''
        Draw first arrival traveltime data into mpl axes a. 
        data of type \ref DataContainer must contain sensorIdx 's' and 'g' and thus numbered internal from [0..n)
    '''

    x = pg.x(data.sensorPositions())
    z = pg.z(data.sensorPositions())

    shots = pg.unique(pg.sort(data('s')))
    geoph = pg.unique(pg.sort(data('g')))

    startOffsetIDX = 0

    if min(min(shots), min(geoph) == 1): 
        startOffsetIDX = 1
    
    a.set_xlim([ min(x), max(x) ])
    a.set_ylim([ max(data('t')), -0.002 ])
    a.figure.show()
    
    for shot in shots:
        gIdx = pg.find(data('s') == shot)
        sensorIdx = [int(i__ - startOffsetIDX) for i__ in data('g')[ gIdx ]]
        a.plot(x[ sensorIdx ], data('t')[ gIdx ], 'x-')
            
    yPixel = a.transData.inverted().transform_point((1, 1))[1]-a.transData.inverted().transform_point((0, 0))[1]
    xPixel = a.transData.inverted().transform_point((1, 1))[0]-a.transData.inverted().transform_point((0, 0))[0]

    # draw shot points
    a.plot(x[ [int(i__ - startOffsetIDX) for i__ in shots] ], np.zeros(len(shots)) + 8.*yPixel, 'gv', markersize = 8)    

    # draw geophone points
    a.plot(x[ [int(i__ - startOffsetIDX) for i__ in geoph] ], np.zeros(len(geoph)) + 3.*yPixel, 'r^', markersize = 8)    

    a.grid()    
    a.set_ylim([ max(data('t')), +16.*yPixel])
    a.set_xlim([ min(x)-5.*xPixel, max(x)+5.*xPixel ])

    a.set_xlabel('x-Coordinate [m]')
    a.set_ylabel('Traveltime [ms]')
# def drawTravelTimeData(...)