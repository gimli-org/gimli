# -*- coding: utf-8 -*-

from xml.dom.minidom import parse

def handleWPTS( wpts ):
    """ Handler for Waypoints in gpx xml-dom """
    w = []

    for wpt in wpts:
        if wpt.hasAttribute('lat'):
            lat = float( wpt.getAttribute('lat') )
        else:
            continue
        if wpt.hasAttribute('lon'):
            lon = float( wpt.getAttribute('lon') )
        else:
            continue

        

        name = wpt.getElementsByTagName( 'name' )[0].childNodes[0].data
        time = wpt.getElementsByTagName( 'time' )[0].childNodes[0].data

        w.append( ( lon, lat, name, time)  )
    return w
#def findWPTS( ... )

def readGPX( filename ):
    """ Extract GPS Waypoint from GPS Exchange Format (GPX). Currently only simple waypoint extraction is supported. """
    
    dom = parse( filename )
    wpts = dom.getElementsByTagName("wpt")

    return handleWPTS( wpts )
# def readGPX( ... )


def readSimpleLatLon( filename ):
    """ """
    w = []
    
    with open( filename, 'r') as fi:
        content = fi.readlines( )
    fi.close()
    
    for line in content:
        vals = line.split()
        if len( vals ) > 1:
            w.append( (float(vals[1]), float(vals[0]), '', 'time')  )
            print w
        #w.append( (float(vals[2]), float(vals[1]), vals[0], 'time')  )
    
    return w
    