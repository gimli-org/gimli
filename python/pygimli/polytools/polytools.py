# -*- coding: utf-8 -*-
""" Tools to create or manage PLC """

import os
from os import system

import pygimli as pg

class Rectangle():
    """
    Simple Rectangle that can be written to a xy file.

    Nothing else.
    """
    def __init__(self, start, size):
        """Initialize rectangle with start position and a given size."""
        self.start = start
        self.size = size
        self.points = []

        self.points.append([self.start[0] - self.size[0] /2.0, self.start[1] - self.size[1] /2.0])
        self.points.append([self.start[0] + self.size[0] /2.0, self.start[1] - self.size[1] /2.0])
        self.points.append([self.start[0] + self.size[0] /2.0, self.start[1] + self.size[1] /2.0])
        self.points.append([self.start[0] - self.size[0] /2.0, self.start[1] + self.size[1] /2.0])

    def writeXY(self, filename, close = False):
        """ Write coordinates to a file """
        fi = open(filename, 'w')
        for p in self.points:
            fi.write(str(p[0]) + "\t" + str(p[1]) + "\n")
        if close:
            p = self.points[0]
            fi.write(str(p[0]) + "\t" + str(p[1]) + "\n")
        fi.close()

    def area(self): 
        """ Return size of the Rectangle"""
        return self.size[0] * self.size[1]

def tetgen(filename, quality=1.2, preserveBoundary=False, verbose=False):
    """
    Create a :term:`Tetgen` :cite:`Si2004` mesh from a PLC.
    
    Forwards to system call tetgen, which must be known to your system.
    
    Parameters
    ----------
    filename: str
       
    quality: float [1.2]
        Refines mesh (to improve mesh quality). [1.1 ... ]
        
    preserveBoundary: bool [False]
        Preserve PLC boundary mesh
        
    verbose: bool [False]
        be verbose
        
    Returns
    -------
    mesh: gimliapi:`GIMLI::Mesh`
    """

    filebody = filename.replace('.poly', '')
    syscal = 'tetgen -pazAC'
    syscal += 'q' + str(quality)

    if not verbose:
        syscal += 'Q'
    else:
        syscal += 'V'

    if preserveBoundary:
        syscal += 'Y'

    syscal += ' ' + filebody + '.poly'

    if verbose:
        print(syscal)

    system(syscal)
    system('meshconvert -it -BD -o ' + filebody + ' ' + filebody + '.1')
    try:
        os.remove(filebody + '.1.node')
        os.remove(filebody + '.1.ele')
        os.remove(filebody + '.1.face')
    except:
        None
    mesh = pg.Mesh(filebody)
    return mesh


def polyAddVIP(filename, pos, marker=0, isRegionMarker=False,
               isHoleMarker=False, maxCellSize=0, verbose=False):
    """
    Add very important point (VIP) to a PLC.
    
    Out of core wrapper for dcfemlib::polytools::polyAddVIP.
    
    Parameters
    ----------
    
    Returns
    -------
    """
    
    syscal = "polyAddVIP -x " + str(pos[0]) + \
                        " -y " + str(pos[1]) + \
                        " -z " + str(pos[2])

    if isHoleMarker:
        syscal += " -H "
    else:
        syscal += " -m " + str(marker)

    if isRegionMarker:
        syscal += " -R "

    if maxCellSize > 0:
        syscal += " -a " + str(maxArea)

    syscal += " " + filename

    if verbose: print(syscal)
    system(syscal)
# def polyAddVIP

def polyAddRectangle(filename, rect, marker=0, depth=0, clean=True):
    """
    Add horizontal plane to a PLC
    
    Out of core wrapper for dcfemlib::polytools::polytools.
    Merge a meshed horizontal Rectangle with given marker[0] to
    a 3D PLC at a given depth [0] clean removes all out of core files
        
    Parameters
    ----------
    
    Returns
    -------
    """
    
    rect.writeXY("__pad.xy", close = True)
    system("polyCreateWorld -d2 -t __pad.xy -C __pad")
    a = rect.area() / 29.0

    system("dctriangle -a " + str(a) + " -q34.0 -S __pad")
    system("polyCreateFacet -o __pad3d -m " + str(marker) + " __pad.bms")
    system("polyTranslate -z " + str(depth) + " __pad3d")

    ### add node to the center of the recangle
#    system("polyAddVIP -x " + str(rect.start[0])
#                    + " -y " + str(rect.start[1])
#                    + " -m -1 __pad3d")
    system("polyMerge " + filename + " __pad3d " + filename)

#    system("polyCreateFacet -o __pad3d -m 1 __pad.bms")
#    system("polyTranslate -z -0.1 __pad3d")
#    system("polyMerge " + filename + " __pad3d " + filename)

    if clean:
        os.remove('__pad.xy')
        os.remove('__pad.poly')
        os.remove('__pad.bms')
        os.remove('__pad3d.poly')
#def polyAddRectangle


def polyCreateWorld(filename, x=None, depth=None, y=None, marker=0,
                    maxCellSize=0, verbose=True):
    """
    Create the PLC of a default world.
    
    Out of core wrapper for dcfemlib::polytools::polyCreateWorld
    
    Parameters
    ----------
    
    Returns
    -------
    """
    
    if depth is None:
        print("Please specify worlds depth.")
        return

    if x is None:
        print("Please specify worlds x dimension.")
        return

    dimension = 3
    z = depth

    if y is None:
        dimension = 2

    syscal = 'polyCreateWorld -d ' + str(dimension) \
                                + ' -x ' + str(x) \
                                + ' -y ' + str(y) \
                                + ' -z ' + str(z) \
                                + ' -m ' + str(marker) \

    if maxCellSize > 0:
        syscal += " -a " + str(maxCellSize)

    syscal = syscal + ' ' + filename

    if verbose:
        print(syscal)

    os.system(syscal)

def polyTranslate(filename, x=0.0, y=0.0, z=0.0, verbose=True):
    """
    Translate (move) a PLC
    
    Out of core wrapper for dcfemlib::polytools.
    Spatial translate (move) the PLC (filename) by x, y and z
        
    Parameters
    ----------
    
    Returns
    -------
    """
    system("polyTranslate " +
           " -x " + str(x) +
           " -y " + str(y) +
           " -z " + str(z) + " " + filename)

