import os
from os import system

import pygimli as g

class Rectangle():
    """
    Simple Rectangle that can be written to a xy file.

    nothing else.
    """
    def __init__(self, start, size):
        self.start = start
        self.size = size
        self.points = []

        self.points.append([ self.start[ 0 ] - self.size[ 0 ] /2.0, self.start[ 1 ] - self.size[ 1 ] /2.0 ])
        self.points.append([ self.start[ 0 ] + self.size[ 0 ] /2.0, self.start[ 1 ] - self.size[ 1 ] /2.0 ])
        self.points.append([ self.start[ 0 ] + self.size[ 0 ] /2.0, self.start[ 1 ] + self.size[ 1 ] /2.0 ])
        self.points.append([ self.start[ 0 ] - self.size[ 0 ] /2.0, self.start[ 1 ] + self.size[ 1 ] /2.0 ])

    def writeXY(self, filename, close = False):
        fi = open(filename, 'w')
        for p in self.points:
            fi.write(str(p[ 0 ]) + "\t" + str(p[ 1 ]) + "\n")
        if close:
            p = self.points[ 0 ]
            fi.write(str(p[ 0 ]) + "\t" + str(p[ 1 ]) + "\n")
        fi.close()
        
    def area(self) : return self.size[ 0 ] * self.size[ 1 ]
    # def writeXY
# class Rectangle


def tetgen(filename, quality = 1.2, preserveBoundary = False, verbose = False):
    
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
    mesh = g.Mesh(filebody)
    return mesh
        

def polyAddVIP(filename, pos, marker=0, isRegionMarker=False,
               isHoleMarker=False, maxCellSize=0, verbose=False):
    '''
        out of core wrapper for dcfemlib::polytools::polyAddVIP
    '''
    syscal = "polyAddVIP -x " + str(pos[ 0 ]) + \
                        " -y " + str(pos[ 1 ]) + \
                        " -z " + str(pos[ 2 ])
      
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

def polyAddRectangle(filename, rect, marker = 0, depth = 0, clean = True):
    '''
        out of core wrapper for dcfemlib::polytools 
        merge a meshed horizontal Rectangle with given marker[0] to a 3D PLC at a given depth [0]
        clean removes all out of core files
    '''
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
    '''
        out of core wrapper for dcfemlib::polytools::polyCreateWorld
    '''
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

# def polyCreateWorld

def polyTranslate(filename, x=0.0, y=0.0, z=0.0, verbose=True):
    """
        out of core wrapper for dcfemlib::polytools 
        Translation the PLC in filename by x,y and z
    """
    system("polyTranslate " + 
           " -x " + str(x) +
           " -y " + str(y) +
           " -z " + str(z) + " " + filename)
    
#def polyTranslate(...)

