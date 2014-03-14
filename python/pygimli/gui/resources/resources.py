# -*- coding: iso-8859-1 -*-

import sys, os
import wx

def iconFileName( filename ):
    globPath = os.path.dirname( __file__ )
    iconpath = os.path.join( globPath, "icons" ) 
    
    if hasattr( sys, "frozen"):
        globPath = os.path.dirname( sys.argv[ 0 ] )
    respath  = os.path.join( globPath, "resources" )
    iconpath = os.path.join( respath, "icons" )
    
    return os.path.join( iconpath, filename )
        
def loadIcon( filename ):
    """Load a bitmap file from the resource/icons subdirectory Returns a
    wx.Bitmap object."""

    #globPath = os.path.dirname( __file__ )
    #iconpath = os.path.join( globPath, "icons" )
    #iconfile = os.path.join( iconPath(), filename )

    #if hasattr( sys, "frozen"):
        #globPath = os.path.dirname( sys.argv[ 0 ] )
    #respath  = os.path.join( globPath, "resources" )
    #iconpath = os.path.join( respath, "icons" )
    #iconfile = os.path.join( iconpath, filename )
    #pass
    #else:
        #pass
    iconfile = iconFileName( filename )
    if not os.path.exists( iconfile ):
        raise IOError('Could not find icon file "%s"; dying'%iconfile)

    bmp = wx.Bitmap( iconfile )
    return bmp

def MakeGray(xxx_todo_changeme, factor, maskColor):
    """
    Make a pixel grayed-out. If the pixel matches the `maskColor`, it won't be
    changed.

    :param `(r,g,b)`: a tuple representing a pixel colour;
    :param `factor`: a graying-out factor;
    :param `maskColor`: a colour mask.
    """
    (r,g,b) = xxx_todo_changeme
    if (r,g,b) != maskColor:
        return [int((230 - x) * factor) + x for x in (r,g,b)]
    else:
        return (r,g,b)
        

def MakeDisabledBitmap(bitmap):
    """
    Convert the given image (in place) to a grayed-out version,
    appropriate for a 'disabled' appearance., 
    
    Taken from aui_utilities and add alpha 

    :param: `bitmap`: the bitmap to gray-out.
    """
    anImage = bitmap.ConvertToImage()    
    factor = 0.7        # 0 < f < 1.  Higher Is Grayer
    
    if anImage.HasMask():
        maskColor = (anImage.GetMaskRed(), anImage.GetMaskGreen(), anImage.GetMaskBlue())
    else:
        maskColor = None
        
    if anImage.HasAlpha():
        alpha = anImage.GetAlphaData()
    else:
        alpha = None

    data = list(map(ord, list(anImage.GetData())))

    for i in range(0, len(data), 3):
        
        pixel = (data[i], data[i+1], data[i+2])
        pixel = MakeGray(pixel, factor, maskColor)

        for x in range(3):
            data[i+x] = pixel[x]

    anImage.SetData(''.join(map(chr, data)))
    
    if alpha:
        anImage.SetAlphaData(alpha)
        
    return anImage.ConvertToBitmap()

def loadXRC( filename ):
    """
    Load a XRC-file from the resource/xrc subdirectory
    Returns a xml object
    """
    xml = wx.xrc.EmptyXmlResource()
    xrcfile = filename + '.xrc'
    
    if hasattr( sys, "frozen"):
        globPath = os.path.dirname( sys.argv[ 0 ] )
        respath  = os.path.join( globPath, "resources" )
        xrcpath = os.path.join( respath, "xrc" )
    else:
        globPath = os.path.dirname( __file__ )
        xrcpath = os.path.join( globPath, "xrc" )

    xrcfile = os.path.join( xrcpath, xrcfile )
        
    if not os.path.exists( xrcfile ):
        raise IOError('Could not find xrc file "%s"; dying' % xrcfile )

    xml.Load( xrcfile )

    return xml
