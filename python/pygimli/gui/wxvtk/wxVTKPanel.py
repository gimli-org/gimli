# -*- coding: iso-8859-1 -*-
#use local copy of wxVTKRenderWindowInteractor until the baseClass.__init__() order is adapted to glcanvas-2.9

import wx

vtk = None

try:
    import vtk
    from wxVTKRenderWindowInteractor import *
except Exception, e:
    print e
    print "no wxVTK installation found"
    print "for win32 you can install from http://cpbotha.net/software/latest-vtk-windows-binaries/"
    wxVTKRenderWindowInteractor = wx.Panel

try:
    import getHandleHack as wxHandle 
except Exception, e:
    pass

class wxVTKPanel( wxVTKRenderWindowInteractor ):
    def __init__( self, parent, ID, *args, **kwargs ):
        wxVTKRenderWindowInteractor.__init__( self, parent, ID, *args, **kwargs )
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground( 0.8, 0.8, 0.8 )
        self.GetRenderWindow().AddRenderer( self.renderer ) 

    def GetHandle( self ):
        '''
            inbuild GetHandle function fails on 2.9.0
        '''
        handle = baseClass.GetHandle( self )
        if not handle:
            return wxHandle.get_window_handle_str( self )
        else:
            return handle

from pygimli.gui.base import AppResource 

class AppResourceWxVTK( AppResource, wxVTKPanel ):
    def __init__( self, parent, rendererSlot, propertyInspectorSlot ):
        AppResource.__init__( self, parent, rendererSlot, propertyInspectorSlot )
        wxVTKPanel.__init__( self, parent, -1 )
        
    def getRendererPanel( self ) : return self
