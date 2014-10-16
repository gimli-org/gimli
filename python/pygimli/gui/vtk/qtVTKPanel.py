# -*- coding: utf-8 -*-
"""
    Use local copy of wxVTKRenderWindowInteractor until the baseClass.__init__() order is adapted to glcanvas-2.9
"""
import sys
import os
import traceback
vtk = None

try:
    import vtk
except ImportError as e:
    traceback.print_exc(file=sys.stdout)
    sys.stderr.write("no VTK installation found")
    sys.stderr.write("for win32 you can install from http://cpbotha.net/software/latest-vtk-windows-binaries/")

#class VTKPanel(wxVTKRenderWindowInteractor):
    #def __init__(self, parent):
        #wxVTKRenderWindowInteractor.__init__(self, parent)
        #self.renderer = vtk.vtkRenderer()
        #self.renderer.SetBackground(0.8, 0.8, 0.8)
        #self.GetRenderWindow().AddRenderer(self.renderer) 

from pygimli.gui.base import AppResource 

class AppResourceVTK(AppResource, vtk.QVTKWidget):
    def __init__(self, parent, rendererSlot, propertyInspectorSlot):
        AppResource.__init__(self, parent, rendererSlot, propertyInspectorSlot)
        vtk.QVTKWidget.__init__(self, parent)
        
    def getRendererPanel(self): 
        """
        """
        return self
