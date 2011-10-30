# -*- coding: iso-8859-1 -*-

from gpsview import *

PluginApplication = GPSViewerApp

MainMenuBarNew_Item = '&GPS View\tCtrl+G'
MainMenuBarNew_ItemHelp = 'Create new GPS viewer'

MainOpenFileSuffix = [ '.gpx' ]
MainOpenFileSlot = GPSViewerApp.openFile
MainOpenWildcard = [ "GPS Exchange Format (*.gpx)" ]