# -*- coding: utf-8 -*-

from fdem import *

PluginApplication = FDEMApp

MainMenuBarNew_Item = '&FDEM\tCtrl+F'
MainMenuBarNew_ItemHelp = 'Create new FDEM'

MainOpenFileSuffix = [ '.fdem' ]
MainOpenFileSlot = FDEMApp.openFile
MainOpenWildcard = [ "fdem file" ]