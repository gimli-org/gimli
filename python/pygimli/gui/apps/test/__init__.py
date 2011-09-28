# -*- coding: iso-8859-1 -*-

from testapp import *

PluginApplication = TestApp

MainMenuBarNew_Item = '&Test\tCtrl+T'
MainMenuBarNew_ItemHelp = 'Create new Test application'

MainOpenFileSuffix = [ '.tt', '.ttt' ]
MainOpenFileSlot = TestApp.openFile
MainOpenWildcard = [ "TestApp test suffix1 (*.tt)", "TestApp test suffix2 (*.ttt)" ]