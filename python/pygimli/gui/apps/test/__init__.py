# -*- coding: utf-8 -*-

import sys
import traceback

try:
    from . testapp import *
except Exception as e:
    traceback.print_exc(file=sys.stdout)
    sys.stderr.write(str(e))
    sys.stderr.write("Fail to __init__ testapp.")

PluginApplication = TestApp

MainMenuBarNew_Item = '&Test\tCtrl+T'
MainMenuBarNew_ItemHelp = 'Create new Test application'

MainOpenFileSuffix = ['.tt', '.ttt']
MainOpenFileSlot = TestApp.openFile
MainOpenWildcard = ["TestApp test suffix1 (*.tt)",
                    "TestApp test suffix2 (*.ttt)"]
