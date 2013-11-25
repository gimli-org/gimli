# -*- coding: utf-8 -*-
"""provide matplotlib wxGTK panel."""
import sys, traceback

try:
	from wxMatplotPanel import wxMatplotPanel, AppResourceWxMPL
	from colorBarWxMPL import ColorBarWxMPL
except Exception:
	print '-'*60
	traceback.print_exc( file = sys.stdout)
	print '-'*60
