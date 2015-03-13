# -*- coding: utf-8 -*-
"""
    provide matplotlib qt widget
"""

import sys
import traceback

try:
    from . qtMatplotPanel import MatplotPanel
    from . qtMatplotPanel import AppResourceMPL
except Exception as e:
    traceback.print_exc(file=sys.stdout)
    sys.stderr.write(str(e))
