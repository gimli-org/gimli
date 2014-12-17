#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This program is part of pygimli
Visit http://www.resistivity.net for further information or the latest version.
"""

import sys
from os import system

try:
    import pygimli as pg

except ImportError:
    sys.stderr.write('ERROR: cannot import the library pygimli.'+ \
                     ' Ensure that pygimli is in your PYTHONPATH')
    sys.exit(1)


def main(argv):
    from optparse import OptionParser

    parser = OptionParser("usage: %prog [options] mesh|mod")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
                      help="be verbose", default=False)

    (options, args) = parser.parse_args()

    print((options, args))

if __name__ == "__main__":
    main(sys.argv[1:])
