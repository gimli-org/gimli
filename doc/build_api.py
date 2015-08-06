#!/usr/bin/env python
"""Script to generate API docs"""
import os, sys
from apigen import ApiDocWriter

if __name__ == '__main__':
    package = 'pygimli'

    try:
        __import__(package)
    except ImportError as e:
        print("Can not import pygimli")
        exit()

    outdir = 'pygimliapi'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    docwriter = ApiDocWriter(package)
    docwriter.package_skip_patterns += [r'\.gui$']
    docwriter.rst_extension = ".rst"
    docwriter.write_api_docs(outdir)
    docwriter.write_index(outdir, "index", relative_to=outdir)
    print('%d files written to %s' % (len(docwriter.written_modules),
                                      os.path.abspath(outdir)))
