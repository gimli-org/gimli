"""
An adapted include directive with an optional preprocessor step.

By default, a python script will be converted that it fits into a literate reST 
with the matplotlib plot directive using :context:.

All the code between "'''" "'''" and #! will be interpreted as reST, the rest is puted in 
the plot directive

Using:
       
.. literate:: python.py

Options
-------

The ``literate`` directive supports the following options:

no options so far

Configuration options
---------------------

no options so far
"""

import sys, os, glob, shutil, imp, warnings, io, re, textwrap
import traceback
   
#import exceptions

from docutils import io, nodes, statemachine, utils
from docutils.parsers.rst import directives
from docutils.parsers.rst.directives.misc import Include as BaseInclude
import sphinx

def preProcessLines(rawtext):
    '''
    '''
    
    def startCodeBlock(text):
            text.append('')
            text.append('.. plot::')
            text.append('    :context:')
            text.append('    :include-source:')
            text.append('')
          
       
    isComment = False
    text = []
    for line in rawtext:
        if '#!/' in line:
            continue
        elif '#!' in line:
            line = line.replace("#!", '').lstrip()
            text.append( line )
            
            startCodeBlock( text )
            continue
        elif (isComment == False) and (("'''" in line) or ('"""' in line)):
            isComment = True
            continue
        elif (isComment == True) and (("'''" in line) or ('"""' in line)):
            isComment = False
            startCodeBlock( text )
            continue
                
        if not isComment:
            text.append('    ' +  line )
        else:
            text.append( line )
        
    fi = open( 'tmp.txt', 'w')
    for l in text:
        fi.write( l + '\n')
    fi.close()
    
    return text
    
#def preProcessLines( ... )

class MyLiterateInclude(BaseInclude):
    """
    Like the standard "Include" directive, but interprets absolute paths
    "correctly", i.e., relative to source directory.
    """
    
    def run(self):
        env = self.state.document.settings.env

        if self.arguments[0].startswith('<') and \
           self.arguments[0].endswith('>'):
            # docutils "standard" includes, do not do path processing
            return BaseInclude.run(self)
        rel_filename, filename = env.relfn2path(self.arguments[0])
        self.arguments[0] = filename
        
        if not self.state.document.settings.file_insertion_enabled:
            raise self.warning('"%s" directive disabled.' % self.name)
        source = self.state_machine.input_lines.source(
            self.lineno - self.state_machine.input_offset - 1)
        source_dir = os.path.dirname(os.path.abspath(source))
        path = directives.path(self.arguments[0])
        if path.startswith('<') and path.endswith('>'):
            path = os.path.join(self.standard_include_path, path[1:-1])
        path = os.path.normpath(os.path.join(source_dir, path))
        path = utils.relative_path(None, path)
        path = nodes.reprunicode(path)
        encoding = self.options.get(
            'encoding', self.state.document.settings.input_encoding)
        tab_width = self.options.get(
            'tab-width', self.state.document.settings.tab_width)
        try:
            self.state.document.settings.record_dependencies.add(path)
            include_file = io.FileInput(
                source_path=path, encoding=encoding,
                error_handler=(self.state.document.settings.\
                               input_encoding_error_handler),
                handle_io_errors=None)
        except IOError as error:
            raise self.severe('Problems with "%s" directive path:\n%s.' %
                      (self.name, error))
        startline = self.options.get('start-line', None)
        endline = self.options.get('end-line', None)
        try:
            if startline or (endline is not None):
                lines = include_file.readlines()
                rawtext = ''.join(lines[startline:endline])
            else:
                rawtext = include_file.read()
        except UnicodeError as error:
            raise self.severe('Problem with "%s" directive:\n%s' %
                              (self.name, ErrorString(error)))
        # start-after/end-before: no restrictions on newlines in match-text,
        # and no restrictions on matching inside lines vs. line boundaries
        after_text = self.options.get('start-after', None)
        if after_text:
            # skip content in rawtext before *and incl.* a matching text
            after_index = rawtext.find(after_text)
            if after_index < 0:
                raise self.severe('Problem with "start-after" option of "%s" '
                                  'directive:\nText not found.' % self.name)
            rawtext = rawtext[after_index + len(after_text):]
        before_text = self.options.get('end-before', None)
        if before_text:
            # skip content in rawtext after *and incl.* a matching text
            before_index = rawtext.find(before_text)
            if before_index < 0:
                raise self.severe('Problem with "end-before" option of "%s" '
                                  'directive:\nText not found.' % self.name)
            rawtext = rawtext[:before_index]
        if 'literal' in self.options:
            # Convert tabs to spaces, if `tab_width` is positive.
            if tab_width >= 0:
                text = rawtext.expandtabs(tab_width)
            else:
                text = rawtext
            literal_block = nodes.literal_block(rawtext, text, source=path)
            literal_block.line = 1
            return [literal_block]
        else:
            include_lines = statemachine.string2lines(
                rawtext, tab_width, convert_whitespace=1)
            
            include_lines = preProcessLines( include_lines )
            
            self.state_machine.insert_input(include_lines, path)
            return []
                       
def setup(app):
    setup.app = app
    setup.config = app.config
    setup.confdir = app.confdir

    app.add_directive('literate', MyLiterateInclude )