# -*- coding: utf-8 -*-
"""A special directive for inline code execution.

The source code inline snipped may be included by:

  2. Included as inline content to the directive::

     .. exec::

        print 1

"""

import sys

from sphinx.util.compat import Directive

class RedirectOutput:
    def __init__(self, console ):
        self.buff = []

        if console == "cout": #sys.stdout
            self.style = 1 #black
        elif console == "cerr": #sys.stderr
            self.style = 2 #red

    def write(self, what):
        ##Do not put a print statement here!!##
        #self.logFile.write( what )
        self.buff.append( what)

class ExecDirective( Directive ):
#class ExecDirective(name, arguments, options, content, lineno,
#                    content_offset, block_text, state, state_machine ):

    print "ExecDirective(Directive):", Directive

    has_content = True

    def run( self ):
        print "run() ExecDirective(Directive):", Directive
        print "name:", self.name
        print "arguments:", self.arguments
        print "options:", self.options
        print "content:", self.content
        print "lineno:", self.lineno
        print "content_offset:", self.content_offset
        print "block_text:", self.block_text
        print "state:", self.state
        print "state_machine:", self.state_machine

        # Now start generating the lines of output
        lines = []

        #lines.extend( ['>>> %s' % self.content] )
        #lines.extend(['::', ''])
        #lines.extend(['>>>%s' % row.rstrip()  for row in self.block_text.split('\n')])

        tmpSysOut = sys.stdout
        tmpSysErr = sys.stderr
        sys.stdout = RedirectOutput( "cout" )
        sys.stderr = RedirectOutput( "cerr" )

        #lines.extend(['::',''])

        # ADD CODING rule somewhere
        # -*- coding: utf-8 -*-
        for row in self.content:
            if len( row ) == 0:
                lines.extend(['\n' ])
                continue

            #if row[0] == '#':
                #lines.extend( [row[1::]+ '\n'] )
                #continue

            lines.extend(['>>> %s' % row.rstrip() ])
            exec( row )

            if sys.stdout.buff:
                for o in sys.stdout.buff:
                    lines.extend(['%s' % o ])
                sys.stdout.buff = []

            if sys.stderr.buff:
                for o in sys.stderr.buff:
                    lines.extend(['%s' % o ])
                    sys.stderr.buff = []

        sys.stdout = tmpSysOut
        sys.stderr = tmpSysErr

        self.state_machine.insert_input( lines, self.state_machine.input_lines.source(0))

        return []

#def plot_directive(name, arguments, options, content, lineno, content_offset, block_text, state, state_machine):

    #print "plot_directive(name, arguments, options, content, lineno,content_offset, block_text, state, state_machine):"

def setup(app):
    #app.add_config_value('todo_include_todos', False, False)

    #app.add_node(execNode,
                #html=(visit_exec_node, depart_exec_node),
                #latex=(visit_exec_node, depart_exec_node),
                #text=(visit_exec_node, depart_exec_node))

    app.add_directive('exec', ExecDirective)
    #app.add_directive('todolist', TodolistDirective)
    #app.connect('doctree-resolved', process_exec_nodes)
    #app.connect('env-purge-doc', purge_todos)