# -*- coding: utf-8 -*-
"""
Fixed and improved version based on "extracting from C++ doxygen documented file Author G.D." and py++ code.

Distributed under the Boost Software License, Version 1.0. (See
accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)

modified by CR

"""
from functools import reduce
import sys
import traceback
import codecs
  
class doxygen_doc_extractor:
    """
    Extracts Doxygen styled documentation from source or generates it from description.
    """
    def __init__(self):
        #for caching source
        self.file_name = None
        self.source = None
        self.hasError = False
    #__init__

    def __call__(self, declaration):
        doc_lines = []
        
        try:
            if self.file_name != declaration.location.file_name:
                self.file_name = declaration.location.file_name
                self.source = open(declaration.location.file_name, encoding='utf-8').readlines()
            find_block_end = False
            
            # search backward until file begin
            for lcount in range(declaration.location.line-2, -1, -1):
                line = self.source[lcount]

                if not find_block_end:
                    try:
                        #CR print line.rstrip()[-2:]
                        if line.rstrip()[-2:] == "*/":
                            find_block_end = True
                    except Exception as e:
                        if not self.hasError:
                            print('*'*100)
                            print(e)
                            print(line)
                            self.hasError = True
                        pass

                if find_block_end:
                    try:
                        if line.lstrip()[:2] == "/*":
                            find_block_end = False
                    except Exception as e:
                        if not self.hasError:
                            self.hasError = True
                            print('*'*100)
                            print(e)
                            print(line)
                        pass

                final_str = self.clear_str(line)

                if not find_block_end and self.is_code(line):
                    break
                if final_str:
                    doc_lines.insert(0, final_str)
        except Exception as e:
            if not self.hasError:
                traceback.print_exc(file=sys.stdout)
                self.hasError = True
                print('*'*100)
                print(e)
                print(self.file_name)
            pass
        finally:
            if len(doc_lines) > 0:
                final_doc_lines = [ line.replace("\n","\\n") for line in doc_lines[:-1] ]
                final_doc_lines.append(doc_lines[-1].replace("\n",""))
                return '\"' + ''.join(final_doc_lines) + '\"'
            else:
                return '\"\"'
    #__call__()

    def clear_str(self, tmp_str):
        """
        Replace */! by space and \brief, @fn, \param, etc
        """
        #CR: add
        def clean ( _str, sym, change2=""):
            return _str.replace(sym, change2)

        #CR: add '\r\n'
        #tmp_str = reduce(clean, [tmp_str, "\r\n", '/','*','!',
                                 #"\\brief","@brief","\\fn","@fn",
                                 #"\\ref","@ref", "\"", "\'", "\\c"])

        for sym in ['/', "\r\n", '/', '*', '!', "\\brief", "@brief",
                    "\\fn", "@fn", "\\ref", "@ref", "\"", "\'", "\\c"]:
            tmp_str = clean(tmp_str, sym)
        
        tmp_str = clean(tmp_str, '\\', '')
        
        #commands list taken form : http://www.stack.nl/~dimitri/doxygen/commands.html
        replacement_list = [
#           "a",
            "addindex",
            "addtogroup",
            "anchor",
            "arg",
            "attention",
            "author",
#           "b",
#           "brief",
            "bug",
#           "c",
            "callgraph",
            "callergraph",
            "category",
            "class",
            ("code","[Code]"),
            "cond",
            "copybrief",
            "copydetails",
            "copydoc",
            "date",
            "def",
            "defgroup",
            "deprecated",
            "details",
            "dir",
            "dontinclude",
            ("dot","[Dot]"),
            "dotfile",
            "e",
            "else",
            "elseif",
            "em",
            ("endcode","[/Code]"),
            "endcond",
            ("enddot","[/Dot]"),
            "endhtmlonly",
            "endif",
            "endlatexonly",
            "endlink",
            "endmanonly",
            "endmsc",
            "endverbatim",
            "endxmlonly",
            "enum",
            "example",
            "exception",
            "extends",
            "f$",
            "f[",
            "f]",
            "f{",
            "f}",
            "file",
#           "fn",
            "headerfile",
            "hideinitializer",
            "htmlinclude",
            "htmlonly",
            "if",
            "ifnot",
            "image",
            "implements",
            "include",
            "includelineno",
            "ingroup",
            "internal",
            "invariant",
            "interface",
            "latexonly",
            "li",
            "line",
            "link",
            "mainpage",
            "manonly",
            "memberof",
            "msc",
#           "n",
            "name",
            "namespace",
            "nosubgrouping",
            "note",
            "overload",
#           "p",
            "package",
            "page",
            "par",
            "paragraph",
            "param",
            "post",
            "pre",
#           "private (PHP only)",
#           "privatesection (PHP only)",
            "property",
#           "protected (PHP only)",
#           "protectedsection (PHP only)",
            "protocol",
#           "public (PHP only)",
#           "publicsection (PHP only)",
#           "ref",
            "relates",
            "relatesalso",
            "remarks",
            "return",
            "retval",
            "sa",
            "section",
            "see",
            "showinitializer",
            "since",
            "skip",
            "skipline",
            "struct",
            "subpage",
            "subsection",
            "subsubsection",
            "test",
            "throw",
            ("todo","TODO"),
            "tparam",
            "typedef",
            "union",
            "until",
            "var",
            "verbatim",
            "verbinclude",
            "version",
            "warning",
            "weakgroup",
            "xmlonly",
            "xrefitem",
#           "$",
#           "@",
#           "\",
#           "&",
#           "~",
#           "<",
#           ">",
#           "#",
#           "%",
            ]

        for command in replacement_list:
            try:
                old,new = command
            except ValueError:
                old = command
                new = command.capitalize()+":"
            tmp_str = clean(tmp_str, "@"+old, new)
            tmp_str = clean(tmp_str, "\\"+old, new)

        return tmp_str.lstrip()
    #clean_str()

    def is_code(self, tmp_str):
        """
        Detects if tmp_str is code or not
        """
        try:
            beg = tmp_str.lstrip()[:2]
            return beg != "//" and beg != "/*"
        except:
            pass
        return False
    #is_code()

#class doxygen_doc_extractor
