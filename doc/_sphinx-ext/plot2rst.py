"""
Taken from: http://tonysyu.github.io/mpltools/index.html
Modified by CR

Generate reStructuredText example from python files.

Generate the rst files for the examples by iterating over the python
example files. Files that generate images should start with 'plot'.

To generate your own examples, add ``'mpltools.sphinx.plot2rst'`` to the list
of ``extensions`` in your Sphinx configuration file. In addition, make sure the
example directory(ies) in `plot2rst_paths` (see below) points to a directory
with examples named `plot_*.py` and include an `index.rst` file.

This code was adapted from scikits-image, which took it from scikits-learn.


Options
-------
The ``plot2rst`` extension accepts the following options:

plot2rst_paths : length-2 tuple, or list of tuples
    Tuple or list of tuples of paths to (python plot, generated rst) files,
    i.e. (source, destination).  Note that both paths are relative to Sphinx
    'source' directory. Defaults to ('../examples', 'auto_examples')

plot2rst_rcparams : dict
    Matplotlib configuration parameters. See
    http://matplotlib.sourceforge.net/users/customizing.html for details.

plot2rst_default_thumb : str
    Path (relative to doc root) of default thumbnail image.

plot2rst_thumb_scale : float
    Scale factor for thumbnail. Defaults to 0.2, which scales the thumbnail to
    1/5th the original size.

plot2rst_plot_tag : str
    When this tag is found in the example file, the current plot is saved and
    tag is replaced with plot path. Defaults to 'PLOT2RST.current_figure'.

plot2rst_index_name : str
    The basename for gallery index file. Each example directory should have an
    index file---typically containing nothing more than a simple header. Note
    that the reStructuredText extension (e.g., 'txt', 'rst') is taken from the
    default extension set for the Sphinx project, so you should not specify
    the extension as part of this name. Defaults to 'index'.

plot2rst_flags : dict
    Flags that can be set in gallery indexes or python example files. See
    Flags_ section below for details.

Flags
-----
You can add flags to example files by added a code comment with the prefix ``#PLOT2RST:`. Flags are specified as key-value pairs; for example::

    #PLOT2RST: auto_plots = False

There are also reStructuredText flags, which can be added like::

    .. plot2rst_gallery_style:: list

Some flags can only be defined in the python example file, while others can only be defined in the gallery index. The following flags are defined:

auto_plots : bool
    If no plot tags are found in the example, `auto_plots` adds all figures to
    the end of the example. Defaults to True.

plot2rst_gallery_style : {'thumbnail' | 'list'}
    Display examples as a thumbnail gallery or as a list of titles. This option
    can also be set at the directory-level by adding::

        .. plot2rst_gallery_style:: list

    to the gallery index. Defaults to 'thumbnail'.


Note: If flags are added at the top of the file, then they are stripped from
the resulting reStructureText output. If they appear after the first text or
code block, then will show up in the generated example.


Suggested CSS definitions
-------------------------

    div.body h2 {
        border-bottom: 1px solid #BBB;
        clear: left;
    }

    /*---- example gallery ----*/

    .gallery.figure {
        float: left;
        width: 200px;
        height: 200px;
    }

    .gallery.figure img{
        display: block;
        margin-left: auto;
        margin-right: auto;
        width: 180px;
    }

    .gallery.figure .caption {
        text-align: center !important;
    }

"""
import re
import sys
import os
import shutil
import token
import tokenize

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import image


CODE_LINK = """
.. only:: html

    **Python source code:** :download:`download <{0}>`
    (generated using |version|)

"""

TOCTREE_TEMPLATE = """
.. toctree::
   :hidden:

   %s

"""

IMAGE_TEMPLATE = """
.. image:: images/%s
    :align: center
    :scale: 75

"""

GALLERY_IMAGE_TEMPLATE = """
.. figure:: %(thumb)s
   :figclass: gallery
   :target: ./%(source)s.html

   :ref:`example_%(link_name)s`

"""

ANIMATION_TEMPLATE = '''
.. raw:: html

   <video controls="controls">
       <source src="%s"
               type="video/mp4" />
       Your browser does not support the video tag.
   </video>
'''


GALLERY_LIST_TEMPLATE = """
* :ref:`example_%(link_name)s`
"""


FLAG_PREFIX = '#PLOT2RST:'

class RedirectOutput:
    def __init__(self, console):
        self.buff = []

        if console == "cout": #sys.stdout
            self.style = 1 #black
        elif console == "cerr": #sys.stderr
            self.style = 2 #red

    def write(self, what):
        ##Do not put a print statement here!!##
        #self.logFile.write( what )
        pass
        #self.buff.append(what)



class Path(str):
    """Path object for manipulating directory and file paths."""

    def __init__(self, path):
        super(Path, self).__init__()
        #super(Path, self).__init__(path)

    @property
    def isdir(self):
        return os.path.isdir(self)

    @property
    def exists(self):
        """Return True if path exists"""
        return os.path.exists(self)

    def pjoin(self, *args):
        """Join paths. `p` prefix prevents confusion with string method."""
        return self.__class__(os.path.join(self, *args))

    def psplit(self):
        """Split paths. `p` prefix prevents confusion with string method."""
        return [self.__class__(p) for p in os.path.split(self)]

    def makedirs(self):
        if not self.exists:
            os.makedirs(self)

    def listdir(self):
        return os.listdir(self)

    def format(self, *args, **kwargs):
        return self.__class__(super(Path, self).format(*args, **kwargs))

    def __add__(self, other):
        return self.__class__(super(Path, self).__add__(other))

    def __iadd__(self, other):
        return self.__add__(other)


def setup(app):
    app.connect('builder-inited', generate_example_galleries)

    app.add_config_value('plot2rst_paths',
                         ('../examples', 'auto_examples'), True)
    app.add_config_value('plot2rst_rcparams', {}, True)
    app.add_config_value('plot2rst_default_thumb', None, True)
    app.add_config_value('plot2rst_thumb_scale', 0.2, True)
    app.add_config_value('plot2rst_plot_tag', 'PLOT2RST.current_figure', True)
    app.add_config_value('plot2rst_anim_tag', 'PLOT2RST.current_anim', True)
    app.add_config_value('plot2rst_index_name', 'index', True)
    app.add_config_value('plot2rst_gallery_style', 'thumbnail', True)
    app.add_config_value('plot2rst_commandTranslator', '', 'html')
    # NOTE: plot2rst_flags gets set with defaults later so that keys that are
    # not set in config file still get set to the desired defaults.
    app.add_config_value('plot2rst_flags', {}, True)


def generate_example_galleries(app):
    cfg = app.builder.config

    default_flags = {'auto_plots': True,
                     'gallery_style': 'thumbnail'}
    default_flags.update(cfg.plot2rst_flags)
    cfg.plot2rst_flags = default_flags

    doc_src = Path(os.path.abspath(app.builder.srcdir)) # path/to/doc/source

    if isinstance(cfg.plot2rst_paths, tuple):
        cfg.plot2rst_paths = [cfg.plot2rst_paths]
    for src_dest in cfg.plot2rst_paths:
        plot_path, rst_path = [Path(p) for p in src_dest]
        example_dir = doc_src.pjoin(plot_path)
        rst_dir = doc_src.pjoin(rst_path)
        generate_examples_and_gallery(example_dir, rst_dir, cfg)


def generate_examples_and_gallery(example_dir, rst_dir, cfg):
    """Generate rst from examples and create gallery to showcase examples."""
    if not example_dir.exists:
        print(("No example directory found at", example_dir))
        return
    rst_dir.makedirs()

    # we create an index.rst with all examples
    gallery_index = open(rst_dir.pjoin('index'+cfg.source_suffix), 'w')

    # Here we don't use an os.walk, but we recurse only twice: flat is
    # better than nested.
    write_gallery(gallery_index, example_dir, rst_dir, cfg)
    for d in sorted(example_dir.listdir()):
        example_sub = example_dir.pjoin(d)
        if example_sub.isdir:
            rst_sub = rst_dir.pjoin(d)
            rst_sub.makedirs()
            write_gallery(gallery_index, example_sub, rst_sub, cfg, depth=1)
    gallery_index.flush()

def write_gallery(gallery_index, src_dir, rst_dir, cfg, depth=0):
    """Generate the rst files for an example directory, i.e. gallery.

    Write rst files from python examples and add example links to gallery.

    Parameters
    ----------
    gallery_index : file
        Index file for plot gallery.
    src_dir : 'str'
        Source directory for python examples.
    rst_dir : 'str'
        Destination directory for rst files generated from python examples.
    cfg : config object
        Sphinx config object created by Sphinx.
    """
    index_name = cfg.plot2rst_index_name + cfg.source_suffix
    gallery_template = src_dir.pjoin(index_name)
    if not os.path.exists(gallery_template):
        print(src_dir)
        print((80*'_'))
        print(('Example directory %s does not have a %s file'
                        % (src_dir, index_name)))
        print('Skipping this directory')
        print((80*'_'))
        return
    flags = get_flags_from_rst(gallery_template)

    gallery_description = open(gallery_template).read()
    gallery_index.write('\n\n%s\n\n' % gallery_description)

    rst_dir.makedirs()
    examples = [fname for fname in sorted(src_dir.listdir(), key=_plots_first)
                      if fname.endswith('.py')]

    ex_names = [ex[:-3] for ex in examples] # strip '.py' extension
    if depth == 0:
        sub_dir = Path('')
    else:
        sub_dir_list = src_dir.psplit()[-depth:]
        sub_dir = Path('/'.join(sub_dir_list) + '/')

    for src_name in examples:
        write_example(src_name, src_dir, rst_dir, cfg)

        link_name = sub_dir.pjoin(src_name)
        link_name = link_name.replace(os.path.sep, '_')
        if link_name.startswith('._'):
            link_name = link_name[2:]

        info = {}
        info['source'] = sub_dir + src_name[:-3]
        info['link_name'] = link_name

        gallery_style = flags.get('plot2rst_gallery_style',
                                  cfg.plot2rst_flags['gallery_style'])
        if gallery_style == 'thumbnail':
            thumb_name = src_name[:-3] + '.png'
            info['thumb'] = sub_dir.pjoin('images/thumb', thumb_name)
            gallery_index.write(GALLERY_IMAGE_TEMPLATE % info)
        elif gallery_style == 'list':
            gallery_index.write(GALLERY_LIST_TEMPLATE % info)

    tempex_names = str()
    for ex_name in ex_names:
        tempex_names += sub_dir + ex_name + "\n   "
    gallery_index.write(TOCTREE_TEMPLATE % (tempex_names))

    # update _templates/gallery.html
    BUILDDIR='_build/'
    os.system("python " + BUILDDIR + "sidebar_gallery.py")


def get_flags_from_rst(rst_file):
    """Return dict of plot2rst flags found in reStructuredText file.

    Flags should have the form:

        .. plot2rst_*:: value

    """
    flags = {}
    with open(rst_file) as f:
        for line in f:
            if line.startswith('.. plot2rst'):
                line = line.lstrip('.. ')
                k, v = [word.strip() for word in line.split('::')]
                flags[k] = v
    return flags


def _plots_first(fname):
    """Decorate filename so that examples with plots are displayed first."""
    if not (fname.startswith('plot') and fname.endswith('.py')):
        return 'zz' + fname
    return fname


def write_example(src_name, src_dir, rst_dir, cfg):
    """Write rst file from a given python example.

    Parameters
    ----------
    src_name : str
        Name of example file.
    src_dir : 'str'
        Source directory for python examples.
    rst_dir : 'str'
        Destination directory for rst files generated from python examples.
    cfg : config object
        Sphinx config object created by Sphinx.
    """
    last_dir = src_dir.psplit()[-1]
    # to avoid leading . in file names, and wrong names in links
    if last_dir == '.' or last_dir == 'examples':
        last_dir = Path('')
    else:
        last_dir += '_'

    src_path = src_dir.pjoin(src_name)
    example_file = rst_dir.pjoin(src_name)
    shutil.copyfile(src_path, example_file)

    image_dir = rst_dir.pjoin('images')
    thumb_dir = image_dir.pjoin('thumb')
    image_dir.makedirs()
    thumb_dir.makedirs()

    # copy static files
    static_src = src_dir.pjoin('static')
    if static_src.exists:
        static_dir = rst_dir.pjoin('static')
        if static_dir.exists:
            shutil.rmtree(static_dir)
        shutil.copytree(static_src, static_dir)

    base_image_name = os.path.splitext(src_name)[0]
    image_path = image_dir.pjoin(base_image_name + '_{0}.png')

    basename, py_ext = os.path.splitext(src_name)
    rst_path = rst_dir.pjoin(basename + cfg.source_suffix)

    if _plots_are_current(src_path, image_path) and rst_path.exists:
        return

    flags = cfg.plot2rst_flags.copy()
    blocks, new_flags = split_code_and_text_blocks(example_file)
    flags.update(new_flags)

    while True:
        head = blocks[0][2]
        if head.startswith('#!') or head.startswith(FLAG_PREFIX):
            blocks.pop(0) # don't add shebangs or flags to rst file.
        else:
            break

    # Note that `process_blocks` executes the source, so plots are now 'active'
    figure_list, rst = process_blocks(blocks, src_path, image_path, cfg)

    rst_link = '.. _example_%s:\n\n' % (last_dir + src_name)
    example_rst = ''.join([rst_link, rst])

    has_inline_plots = any(cfg.plot2rst_plot_tag in b[2] for b in blocks)

    if not has_inline_plots and flags['auto_plots']:
        # Show all plots at the end of the example
        if len(plt.get_fignums()) > 0:
            figure_list = save_all_figures(image_path)
            img_blocks = [IMAGE_TEMPLATE % f.lstrip('/') for f in figure_list]
            example_rst += ''.join(img_blocks)
    plt.close('all')

    example_rst += CODE_LINK.format(src_name)

    f = open(rst_path,'w')
    f.write(example_rst)
    f.flush()

    thumb_path = thumb_dir.pjoin(src_name[:-3] + '.png')
    if figure_list:
        first_image_file = image_dir.pjoin(figure_list[0].lstrip('/'))
        if first_image_file.exists:
            thumb_scale = cfg.plot2rst_thumb_scale
            image.thumbnail(first_image_file, thumb_path, thumb_scale)

    if not thumb_path.exists:
        if cfg.plot2rst_default_thumb is None:
            print("WARNING: No plots found and default thumbnail not defined.")
            print("Specify 'plot2rst_default_thumb' in Sphinx config file.")
        else:
            shutil.copy(cfg.plot2rst_default_thumb, thumb_path)


def _plots_are_current(src_path, image_path):
    first_image_file = Path(image_path.format(1))
    needs_replot = (not first_image_file.exists or
                    _mod_time(first_image_file) <= _mod_time(src_path))
    return not needs_replot


def _mod_time(file_path):
    return os.stat(file_path).st_mtime


def split_code_and_text_blocks(source_file):
    """Return list with source file separated into code and text blocks.

    Returns
    -------
    blocks : list of (label, (start, end+1), content)
        List where each element is a tuple with the label ('text' or 'code'),
        the (start, end+1) line numbers, and content string of block.
    flags : dict
        Option flags for plot2rst that were found in the source file.
    """
    block_edges, idx_first_text_block, flags = analyze_blocks(source_file)

    with open(source_file) as f:
        source_lines = f.readlines()

    if idx_first_text_block is None:
        blocks = [('code', (1, len(source_lines)), ''.join(source_lines))]
        return blocks, flags

    # Every other block should be a text block
    idx_text_block = np.arange(idx_first_text_block, len(block_edges), 2)
    blocks = []
    slice_ranges = list(zip(block_edges[:-1], block_edges[1:]))
    for i, (start, end) in enumerate(slice_ranges):
        block_label = 'text' if i in idx_text_block else 'code'
        # subtract 1 from indices b/c line numbers start at 1, not 0
        content = ''.join(source_lines[start-1:end-1])
        blocks.append((block_label, (start, end), content))
    return blocks, flags


def analyze_blocks(source_file):
    """Return starting line numbers of code and text blocks

    Returns
    -------
    block_edges : list of int
        Line number for the start of each block. Note the
    idx_first_text_block : {0 | 1}
        0 if first block is text then, else 1 (second block better be text).
    flags : dict
        Option flags for plot2rst that were found in the source file.
    """
    flags = {}
    block_edges = []
    with open(source_file) as f:
        token_iter = tokenize.generate_tokens(f.readline)

        for token_tuple in token_iter:
            t_id, t_str, (srow, scol), (erow, ecol), src_line = token_tuple
            tok_name = token.tok_name[t_id]
            if tok_name == 'STRING' and scol == 0:
                # Add one point to line after text (for later slicing)
                block_edges.extend((srow, erow+1))
            elif tok_name == 'COMMENT' and t_str.startswith(FLAG_PREFIX):
                flag_args = t_str.lstrip(FLAG_PREFIX).split('=')
                if not len(flag_args) == 2:
                    raise ValueError("Flags must be key-value pairs.")
                key = flag_args[0].strip()
                flags[key] = eval(flag_args[1])
    idx_first_text_block = 0
    if not block_edges: # no text blocks
        idx_first_text_block = None
    else:
        # when example doesn't start with text block.
        if not block_edges[0] == 1:
            block_edges.insert(0, 1)
            idx_first_text_block = 1
        # when example doesn't end with text block.
        if not block_edges[-1] == erow: # iffy: I'm using end state of loop
            block_edges.append(erow)
    return block_edges, idx_first_text_block, flags


def process_blocks(blocks, src_path, image_path, cfg):
    """Run source, save plots as images, and convert blocks to rst.

    Parameters
    ----------
    blocks : list of block tuples
        Code and text blocks from example. See `split_code_and_text_blocks`.
    src_path : str
        Path to example file.
    image_path : str
        Path where plots are saved (format string which accepts figure number).
    cfg : config object
        Sphinx config object created by Sphinx.

    Returns
    -------
    figure_list : list
        List of figure names saved by the example.
    rst_text : str
        Text with code wrapped code-block directives.
    """
    src_dir, src_name = src_path.psplit()
    if not src_name.startswith('plot'):
        convert_func = dict(code=codestr2rst, text=docstr2rst)
        rst_blocks = [convert_func[blabel](bcontent)
                      for i, (blabel, brange, bcontent) in enumerate(blocks)]
        return [], '\n'.join(rst_blocks)

    # index of blocks which have inline plots
    inline_tag = cfg.plot2rst_plot_tag
    idx_inline_plot = [i for i, b in enumerate(blocks)
                       if inline_tag in b[2]]

    image_dir, image_fmt_str = image_path.psplit()

    figure_list = []
    plt.rcdefaults()
    plt.rcParams.update(cfg.plot2rst_rcparams)
    plt.close('all')

    example_globals = {}
    rst_blocks = []
    fig_num = 1
    anim_num = 1
    lastCoutBuff = []

    tmpSysOut = sys.stdout
    tmpSysErr = sys.stderr
            
    print("Processing:", src_path)
    plt.ion()
    for i, (blabel, brange, bcontent) in enumerate(blocks):
        if blabel == 'code':

            print('-'*100)
            print(bcontent, example_globals)

            #tmpSysOut = sys.stdout
            #tmpSysErr = sys.stderr
            
            sys.stdout = RedirectOutput("cout")
            sys.stderr = RedirectOutput("cerr")

            exec(bcontent, example_globals)
            rst_blocks.append(codestr2rst(bcontent))

            if len(sys.stdout.buff) > 0:
                lastCoutBuff = sys.stdout.buff
            coutbuff = sys.stdout.buff
            cerrbuff = sys.stderr.buff

            sys.stdout = tmpSysOut
            sys.stderr = tmpSysErr

            print('#'*100)
            print("coutbuf:", coutbuff)
            print("cerrbuf:", cerrbuff)

            if len(coutbuff) > 0:
                rst_blocks.append(printcout2rst(coutbuff))

            if len(cerrbuff) > 0:
                rst_blocks.append(printcerr2rst(cerrbuff))
            print('-'*100)

        else:
            if i in idx_inline_plot:
                plt.savefig(image_path.format(fig_num))
                figure_name = image_fmt_str.format(fig_num)
                fig_num += 1
                figure_list.append(figure_name)
                figure_link = os.path.join('images', figure_name)
                bcontent = bcontent.replace(inline_tag, figure_link)


            if '.. lastcout::' in bcontent:
                bcontent = bcontent.replace('.. lastcout::',
                                            printcout2rst(lastCoutBuff))
                lastCoutBuff = []

            if '.. animate::' in bcontent:
                bcontent = bcontent.replace('"""', '')
                vals = bcontent.split()
                #print(vals)
                animator = vals[2]
                anim_name = image_path.format(anim_num).replace('.png', '.mp4')
                args = vals[3:]
                print(args)
                exec(animator + '.save("'+anim_name+'", ' + \
                    ','.join(args) + ' ,extra_args=["-vcodec", "libx264"])', example_globals)
                rst_blocks.append(ANIMATION_TEMPLATE % (anim_name))
                anim_num += 1
                continue

            rst_blocks.append(docstr2rst(bcontent, cfg))
    plt.ioff()
    return figure_list, '\n'.join(rst_blocks)


def printcerr2rst(outbuff):
    code_directive = ".. error:: \n"
    indented_block = '\t\n'
    for t in outbuff:
        if t == '\n':
            indented_block += '\t\n'
        else:
            indented_block += '\t\t> ' + t.replace('\n', '\n\t')
    return code_directive + indented_block + "*"

def printcout2rst(outbuff):
    """Return reStructuredText code block from print out string"""
    code_directive = "..\n\n"
    indented_block = '\t\n'
    for t in outbuff:
        if t == '\n':
            indented_block += '\t\n\n'
        elif t == ' ':
            indented_block += ' '
        else:
            indented_block += '\t\t*' + t.replace('\n', '\n\t') + '*'

    return code_directive + indented_block

def codestr2rst(codestr):
    """Return reStructuredText code block from code string"""
    code_directive = ".. code-block:: python\n\n"
    indented_block = '\t' + codestr.replace('\n', '\n\t')
    return code_directive + indented_block


def docstr2rst(docstr, cfg):
    """Return reStructuredText from docstring"""
    idx_whitespace = len(docstr.rstrip()) - len(docstr)
    whitespace = docstr[idx_whitespace:]
    #CR print "##################################################################"
    #CR print docstr
    #CR print "##################################################################"
    #CR print whitespace
    #CR print "##################################################################"
    #CR eval() eats latex commands like \a \t so we need to replace them first

    mathDictionary = {}
    commandDictionary = {}
    command1Dictionary = {}

    if cfg.plot2rst_commandTranslator:
        #print(cfg.plot2rst_commandTranslator)
        mathDictionary = cfg.plot2rst_commandTranslator['mathDictionary']
        commandDictionary = cfg.plot2rst_commandTranslator['commandDictionary']
        command1Dictionary = cfg.plot2rst_commandTranslator['command1Dictionary']

    current = docstr
    for x in mathDictionary:
        current = re.sub('\\\\' + x + '(?!\w)',
                         '\\operatorname{' + mathDictionary[x] + '}', current)
    for x in commandDictionary:
        current = re.sub('\\\\' + x + '(?!\w)', commandDictionary[x], current)
    for x in command1Dictionary:
        n=1
        while n > 0:
            sub = re.search('\\\\'+ x +'{([A-Za-z]*)}', current)
            if sub:
                (current,n) = re.subn('\\\\' + x + '{([A-Za-z]*)}',
                             command1Dictionary[x].replace('#1',sub.group(1)),
                             current,count=1)
                #print(x, current, n, sub.group(1))
            else:
                n = 0
    docstr = current

    #print(docstr)
    #sys.exit()
    return eval(docstr.replace("\\","\\\\")) + whitespace



def save_all_figures(image_path):
    """Save all matplotlib figures.

    Parameters
    ----------
    image_path : str
        Path where plots are saved (format string which accepts figure number).
    """
    figure_list = []
    image_dir, image_fmt_str = image_path.psplit()
    fig_mngr = matplotlib._pylab_helpers.Gcf.get_all_fig_managers()
    for fig_num in (m.num for m in fig_mngr):
        # Set the fig_num figure as the current figure as we can't
        # save a figure that's not the current figure.
        plt.figure(fig_num)
        plt.savefig(image_path.format(fig_num))
        figure_list.append(image_fmt_str.format(fig_num))
    return figure_list

