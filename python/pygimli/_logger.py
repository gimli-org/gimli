# -*- coding: utf-8 -*-
"""pyGIMLi - Stuff for logging

    use with
    pg.info()
    pg.warn()
    pg.debug()
    pg.error() # non critical
    pg.critical() # raises exception

"""

import sys
import logging
import inspect
import traceback

from . import core

__ANSICOLORS__ = {'r':'\033[0;31;49m', #normal, #FG red; #BG black
                  'g':'\033[0;32;49m',
                  'y':'\033[0;33;49m',
                  'b':'\033[0;34;49m',
                  'm':'\033[0;35;49m',
                   'DARK_GREY':'\033[1;30m',
                   'NC':'\033[0m', # No Color
                    }

def _msg(*args):
    msg = ''
    for i, arg in enumerate(args):
        msg += str(arg)
        if i < len(args)-1:
            msg += ' '
    return msg
    
def _(*args, c=None):
    # will probably only for linux or any msys like shell
    if c is None:
        return _msg(*args)
    elif '\033' in c:
        return c + _msg(*args) + __ANSICOLORS__['NC']
    try:
        return __ANSICOLORS__[c] + _msg(*args) + __ANSICOLORS__['NC']
    except:
        return '\033[' + c + 'm' + _msg(*args) + __ANSICOLORS__['NC']

def _get_class_from_frame(fr):
    args, _, _, value_dict = inspect.getargvalues(fr)
    if len(args) and args[0] == 'self':
        instance = value_dict.get('self', None)
        if instance:
            return getattr(instance, '__class__', None)
    return None

def whereAmI(nr=2):
    clsName = _get_class_from_frame(inspect.stack()[nr][0])
    method = inspect.stack()[nr][3]
    return str(clsName) + '.' + method

def p(*args, c='y'):
    deprecated(hint='use _d instead')
    print(_(*args, c=c))

def _g(*args):
    _d(*args, c='g')

def _y(*args):
    _d(*args, c='y')

def _r(*args):
    _d(*args, c='r')

def _d(*args, c='y'):
    """Simplistic colored debug msg"""
    print(_(whereAmI(), ':', *args, c=c))


class ColorFormatter(logging.Formatter):
    def __init__(self, fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s'):
        logging.Formatter.__init__(self, fmt, "%d/%m/%y - %H:%M:%S")
        self._formatter = {
            logging.INFO: logging.PercentStyle('%(asctime)s - %(name)s - ' + _('%(levelname)s', c='g') + ' - %(message)s'),
            logging.ERROR: logging.PercentStyle('%(asctime)s - %(name)s - ' + _('%(levelname)s', c='r') + ' - %(message)s'),
            logging.WARNING: logging.PercentStyle('%(asctime)s - %(name)s - ' + _('%(levelname)s', c='y') + ' - %(message)s'),
            logging.CRITICAL: logging.PercentStyle('%(asctime)s - %(name)s - ' + _('%(levelname)s', c='5;31;1;49') + ' - %(message)s'),
            logging.DEBUG: logging.PercentStyle('%(asctime)s - %(name)s - ' + _('%(levelname)s', c='m') + ' - %(message)s'),
            'DEFAULT': logging.PercentStyle('%(asctime)s - %(name)s - %(levelname)s - %(message)s'),
        }

    def format(self, record):
        self._style = self._formatter.get(record.levelno, self._formatter['DEFAULT'])
        return logging.Formatter.format(self, record)

logger = logging.getLogger('pyGIMLi')
streamHandler = logging.StreamHandler()
streamHandler.setFormatter(ColorFormatter())
logger.root.addHandler(streamHandler)


def addLogLevel(value, name):
    """Add a new log level to the :mod:`logging` module.
    
    Parameters
    ----------
    value: int
        log level number.
    name: str
        Name for the log level.
    """
    logging.addLevelName(value, name)
    setattr(logging, name, value)

# CRITICAL = 50
# ERROR = 40
# WARNING = 30
# INFO = 20
VERBOSE = 15
# DEBUG = 10
# NOTSET = 0
addLogLevel(VERBOSE, 'VERBOSE')
def __logVerbose(msg, *args, **kwargs):
    if logger.isEnabledFor(VERBOSE):
        logger._log(VERBOSE, msg, args, **kwargs)
logger.verbose = __logVerbose

__verbose_level__ = 0

def setVerbose(v):
    level = logging.INFO
    if v:
        __verbose_level__ = 1
        level = logging.VERBOSE
    else:
        __verbose_level__ = 0
        level = logging.INFO
    logger.setLevel(level)

def v(funct):
    """Decorator to enable verbose messages for the scope of a function.

    Examples
    --------
    >>> import pygimli as pg                                    
    >>> @pg.v                                                  # doctest: +SKIP
    >>> def foo():                                             # doctest: +SKIP
    ...     pg.verbose('foo')                                  # doctest: +SKIP 
    >>> foo()                                                  # doctest: +SKIP
    >>> def bar(d):
    ...     pg.verbose('bar', d)
    >>> bar('verbose should be off')
    >>> pg.setVerbose(1)
    >>> bar('verbose should be on (1)')
    >>> pg.setVerbose(0)
    >>> pg.v(bar)('verbose should be on (2)')
    """
    def wrapper(*args, **kwargs):
        o = logger.level
        logger.setLevel(logging.VERBOSE)
        rv = funct(*args, **kwargs)
        logger.setLevel(o)
        return rv
    return wrapper

def d(funct):
    """Decorator to enable debug messages for the scope of a function.
    
    Examples
    --------
    >>> import pygimli as pg                                   
    >>> @pg.d                                                  # doctest: +SKIP
    >>> def foo():                                             # doctest: +SKIP
    ...     pg.debug('foo')                                    # doctest: +SKIP
    >>> foo()                                                  # doctest: +SKIP
    >>> def bar(d):
    ...     pg.debug('bar', d)
    >>> bar('debug should be off')
    >>> pg.setDebug(1)
    >>> bar('debug should be on (1)')
    >>> pg.setDebug(0)
    >>> pg.d(bar)('debug should be on (2)')
    """
    def wrapper(*args, **kwargs):
        o = logger.level
        logger.setLevel(logging.DEBUG)
        rv = funct(*args, **kwargs)
        logger.setLevel(o)
        return rv
    return wrapper

def verbose():
    return __verbose_level__

def setDebug(d):
    level = logging.INFO
    if d:
        core._pygimli_.setDebug(True)
        level = logging.DEBUG
    else:
        core._pygimli_.setDebug(False)
        level = logging.INFO

    logger.setLevel(level)
    logging.getLogger('Core').setLevel(level)
    logging.basicConfig(level=level,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S',
                        #filename='pygimli.log'
                    )

if '--debug' in sys.argv or '-d' in sys.argv:
    setDebug(True)
else:
    setDebug(False)

if '--verbose' in sys.argv or '-v' in sys.argv:
    setVerbose(True)
else:
    setVerbose(False)


def info(*args):
    logger.info(_msg(*args))

def warn(*args):
    logger.warning(_msg(*args))

def error(*args):
    logger.error(whereAmI() + "\n" + _msg(*args))

def debug(*args, withTrace=False):
    """
    Parameters
    ----------
    """
    if withTrace:
        traceback.print_exc()
    logger.debug(_msg(*args))

def verbose(*args):
    logger.verbose(_msg(*args))

def critical(*args):
    logger.critical(whereAmI() + "\n" + _msg(*args))
    raise Exception(_msg(*args))

def deprecated(msg='', hint=''):
    logger.warning(whereAmI() + "\n" + msg + ", is deprecated, please use:" + hint + " instead.")

def renameKwarg(old, new, kwargs):
    if old in kwargs:
        logger.warning("Keyword argument name changed from '" + old + \
                 "' to '" + new + "'")
        kwargs[new] = kwargs.pop(old)

def warnNonEmptyArgs(kwargs):
    if len(kwargs) > 0:
        logger.warning(whereAmI() + "Unrecognized keyword arguments for method:" + _msg(kwargs))
