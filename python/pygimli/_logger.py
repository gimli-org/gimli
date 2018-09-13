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
    # will probably only for linux
    if c is None:
        return _msg(*args)
    elif '\033' in c:
        return c + _msg(*args) + __ANSICOLORS__['NC']
    try:
        return __ANSICOLORS__[c] + _msg(*args) + __ANSICOLORS__['NC']
    except:
        return '\033[' + c + 'm' + _msg(*args) + __ANSICOLORS__['NC']


class ColorFormatter(logging.Formatter):
    def __init__(self, fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s'):
        logging.Formatter.__init__(self, fmt)
        self._formatter = {
            logging.INFO: logging.PercentStyle('%(asctime)s - %(name)s - %(levelname)s - %(message)s'),
            logging.ERROR: logging.PercentStyle('%(asctime)s - %(name)s - ' + _('%(levelname)s', c='r') + '- %(message)s'),
            logging.WARNING: logging.PercentStyle('%(asctime)s - %(name)s - ' + _('%(levelname)s', c='y') + '- %(message)s'),
            logging.CRITICAL: logging.PercentStyle('%(asctime)s - %(name)s - ' + _('%(levelname)s', c='5;31;1;49') + '- %(message)s'),
            logging.DEBUG: logging.PercentStyle('%(asctime)s - %(name)s - ' + _('%(levelname)s', c='m') + '- %(message)s'),
            'DEFAULT': logging.PercentStyle('%(asctime)s - %(name)s - %(levelname)s - %(message)s'),
        }

    def format(self, record):
        self._style = self._formatter.get(record.levelno, self._formatter['DEFAULT'])
        return logging.Formatter.format(self, record)

logger = logging.getLogger('pyGIMLi')
streamHandler = logging.StreamHandler()
streamHandler.setFormatter(ColorFormatter())
logger.root.addHandler(streamHandler)


def setDebug(d):
    level = logging.INFO
    if d:
        core._pygimli_.setDebug(True)
        level = logging.DEBUG
        logger.debug("Set debug mode: on")
    else:
        core._pygimli_.setDebug(False)
        level = logging.INFO
        logger.debug("Set debug mode: off")
    
    logger.setLevel(level)
    logging.getLogger('Core').setLevel(level)
    logging.basicConfig(level=level,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%m/%d/%Y %H:%M:%S',
                    #filename='pygimli.log'
                    )

if '--debug' in sys.argv:
    setDebug(True)
else:
    setDebug(False)

def info(*args):
    logger.info(_msg(*args))

def warn(*args):
    logger.warning(_msg(*args))

def error(*args):
    caller = sys._getframe(1).f_code.co_name
    logger.error(caller + "\n" + _msg(*args))

def debug(*args):
    logger.debug(_msg(*args))

def critical(*args):
    logger.critical(_msg(*args))
    raise Exception(msg)

def deprecated(msg='', hint=''):
    caller = sys._getframe(1).f_code.co_name
    logger.warning(caller + "\n" + msg + ", is deprecated, please use:" + hint + " instead.")

def renameKwarg(old, new, kwargs):
    if old in kwargs:
        logger.warning("Keyword argument name changed from '" + old + \
                 "' to '" + new + "'")
        kwargs[new] = kwargs.pop(old)

def warnNonEmptyArgs(kwargs):
    if len(kwargs) > 0:
        caller = sys._getframe(1).f_code.co_name
        logger.warning("Unrecognized keyword arguments for method: '" + caller
                       + "' "  + _msg(kwargs))

