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

logger = logging.getLogger('pyGIMLi')

def setDebug(d):
    if d:
        core._pygimli_.setDebug(True)
        logger.setLevel(logging.DEBUG)
        logging.getLogger('Core').setLevel(logging.DEBUG)
        logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%m/%d/%Y %H:%M:%S',
                    #filename='pygimli.log'
                    )
        logger.debug("Set debug mode: on")
    else:
        core._pygimli_.setDebug(False)
        logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%m/%d/%Y %H:%M:%S',
                    #filename='pygimli.log'
                    )
        logger.debug("Set debug mode: off")
        logging.getLogger('Core').setLevel(logging.INFO)
        logger.setLevel(logging.INFO)

if '--debug' in sys.argv:
    setDebug(True)
else:
    setDebug(False)

def _msg(*args):
    msg = ''
    for i, arg in enumerate(args):
        msg += str(arg)
        if i < len(args)-1:
            msg += ' '
    return msg

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

