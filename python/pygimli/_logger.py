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

def info(msg):
    logger.info(msg)

def warn(msg):
    logger.warn(msg)

def error(msg):
    logger.error(msg)

def debug(msg):
    logger.debug(msg)

def critical(msg):
    logger.critical(msg)
    raise Exception(msg)

def deprecated(msg, hint):
    logger.warn(msg + ", is deprecated, use:" + hint + " instead.")

def renameKwarg(old, new, kwargs):
    if old in kwargs:
        logger.warn("Keyword argument name changed from '" + old + \
                 "' to '" + new + "'")
        kwargs[new] = kwargs.pop(old)

def warnNonEmptyArgs(kwargs):
    if len(kwargs) > 0:
        logger.warn("unrecognized keyword arguments", kwargs)

