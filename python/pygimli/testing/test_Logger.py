import logging

import pygimli as pg

log = logging.getLogger('pyGIMLi')

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%m/%d/%Y %H:%M:%S',
                    #filename='example.log'
                    )

log.info("Start test.")
log.warn("Start test.")
log.debug("Start test.")
log.error("Start test.")
log.critical("Start test.")

# test core log
pg.setThreadCount(2)

# call core log
pg.log(pg.Debug, "core log")

log.exception("Exception")
