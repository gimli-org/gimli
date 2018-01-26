import logging

import pygimli as pg

#log = logging.getLogger('pyGIMLi')

#logging.basicConfig(level=logging.DEBUG,
                    #format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    #datefmt='%m/%d/%Y %H:%M:%S',
                    ##filename='example.log'
                    #)
pg.version()

# test pygimli log
pg.info("Start numeric log test." + str(pg.log(pg.RVector(1, 1.))))
pg.warn("Start warning test.")
pg.error("Start error test.")
#pg.critical("Start critical test.")

pg.debug("debug 0")

pg.setDebug(1)
pg.debug("debug ON")
pg.setThreadCount(2)

# should not printed out
pg.setDebug(0)
pg.debug("debug OFF")
pg.setThreadCount(2)

# test core log (should not be used outside the core)
pg.log(pg.Info, "core log ")
pg.log(pg.Warning, "core log ")
pg.log(pg.Error, "core log ")
pg.log(pg.Critical, "core log ")

#pg.logger.exception("Exception")

