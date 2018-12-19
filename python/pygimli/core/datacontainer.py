# -*- coding: utf-8 -*-
"""
Extensions to the core DataContainer class[es].
"""

from .._logger import deprecated, info, warn, critical
from ._pygimli_ import (RVector3, RVector, DataContainer, DataContainerERT)


def __DataContainerERT_setSensors(self, sensors):
    """Set Sensor positions.

    Set all sensor positions. 
    This is just syntactic sugar for setSensorPositions.

    Parameters
    ----------
    sensors: iterable
        Iterable that can be converted into a pg.Pos.
    
    Tests
    -----
    >>> import pygimli as pg
    >>> d = pg.DataContainerERT()
    >>> d.setSensors(pg.utils.grange(0.0, 3, n=4))
    >>> assert d.sensorCount() == 4
    """
    for s in sensors:
        if isinstance(s, float) or isinstance(s, int):
            self.createSensor(RVector3(s, 0.0))
        else:
            self.createSensor(RVector3(s))
    
DataContainer.setSensors = __DataContainerERT_setSensors

def __DataContainerERT_addFourPointData(self, *args):
    """Add a new data point to the end of the dataContainer.
    
    Add a new 4 point measurement to the end of the dataContainer and increase
    the data size by one. The index of the new data point is returned.
     
    Parameters
    ----------
    *args: [int]
        At least for index values for A, B, M and N.

    Returns
    -------
    ret: int
        Index of this new data point.
    
    Examples
    --------
    >>> import pygimli as pg
    >>> d = pg.DataContainerERT()
    >>> d.setSensors(pg.utils.grange(0, 3, n=4))
    >>> d.addFourPointData(0,1,2,3)
    0
    >>> d.addFourPointData([3,2,1,0])
    1
    >>> print(d)
    Data: Sensors: 4 data: 2
    """
    try:
        if len(args) == 1:
            return self.createFourPointData(self.size(), 
                                            args[0][0], args[0][1], 
                                            args[0][2], args[0][3])
        else:
            return self.createFourPointData(self.size(), 
                                            args[0], args[1], 
                                            args[2], args[3])
    except:
        print("args:", args)
        critical("Can't interpret arguments:", *args)
    
DataContainerERT.addFourPointData = __DataContainerERT_addFourPointData
