# -*- coding: utf-8 -*-
"""
Extensions to the core DataContainer class[es].
"""

from .._logger import deprecated, info, warn, critical
from ._pygimli_ import (RVector3, RVector, DataContainer, DataContainerERT)


def __DataContainer_setSensors(self, sensors):
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
    for i, s in enumerate(sensors):
        nS = s
        if isinstance(s, float) or isinstance(s, int):
            nS = RVector3(s, 0.0)
            
        if i > self.sensorCount():
            self.createSensor(nS)
        else:
            self.setSensorPosition(i, nS)

    
DataContainer.setSensors = __DataContainer_setSensors


def __DataContainerERT_addFourPointData(self, *args, **kwargs):
    """Add a new data point to the end of the dataContainer.
    
    Add a new 4 point measurement to the end of the dataContainer and increase
    the data size by one. The index of the new data point is returned.
     
    Parameters
    ----------
    *args: [int]
        At least for index values for A, B, M and N.
    **args: dict
        Values for the actual data configuration.

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
    >>> d.addFourPointData([3,2,1,0], rhoa=1.0)
    1
    >>> print(d)
    Data: Sensors: 4 data: 2
    >>> print(d('rhoa'))
    2 [0.0, 1.0]
    """
    try:
        if len(args) == 1:
            idx =  self.createFourPointData(self.size(), 
                                            args[0][0], args[0][1], 
                                            args[0][2], args[0][3])
        else:
            idx = self.createFourPointData(self.size(), 
                                            args[0], args[1], 
                                            args[2], args[3])
    except:
        print("args:", args)
        critical("Can't interpret arguments:", *args)

    for k, v in kwargs.items():
        if not self.haveData(k):
            self.add(k)
        self.ref(k)[idx] = v
    return idx
    
DataContainerERT.addFourPointData = __DataContainerERT_addFourPointData
