# -*- coding: utf-8 -*-
"""
Extensions to the core DataContainer class[es].
"""
import numpy as np
from . logger import deprecated, info, warn, critical, verbose
from .core import (RVector3, RVector, IndexArray,
                   DataContainer, DataContainerERT)

def __DataContainer_str(self):
    return "Data: Sensors: " + str(self.sensorCount()) + " data: " + \
        str(self.size()) + ", nonzero entries: " + \
        str([d for d in self.dataMap().keys() if self.isSensorIndex(d) or
             self.haveData(d)])
DataContainer.__repr__ =__DataContainer_str
DataContainer.__str__ =__DataContainer_str


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

def __DC_setVal(self, key, val):
    if isinstance(key, (tuple, list)):
        assert len(key) == len(val), \
            "The number of key must match the number of items in the " + \
            "outer item list"
        for subkey, item in zip(key, val):
            self[subkey] = item
        return

    if len(val) >  self.size():
        verbose("DataContainer resized to:", len(val))
        self.resize(len(val))
    self.set(key, val)
DataContainer.__setitem__ = __DC_setVal

def __DC_getVal(self, key):
    if isinstance(key, (tuple, list)):
        return np.vstack([self[subkey] for subkey in key]).T
    if self.isSensorIndex(key):
        return np.array(self(key), dtype=int)
    #return self(key).array() // d['a'][2] = 0.0, would be impossible
    return self(key)
DataContainer.__getitem__ = __DC_getVal


def __DataContainerERT_addFourPointData(self, *args,
                                        indexAsSensors=False, **kwargs):
    """Add a new data point to the end of the dataContainer.

    Add a new 4 point measurement to the end of the dataContainer and increase
    the data size by one. The index of the new data point is returned.

    Parameters
    ----------
    *args: [int]
        At least four index values for A, B, M and N.
    indexAsSensors: bool [False]
        The indices A, B, M and N are additionally interpreted as sensor position in [m, 0, 0].
    **kwargs: dict
        Named values for the data configuration.

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
    Data: Sensors: 4 data: 2, nonzero entries: ['a', 'b', 'm', 'n', 'rhoa', 'valid']
    >>> print(d('rhoa'))
    2 [0.0, 1.0]
    """
    try:
        if len(args) == 1:
            a, b, m, n = args[0][:]
        else:
            [a, b, m, n] = args

        if indexAsSensors:
            a = self.createSensor([float(a), 0.0, 0.0])
            b = self.createSensor([float(b), 0.0, 0.0])
            m = self.createSensor([float(m), 0.0, 0.0])
            n = self.createSensor([float(n), 0.0, 0.0])
        idx = self.createFourPointData(self.size(), a, b, m, n)


    except Exception as e:
        print(e)
        print("args:", args, len(args))
        critical("Can't interpret arguments:", *args)

    for k, v in kwargs.items():
        if not self.haveData(k):
            self.add(k)
        self.ref(k)[idx] = v
    return idx

DataContainerERT.addFourPointData = __DataContainerERT_addFourPointData
