#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Caching manager with function decorator.

Does not have many functions yet but shows the way for more.

Input supports python base types and all pg.core objects with .hash() method.
Output supports DataContainerERT, ...

TODO:

  *   Output types:
        numpy.ndarray, pg.Mesh. pg.Vector, pg.Matrix


To use just add the decorator.

@pg.cache
def myLongRunningStuff(*args, **kwargs):
    #...
    return results
"""
import os 
import hashlib
import sys
import json
import time

import pygimli as pg

class Cache(object):
    def __init__(self, hashValue):
        self._value = None
        self._hash = hashValue
        self._name = CacheManager().cachingPath(str(self._hash))
        self._info = None
        self.restore()

    @property
    def info(self):
        if self._info is None:
            self._info = {'type': '',
                          'file': '',
                          'date': 0,
                          'dur': 0.0,
                          'restored': 0,
                          'codeinfo': '',
                          'args': '',
                          'kwargs': {},
                          }
        return self._info

    @info.setter
    def info(self, i):
        self._info = i

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, v):
        self.info['type'] = str(type(v).__name__)
        
        # if len(self.info['type']) != 1:
        #     pg.error('only single return caches supported for now.')    
        #     return 

        self.info['file'] = self._name
        
        self.updateCacheInfo()

        v.save(self._name)
        self._value = v
        pg.info('Cache stored:', self._name)

    def updateCacheInfo(self):
        with open(self._name + '.json', 'w') as of:
            json.dump(self.info, of, sort_keys=False,
                      indent=4, separators=(',', ': '))   

    def restore(self):
        if os.path.exists(self._name):
            
            try:
                with open(self._name + '.json', 'r') as file:
                    self.info = json.load(file)
                
                # if len(self.info['type']) != 1:
                #     pg.error('only single return caches supported for now.')    

                if self.info['type'] == 'DataContainerERT':
                    self._value = pg.DataContainerERT(self.info['file'])
                    
                if self.value is not None:
                    self.info['restored'] = self.info['restored'] + 1
                    self.updateCacheInfo()
                    pg.info('Cache restored ({1}s x {0}): {2}'.format(self.info['restored'], 
                                                                    round(self.info['dur'], 1),
                                                                    self._name))
            except Exception as e:
                print(e)
                pg.error('Cache restoring failed.')

class CacheManager(object):
    __instance = None
    __has_init = False

    def __new__(cls):
        if cls.__instance is None:
            cls.__instance = object.__new__(cls)
        return cls.__instance

    def __init__(self):
        if not self.__has_init:
            self._caches = {}
            self.__has_init = True

    @staticmethod
    def instance(cls):
        return cls.__instance__
        
    def cachingPath(self, fName):
        """Create a path name for the cache"""        
        if not os.path.exists('.cache'):
            os.mkdir('.cache')
        return os.path.join('.cache/', fName)

    def functInfo(self, funct):
        """Return unique info string about the called function."""
        return funct.__code__.co_filename + ":" + funct.__qualname__ + ":" + pg.versionStr()

    def hash(self, funct, *args, **kwargs):
        """"Create a hash value"""
        functInfo = self.functInfo(funct)
        funcHash = int(hashlib.sha224(functInfo.encode()).hexdigest()[:16], 16)
        codeHash = int(hashlib.sha224(funct.__code__.co_code).hexdigest()[:16], 16)

        argHash = 0
        for a in args:
            argHash = argHash ^ hash(a)
        
        for k, v in kwargs.items():
            argHash = argHash ^ hash(v)
                
        return funcHash ^ codeHash ^ argHash

    def cache(self, funct, *args, **kwargs):
        """ Create a unique cache """
        hashVal = self.hash(funct, *args, **kwargs)

        cache = Cache(hashVal)
        cache.info['codeinfo'] = self.functInfo(funct)
        cache.info['args'] = str(args)
        cache.info['kwargs'] = str(kwargs)

        return cache


def cache(funct):
    """Cache decorator."""
    pg.tic()
    def wrapper(*args, **kwargs):
        
        cache = CacheManager().cache(funct, *args, **kwargs)
        # pg.toc('restore')
        if cache.value is not None:
            return cache.value
        else:
            rv = funct(*args, **kwargs)
            cache.info['date'] = time.time()
            cache.info['dur'] = pg.dur()
            cache.value = rv
            # pg.toc('store')
            return rv
    # pg.toc('cache')
    return wrapper
