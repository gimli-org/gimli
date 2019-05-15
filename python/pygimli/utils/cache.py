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
import sys
import os 
import inspect
import hashlib
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
                          'version': '',
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

        if self.info['type'] == 'Mesh':
            v.saveBinaryV2(self._name)
        else:
            v.save(self._name)

        self._value = v
        pg.info('Cache stored:', self._name)

    def updateCacheInfo(self):
        with open(self._name + '.json', 'w') as of:
            json.dump(self.info, of, sort_keys=False,
                      indent=4, separators=(',', ': '))   

    def restore(self):
        if os.path.exists(self._name + '.json'):
            
            try:
                with open(self._name + '.json', 'r') as file:
                    self.info = json.load(file)
                
                # if len(self.info['type']) != 1:
                #     pg.error('only single return caches supported for now.')    

                if self.info['type'] == 'DataContainerERT':
                    self._value = pg.DataContainerERT(self.info['file'])
                elif self.info['type'] == 'RVector':
                    self._value = pg.RVector(self.info['file'])
                elif self.info['type'] == 'Mesh':
                    pg.tic()
                    self._value = pg.Mesh()
                    self._value.loadBinaryV2(self.info['file'] + '.bms')
                    pg.debug("Restoring cache took:", pg.dur(), "s")    
            
                if self.value is not None:
                    self.info['restored'] = self.info['restored'] + 1
                    self.updateCacheInfo()
                    pg.info('Cache {3} restored ({1}s x {0}): {2}'.format(self.info['restored'], 
                                                                    round(self.info['dur'], 1),
                                                                    self._name,
                                                                    self.info['codeinfo']))
                else:
                    pg.warn('Could not restore cache of type {0}.'.format(self.info['type']))
    
                pg.debug("Restoring cache took:", pg.dur(), "s")
            except Exception as e:
                import traceback
                traceback.print_exc(file=sys.stdout)
                print(self.info)
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
        return funct.__code__.co_filename + ":" + funct.__qualname__

    def strhash(self, string):
        return int(hashlib.sha224(string.encode()).hexdigest()[:16], 16)

    def hash(self, funct, *args, **kwargs):
        """"Create a hash value"""
        pg.tic()
        functInfo = self.functInfo(funct)
        funcHash = self.strhash(functInfo)
        versionHash = self.strhash(pg.versionStr())
        codeHash = self.strhash(inspect.getsource(funct))

        argHash = 0
        for a in args:
            if isinstance(a, str):
                argHash = argHash ^ self.strhash(a)
            elif isinstance(a, list):
                for item in a:
                    if isinstance(item, str):
                        argHash = argHash ^ self.strhash(item)
                    else:
                        argHash = argHash ^ hash(item)
            else:
                argHash = argHash ^ hash(a)
        
        for k, v in kwargs.items():
            if isinstance(v, str):
                argHash = argHash ^ self.strhash(v)
            else:
                argHash = argHash ^ hash(v)
                
        pg.debug("Hashing took:", pg.dur(), "s")
        return funcHash ^ versionHash ^ codeHash ^ argHash

    def cache(self, funct, *args, **kwargs):
        """ Create a unique cache """
        hashVal = self.hash(funct, *args, **kwargs)

        cache = Cache(hashVal)
        cache.info['codeinfo'] = self.functInfo(funct)
        cache.info['version'] = pg.versionStr()
        cache.info['args'] = str(args)
        cache.info['kwargs'] = str(kwargs)

        return cache


def cache(funct):
    """Cache decorator."""
    def wrapper(*args, **kwargs):
        cache = CacheManager().cache(funct, *args, **kwargs)
        if cache.value is not None:
            return cache.value
        else:
            rv = funct(*args, **kwargs)
            cache.info['date'] = time.time()
            cache.info['dur'] = pg.dur()
            try:
                cache.value = rv
            except Exception as e:
                print(e)
                pg.warn("Can't cache:", rv)
            return rv
    return wrapper
