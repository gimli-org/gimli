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

import pygimli as pg

class Cache(object):
    def __init__(self, hashValue):
        self._value = None
        self._hash = hashValue
        self._name = CacheManager().cachingPath(str(self._hash))
        self.restore()
        
    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, v):
        pack = {'type': str(type(v).__name__),
                'file': self._name }
        
        with open(self._name + '.json', 'w') as of:
            json.dump(pack, of)   
        
        v.save(self._name)
        self._value = v
        pg.info('Cache stored:', self._name)

    def restore(self):
        if os.path.exists(self._name):
            
            pg.info('Restoring cache:', self._name)

            try:
                with open(self._name + '.json', 'r') as file:
                    pack = json.load(file)
                
                if pack['type'] == 'DataContainerERT':
                    self._value = pg.DataContainerERT(pack['file'])
            except Exception as e:
                print(e)
                pg.error('Cache restoring failed.')

    def load(self, fileName):
        try:
            self._value = np.load(CacheManager().cachingPath(fileName))
        except:
            pass

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

    def hash(self, funct, *args, **kwargs):
        """"Create a hash value"""
        funcIdent = funct.__code__.co_filename + ":" + funct.__qualname__ 
        funcHash = int(hashlib.sha224(funcIdent.encode()).hexdigest()[:16], 16)
        codeHash = int(hashlib.sha224(funct.__code__.co_code).hexdigest()[:16], 16)

        argHash = 0
        for a in args:
            argHash = argHash ^ hash(a)
        
        for k, v in kwargs.items():
            pg.hashCombine(argHash, hash(v))
                
        return funcHash ^ codeHash ^ argHash

    def cache(self, funct, *args, **kwargs):
        """ Create a unique cache """
        hashVal = self.hash(funct, *args, **kwargs)

        cache = Cache(hashVal)
        return cache


def cache(funct):
    """Cache decorator."""
    def wrapper(*args, **kwargs):
        cache = CacheManager().cache(funct, *args, **kwargs)
        if cache.value:
            return cache.value
        else:
            rv = funct(*args, **kwargs)
            cache.value = rv
            return rv
    return wrapper
