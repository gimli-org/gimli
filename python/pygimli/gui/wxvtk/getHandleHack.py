# -*- coding: iso-8859-1 -*-
# Copyright (C) 2007 Maryland Robotics Club 
# Copyright (C) 2007 Joseph Lisee <jlisee@umd.edu> 
# All rights reserved. 
# 
# Author: Joseph Lisee <jlisee@umd.edu> 
#
# modified for pybert's vtk import to fit python-2.9@linux
# Standard python imports 
import ctypes as __ctypes 
import os.path as __path 

# Library Imports 
import wx as __wx 

# Load library 
__lib_path = __path.join(__path.dirname( __file__ ), 'libgetHandleHack.so') 
__lib = __ctypes.cdll.LoadLibrary(__path.abspath(__lib_path)) 
# Look up function 
__get_window_handle = __lib.get_window_handle 
# Set types to match prototype: char* get_window_handle_str(void*) 
__get_window_handle.restype = __ctypes.c_long 
__get_window_handle.argtypes = [__ctypes.c_void_p] 

def get_window_handle_str(window): 
    if not isinstance(window, __wx.Window): 
        raise TypeError('Must be called with a wx.Window or a subclass') 
    
    handle = __get_window_handle( int(window.this) )
    if handle == -2:
        print("not window found")
        handle = 0
    elif handle == -1:
        print("not widget found")
        handle = 0
        
    return str( handle )