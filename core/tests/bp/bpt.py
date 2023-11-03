import os
os.add_dll_directory("C:\\Users\\Carsten\\src\\gimli\\gimli\\core\\tests\\bp\\")
import _bpt_

print(dir(_bpt_))

try:
    _bpt_.t()
except Exception as e:
    print("catched:", e) # Exception will not be catched but throws the malloc problem

_bpt_.A()
_bpt_.A(2)