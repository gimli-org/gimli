import _bpt_

print(dir(_bpt_))

try:
    _bpt_.t()
except Exception as e:
    print("catched:", e) # Exception will not be catched but throws the malloc problem
