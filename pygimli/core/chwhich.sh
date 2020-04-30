#!/usr/bin/env bash
PATH=.:$PATH  # make sure it takes the one here which is installed
for dll in *.dll *.pyd;do 
objdump.exe -p $dll |grep DLL
done|sort|uniq|tail -n +3|
while read eins zwei dll;do 
which $dll
done|sort
