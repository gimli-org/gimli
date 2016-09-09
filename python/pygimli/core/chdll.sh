#!/usr/bin/env bash
for dll in *.dll *.pyd;do objdump.exe -p $dll |grep DLL;done|sort|uniq
