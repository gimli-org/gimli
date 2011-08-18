#!/bin/bash
MAKEOPTS=-j7
mkdir -p generated

[ -f pygimli.cache ] && rm pygimli.cache

if [ $OSTYPE = "msys" ]; then
	MAKEFILE=Makefile.msys
	MAKEOPTS=-j2
	echo "build for mingw"
else
	MAKEFILE=Makefile.linux
	echo "build for linux gcc"
fi
#	make -f $MAKEFILE clean

if [ $# -lt 1 ]; then 
	python generate_pygimli_code.py

#	[ -f generated/__call_policies.pypp.hpp ] || cp tmp/__call_policies.pypp.hpp generated/
	make -f $MAKEFILE $MAKEOPTS 
elif [ $# -gt 0 ]; then
        if [ "$1" = "test" ]; then
		echo "build testsuite"
		python generate_pygimli_code.py test
#		cp tmp/__call_policies.pypp.hpp generated/
		DEFINES='-D PYTEST' make -f $MAKEFILE $MAKEOPTS
		exit
	else
		echo "unknown command" $1	
	fi
fi

