TARGETNAME=_bpt_
## LINUX settings
#
#
#

### MACOS settings
# PYINC=-I $(ENV)/include/python3.6m/
# BPLIB=$(ENV)/lib/libboost_python36.dylib
# BPINC=-I $(ENV)/include
# PYLIB=-undefined dynamic_lookup
# CXX=clang++ 
# AR=clang++ -bundle


#cd C:/Users/carsten/src/gimli/build/core/python && C:/Users/carsten/msys64/mingw64/bin/cmake.exe -E rm -f CMakeFiles/_pygimli_.dir/objects.a
#cd C:/Users/carsten/src/gimli/build/core/python && C:/Users/carsten/msys64/mingw64/bin/ar.exe cr CMakeFiles/_pygimli_.dir/objects.a @CMakeFiles/_pygimli_.dir/objects1.rsp
#cd C:/Users/carsten/src/gimli/build/core/python && C:/Users/carsten/msys64/mingw64/bin/c++.exe -O3 -DNDEBUG -Wno-deprecated-declarations -Wno-class-memaccess -Wno-class-memaccess -Wno-parentheses -Wno-attributes -Wno-deprecated-declarations -Wno-parentheses -Wno-attributes -shared -o C:/Users/carsten/src/gimli/gimli/pygimli/core/_pygimli_.pyd -Wl,--major-image-version,0,--minor-image-version,0 -Wl,--whole-archive CMakeFiles/_pygimli_.dir/objects.a -Wl,--no-whole-archive @CMakeFiles/_pygimli_.dir/linklibs.rsp
# C:/Users/carsten/msys64/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/10.2.0/../../../../x86_64-w64-mingw32/bin/ld.exe: CMakeFiles/_pygimli_.dir/objects.a(type.pypp.cpp.obj):type.pypp.cpp:(.text.startup+0x7): undefined reference to `__imp__Py_NoneStruct'
# collect2.exe: error: ld returned 1 exit status

PYV=38
ANACONDA=../../../../../../anaconda3_2021_11

GCCV=11.2.0

#add_definitions(/DBOOST_PYTHON_STATIC_LIB)

### WIN32/MINGW64 settings
TARGET=$(TARGETNAME).pyd
#BOOST_ROOT=../../../../thirdParty/dist-GNU-$(GCCV)-64/boost_1_68_0-gcc-$(GCCV)-64-py$(PYV)
PYTHON_ROOT=$(ANACONDA)/envs/py$(PYV)
BOOST_ROOT=$(ANACONDA)/envs/py$(PYV)/Library

PYINC=-I $(PYTHON_ROOT)/include
PYLIB=$(PYTHON_ROOT)/libs/python$(PYV).lib

BPINC=-I $(BOOST_ROOT)/include
#BPAR=$(BOOST_ROOT)/lib/libboost_python$(PYV)-mt.dll
#BPLIB=$(BPAR).a
#BPLIB=$(BOOST_ROOT)/bin/boost_python39.dll
#BPLIB=$(BOOST_ROOT)/lib/boost_python39.lib
#BPLIB=$(BOOST_ROOT)/lib/boost_python39-vc141-mt-x64-1_73.lib
#BPLIB=$(BOOST_ROOT)/lib/libboost_python39.lib
#BPLIB=$(BOOST_ROOT)/lib/libboost_python39.lib
#BPLIB=-L . -lboost_python39

#BPCFLAGS=-D BOOST_PYTHON_STATIC_LIB
#BPLIB=./libboost_python39.lib

BPCFLAGS=
BPLIB=./boost_python39.dll

#pacman -S mingw-w64-x86_64-tools
#gendef /c/Users/Carsten/anaconda3_2021_11/envs/py39/Library/bin/boost_python39.dll
#dlltool.exe -d boost_python39.def -D /c/Users/Carsten/anaconda3_2021_11/envs/py39/Library/bin/boost_python39.dll -l libboost_python39.a
#BPLIB=libboost_python39.a

GIMLIINC=-I ../../src
GIMLILIB=-L ../../../../build/lib/ -lgimli


CXX=g++ -Wno-attributes -Wno-deprecated-declarations
AR=g++ -s -shared 

default: test

#c++.exe -DNDEBUG -shared -o _bpt_.pyd -Wl,--whole-archive _bpt_.a -Wl,--no-whole-archive $(PYLIB) $(BPLIB)
_bpt_.pyd: bpt.o
	ar.exe r _bpt_.a bpt.o
	g++.exe -DNDEBUG -shared -o _bpt_.pyd -Wl,--whole-archive _bpt_.a -Wl,--no-whole-archive $(PYLIB) $(BPLIB)
	
bpt.o: bpt.cpp
	$(CXX) -fPIC $(BPCFLAGS) -c bpt.cpp $(PYINC) $(BPINC)

_cpt_.pyd: cpt.o
	ar.exe r _cpt_.a cpt.o
	g++.exe -DNDEBUG -shared -o _cpt_.pyd -Wl,--whole-archive _cpt_.a -Wl,--no-whole-archive $(PYLIB)
	g++.exe -DNDEBUG -o cpt.exe -Wl,--whole-archive cpt.o -Wl,--no-whole-archive $(PYLIB)

cpt.o: cpt.cpp
	$(CXX) -fPIC $(BPCFLAGS) -c cpt.cpp $(PYINC)


# $(TARGET): bpt_gimli.o
# 	ar.exe r _bpt_.a bpt.o
# 	g++.exe -DNDEBUG -shared -o _bpt_.pyd -Wl,--whole-archive _bpt_.a -Wl,--no-whole-archive $(PYLIB) -Wl,-Bstatic $(BPLIB)

	   		
# ar.exe r _bpt__gimli.a bpt_gimli.o
# c++.exe -DNDEBUG -shared -o _bpt__gimli.pyd -Wl,--whole-archive _bpt__gimli.a -Wl,--no-whole-archive $(BPLIB) $(PYLIB) $(GIMLILIB)
# $(TARGET): bpt.o 
# 	$(AR) -o $(TARGET) *.o $(BPLIB) $(PYLIB)


bpt_gimli.o: bpt_gimli.cpp
	$(CXX) -fPIC -c bpt_gimli.cpp $(PYINC) $(BPINC) $(GIMLIINC)


#cp $(BPAR) .
test: $(TARGET)
	python -u bpt.py
	python -u bpt_gimli.py

clean:
	rm -f *.a *.pyd *.o *.so *~ *.dll
