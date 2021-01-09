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


### WIN32/MINGW64 settings
TARGET=$(TARGETNAME).pyd
BOOST_ROOT=../../../../thirdParty/dist-GNU-10.2.0-64/boost_1_68_0-gcc-10.2.0-64-py38
PYTHON_ROOT=../../../../../../anaconda3

PYINC=-I $(PYTHON_ROOT)/include
PYLIB=$(PYTHON_ROOT)/libs/python38.lib
BPINC=-I $(BOOST_ROOT)/include
BPAR=$(BOOST_ROOT)/lib/libboost_python38-mt.dll
BPLIB=$(BPAR).a

CXX=g++ -Wno-attributes -Wno-deprecated-declarations
AR=g++ -s -shared 

default: test

$(TARGET): bpt.o 
	ar.exe r _bpt_.a bpt.o
	c++.exe -DNDEBUG -shared -o _bpt_.pyd -Wl,--whole-archive _bpt_.a -Wl,--no-whole-archive $(BPLIB) $(PYLIB)
# $(TARGET): bpt.o 
# 	$(AR) -o $(TARGET) *.o $(BPLIB) $(PYLIB)

bpt.o: bpt.cpp
	$(CXX) -fPIC -c bpt.cpp $(PYINC) $(BPINC)

test: $(TARGET)
	cp $(BPAR) .
	python -u bpt.py

clean:
	rm -f *.a *.pyd *.o *.so *~ *.dll
