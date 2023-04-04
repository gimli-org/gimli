TARGET := _bpt_

#######################
## MACOS variant
#######################

ENV := /opt/homebrew/Caskroom/mambaforge/base
PYINC := -I $(ENV)/include/python3.9/
PYLIB := $(ENV)/lib/libpython3.9.dylib
BPINC := -I $(ENV)/include
#BPLIB=$(ENV)/lib/libboost_python39.a
BPLIB=$(ENV)/lib/libboost_python39.dylib

## compiles and runs
#CXX=clang++  
#AR=clang++ 

## compiles and runs
CXX=c++ 
AR=c++ 

### leads to undefined symbol in flat namespace for a bp object
#CXX := /opt/homebrew/bin/g++-12 -stdlib=libstdc++ # compiles but undef. symbol
#CXX := /opt/homebrew/bin/g++-12 -stdlib=libc++ # does not compile
#AR := /opt/homebrew/bin/g++-12  
CFLAGS := -fPIC -Wno-deprecated-declarations -std=c++11
#CFLAGS := $(CFLAGS) -DBOOST_PYTHON_STATIC_LIB  
LFLAGS := -bundle -undefined dynamic_lookup 


#######################
## LINUX variant
#######################
#ENV := /usr
#PYINC := -I $(ENV)/include/python3.10/
#PYLIB := $(ENV)/lib64/libpython3.10.so
#BPINC := -I $(ENV)/include
#BPLIB := $(ENV)/lib64/libboost_python310.so $(PYLIB)
#CXX := g++ 
#AR := g++ 
#LFLAGS := -shared

GINC := -I ../../src -I ./src/generated -I $(ENV)/lib/python3.9/site-packages/numpy/core/include/ -Wno-undefined-var-template
GINC := $(GINC) -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX13.1.sdk
GINC := $(GINC) -isystem /opt/homebrew/Caskroom/mambaforge/base/include
GINC := $(GINC) -MD -MT
GINC := $(GINC) -DBOOST_PYTHON_DYN_LINK -DBOOST_PYTHON_NO_LIB -DBOOST_PYTHON_NO_PY_SIGNATURES -DBOOST_PYTHON_USE_GCC_SYMBOL_VISIBILITY -DPYGIMLI -D_pygimli__EXPORTS 

CXX := /Library/Developer/CommandLineTools/usr/bin/c++ 
GINC := -DBOOST_PYTHON_DYN_LINK -DBOOST_PYTHON_NO_LIB -DBOOST_PYTHON_NO_PY_SIGNATURES -DBOOST_PYTHON_USE_GCC_SYMBOL_VISIBILITY -DHAVE_CONFIG_CMAKE_H -DPYGIMLI -D_pygimli__EXPORTS -I/Users/fwagner/git/gimli/build/core/python/generated -I/Users/fwagner/git/gimli/source/core/src -I/Users/fwagner/git/gimli/build -I/opt/homebrew/Caskroom/mambaforge/base/include/python3.9 -I/Users/fwagner/git/gimli/build/core/python -I/opt/homebrew/Caskroom/mambaforge/base/lib/python3.9/site-packages/numpy/core/include -isystem /opt/homebrew/Caskroom/mambaforge/base/include -Wno-comment -Wno-misleading-indentation -Wno-class-memaccess -Wno-vla -Wno-maybe-uninitialized -Wno-unused-variable -Wno-unused-but-set-variable -Wno-int-in-bool-context -Wno-unused-local-typedefs -Wno-deprecated-declarations -Wno-unknown-warning-option -Wno-undefined-var-template -DBOOST_BIND_GLOBAL_PLACEHOLDERS -Wno-infinite-recursion -Wno-infinite-recursion -O3 -DNDEBUG -O2 -ansi -pedantic -fno-omit-frame-pointer -Wall -Wno-overloaded-virtual -Qunused-arguments -std=gnu++11 -arch arm64 -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX13.1.sdk -fPIC -MD -MT -MF 

GLIB := libgimli.dylib

#GLIB :=  -Wno-comment -Wno-misleading-indentation -Wno-class-memaccess -Wno-vla -Wno-maybe-uninitialized -Wno-unused-variable -Wno-unused-but-set-variable -Wno-int-in-bool-context -Wno-unused-local-typedefs -Wno-deprecated-declarations -Wno-unknown-warning-option -Wno-undefined-var-template -DBOOST_BIND_GLOBAL_PLACEHOLDERS -Wno-infinite-recursion -Wno-infinite-recursion -O3 -DNDEBUG -O2 -ansi -pedantic -fno-omit-frame-pointer -Wall -Wno-overloaded-virtual -Qunused-arguments -arch arm64 -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX13.1.sdk -bundle -Wl,-headerpad_max_install_names
GLIB := 

LFLAGS := libgimli.dylib -bundle -undefined dynamic_lookup $(PYLIB)

#LFLAGS := -Wl,-rpath,/Users/fwagner/git/gimli/build/lib -Wl,-rpath,/opt/homebrew/Caskroom/mambaforge/base/lib libgimli.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libpython3.9.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libboost_python39.dylib -bundle -undefined dynamic_lookup -bundle -undefined dynamic_lookup /Users/fwagner/git/gimli/thirdParty/dist-AppleClang-14.0.0.14000029-64/lib/libtriangle.a /opt/homebrew/Caskroom/mambaforge/base/lib/libcholmod.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libcamd.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libamd.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libcolamd.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libccolamd.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libsuitesparseconfig.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libopenblas.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libumfpack.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libamd.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libcolamd.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libccolamd.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libsuitesparseconfig.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libopenblas.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libumfpack.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libpython3.9.dylib 

SRC := $(wildcard ./src/generated/*.cpp)
OBJ := $(patsubst ./src/generated/%.cpp,./src/generated/%.o,$(SRC))

default:
	$(CXX) $(CFLAGS) -c bpt.cpp $(PYINC) $(BPINC)
	$(AR) -o $(TARGET).so *.o $(BPLIB) $(LFLAGS)

%.o: %.cpp
	$(CXX) $(CFLAGS) $(PYINC) $(GINC) $(BPINC) -I ./src -c $< -o $@

pgtest: $(OBJ)
	$(AR) -o _pygimli_.so $(OBJ) $(BPLIB) $(GLIB) $(LFLAGS)

test:
	python bpt.py

clean:
	rm -f *.o *.so *~ *.dylib $(OBJ)

#/Library/Developer/CommandLineTools/usr/bin/c++ -DBOOST_PYTHON_DYN_LINK -DBOOST_PYTHON_NO_LIB -DBOOST_PYTHON_NO_PY_SIGNATURES -DBOOST_PYTHON_USE_GCC_SYMBOL_VISIBILITY -DHAVE_CONFIG_CMAKE_H -DPYGIMLI -D_pygimli__EXPORTS 
#-I/Users/fwagner/git/gimli/build/core/python/generated 
#-I/Users/fwagner/git/gimli/source/core/src 
#-I/Users/fwagner/git/gimli/build 
#-I/opt/homebrew/Caskroom/mambaforge/base/include/python3.9 
#-I/Users/fwagner/git/gimli/build/core/python 
#-I/opt/homebrew/Caskroom/mambaforge/base/lib/python3.9/site-packages/numpy/core/include 
#-isystem /opt/homebrew/Caskroom/mambaforge/base/include 
#-Wno-comment -Wno-misleading-indentation -Wno-class-memaccess -Wno-vla -Wno-maybe-uninitialized -Wno-unused-variable -Wno-unused-but-set-variable -Wno-int-in-bool-context -Wno-unused-local-typedefs -Wno-deprecated-declarations -Wno-unknown-warning-option -Wno-undefined-var-template -DBOOST_BIND_GLOBAL_PLACEHOLDERS -Wno-infinite-recursion -Wno-infinite-recursion 
#-O3 -DNDEBUG -O2 -ansi -pedantic -fno-omit-frame-pointer -Wall -Wno-overloaded-virtual -Qunused-arguments 
#-std=gnu++11 -arch arm64 
#-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX13.1.sdk -fPIC -MD -MT 
#core/python/CMakeFiles/_pygimli_.dir/generated/BVector.pypp.cpp.o -MF CMakeFiles/_pygimli_.dir/generated/BVector.pypp.cpp.o.d -o CMakeFiles/_pygimli_.dir/generated/BVector.pypp.cpp.o -c /Users/fwagner/git/gimli/build/core/python/generated/BVector.pypp.cpp
#-Wl,-rpath,/Users/fwagner/git/gimli/build/lib -Wl,-rpath,/opt/homebrew/Caskroom/mambaforge/base/lib ../../lib/libgimli.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libpython3.9.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libboost_python39.dylib -bundle -undefined dynamic_lookup -bundle -undefined dynamic_lookup /Users/fwagner/git/gimli/thirdParty/dist-AppleClang-14.0.0.14000029-64/lib/libtriangle.a /opt/homebrew/Caskroom/mambaforge/base/lib/libcholmod.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libcamd.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libamd.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libcolamd.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libccolamd.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libsuitesparseconfig.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libopenblas.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libumfpack.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libamd.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libcolamd.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libccolamd.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libsuitesparseconfig.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libopenblas.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libumfpack.dylib /opt/homebrew/Caskroom/mambaforge/base/lib/libpython3.9.dylib 

#/Library/Developer/CommandLineTools/usr/bin/c++  -Wno-comment -Wno-misleading-indentation -Wno-class-memaccess -Wno-vla -Wno-maybe-uninitialized -Wno-unused-variable -Wno-unused-but-set-variable -Wno-int-in-bool-context -Wno-unused-local-typedefs -Wno-deprecated-declarations -Wno-unknown-warning-option -Wno-undefined-var-template -DBOOST_BIND_GLOBAL_PLACEHOLDERS -Wno-infinite-recursion -Wno-infinite-recursion -O3 -DNDEBUG -O2 -ansi -pedantic -fno-omit-frame-pointer -Wall -Wno-overloaded-virtual -Qunused-arguments -arch arm64 -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX13.1.sdk -bundle -Wl,-headerpad_max_install_names -o /Users/fwagner/git/gimli/source/pygimli/core/_pygimli_.so
