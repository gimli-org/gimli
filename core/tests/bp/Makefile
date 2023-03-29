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
CXX=g++ 
AR=g++ 

### leads to undefined symbol in flat namespace for a bp object
#CXX := /opt/homebrew/bin/g++-12 -stdlib=libstdc++ # compiles but undef. symbol
#CXX := /opt/homebrew/bin/g++-12 -stdlib=libc++ # does not compile
#AR := /opt/homebrew/bin/g++-12  
 
CFLAGS := -fPIC -Wno-deprecated-declarations 
#CFLAGS := $(CFLAGS) -DBOOST_PYTHON_STATIC_LIB  
LFLAGS := -bundle -undefined dynamic_lookup 


#######################
## LINUX variant
#######################
ENV := /usr
PYINC := -I $(ENV)/include/python3.10/
PYLIB := $(ENV)/lib64/libpython3.10.so
BPINC := -I $(ENV)/include
BPLIB := $(ENV)/lib64/libboost_python310.so $(PYLIB)
CXX := g++ 
AR := g++ 
LFLAGS := -shared


SRC := $(wildcard ./src/generated/*.cpp)
OBJ := $(patsubst ./src/generated/%.cpp,./src/generated/%.o,$(SRC))

default:
	$(CXX) $(CFLAGS) -c bpt.cpp $(PYINC) $(BPINC)
	$(AR) -o $(TARGET).so *.o $(BPLIB) $(LFLAGS)

%.o: %.cpp
	$(CXX) -fPIC $(PYINC) $(BPINC) -I ./src -c $< -o $@

pgtest: $(OBJ)
	$(AR) -o pgtest.dylib $(OBJ) $(BPLIB)

test:
	python bpt.py

clean:
	rm -f *.o *.so *~ *.dylib $(OBJ)
