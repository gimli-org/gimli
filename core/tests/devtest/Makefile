INC=-I ~/.local/include/eigen3/ -I ../../src

LIB=-L ../../../../build/lib/ -lgimli

all: test

test: test.cpp
	g++ test.cpp $(INC) -o test $(LIB) 