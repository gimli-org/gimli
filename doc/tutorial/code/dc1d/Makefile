GIMLIINC=-I ../../../../../trunk/src
GIMLILIB=-L$(HOME)/lib -lgimli

default: dc1dsmooth dc1dblock dc1dsmooth-struct dc1dsmooth-region

dc1dsmooth: dc1dsmooth.cpp
	g++ -g -O3 dc1dsmooth.cpp -o dc1dsmooth $(GIMLIINC) $(GIMLILIB)
	
dc1dblock: dc1dblock.cpp
	g++ -g -O3 dc1dblock.cpp -o dc1dblock $(GIMLIINC) $(GIMLILIB)

dc1dsmooth-struct: dc1dsmooth-struct.cpp
	g++ -g -O3 dc1dsmooth-struct.cpp -o dc1dsmooth-struct $(GIMLIINC) $(GIMLILIB)

dc1dsmooth-region: dc1dsmooth-region.cpp
	g++ -g -O3 dc1dsmooth-region.cpp -o dc1dsmooth-region $(GIMLIINC) $(GIMLILIB)

clean:
	rm -rf *~ dc1dsmooth dc1dblock dc1dsmooth-struct dc1dsmooth-region
	
check:
	make clean default
	./dc1dblock sond1-100.dat