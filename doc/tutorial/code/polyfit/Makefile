GIMLIINC=-I ../../../../../trunk/src
GIMLILIB=-L$(HOME)/lib -lgimli

default: polyfit0 polyfit1

polyfit0: polyfit0.cpp
	g++ -g -O3 polyfit0.cpp -o polyfit0 $(GIMLIINC) $(GIMLILIB)
	
polyfit1: polyfit1.cpp
	g++ -g -O3 polyfit1.cpp -o polyfit1 $(GIMLIINC) $(GIMLILIB)

clean:
	rm -rf *~ polyfit0 polyfit1
	