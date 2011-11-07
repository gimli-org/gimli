GIMLIINC=-I ../../../../../trunk/src
GIMLILIB=-L../../../../../trunk/src/.libs -lgimli -lcolamd

default: dcem1dinv

dcem1dinv: dcem1dinv.cpp
	g++ -g -O3 dcem1dinv.cpp -o dcem1dinv $(GIMLIINC) $(GIMLILIB)

clean:
	rm -rf *~ dcem1dinv
	