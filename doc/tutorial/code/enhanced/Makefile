GIMLIINC=-I ../../../../src
GIMLILIB=-L ../../../../src/.libs -lgimli -L../../../../external/lib -lcolamd
#GIMLILIB=-L$(HOME)/lib -lgimli

default: mt1dinv0 ttinvoffset

mt1dinv0: mt1dinv0.cpp
	g++ -g -O3 mt1dinv0.cpp -o mt1dinv0 $(GIMLIINC) $(GIMLILIB)

ttinvoffset: ttinvoffset.cpp
	g++ -g -O3 ttinvoffset.cpp -o ttinvoffset $(GIMLIINC) $(GIMLILIB)

clean:
	rm -rf *~ mt1dinv0 ttinvoffset
