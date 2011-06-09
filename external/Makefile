TRIANGLE=triangle.zip
TRIANGLE_URL=http://www.netlib.org/voronoi/$(TRIANGLE)

SUITESPARSE=SuiteSparse-3.6.1.tar.gz
SUITESPARSE_URL=http://www.cise.ufl.edu/research/sparse/SuiteSparse/$(SUITESPARSE)


default: prep triangle suitesparse

prep:
	mkdir -p lib

$(TRIANGLE):
	wget -nc -nd $(TRIANGLE_URL)
	unzip -d triangle $(TRIANGLE)

triangle: prep $(TRIANGLE)
	echo "building triangle for ", ${OSTYPE}
 
	if [ $(OSTYPE) == "msys" ]; then \
		sed -i -e 's/-DLINUX/-DCPU86/g' triangle/makefile && \
		sed -i -e 's/CC = cc/CC = gcc -fPIC/g' triangle/makefile; \
	else \
		sed -i -e 's/CC = cc/CC = gcc -fPIC/g' triangle/makefile; \
	fi
	
	make -C triangle/ trilibrary
	rm -f lib/libtriangle.a
	ar cqs lib/libtriangle.a triangle/triangle.o

$(SUITESPARSE):
	wget -nc -nd $(SUITESPARSE_URL)
	tar -xzvf $(SUITESPARSE)

suitesparse: prep $(SUITESPARSE)

	echo "building suitesparse for ", $(OSTYPE)

	ln -fs SuiteSparse/LDL
	ln -fs SuiteSparse/AMD
	ln -fs SuiteSparse/CAMD
	ln -fs SuiteSparse/COLAMD
	ln -fs SuiteSparse/CCOLAMD
	ln -fs SuiteSparse/CHOLMOD

	if [ $(OSTYPE) == "msys" ]; then \
		sed -i -e 's/CC = cc/CC = gcc /g' SuiteSparse/UFconfig/UFconfig.mk; \
	else \
		sed -i -e 's/CC = cc/CC = gcc -fPIC/g' SuiteSparse/UFconfig/UFconfig.mk; \
	fi

	sed -i -e 's/CHOLMOD_CONFIG =/CHOLMOD_CONFIG = -DNPARTITION/g' SuiteSparse/UFconfig/UFconfig.mk; 

	for i in AMD CAMD COLAMD CCOLAMD LDL CHOLMOD; do \
		$(MAKE) -C $$i/Lib clean; \
		$(MAKE) -C $$i/Lib; \
		cp SuiteSparse/$$i/Lib/*.a lib/; \
	done

clean:
	echo "Cleaning ..."
	rm -rf lib triangle AMD CAMD COLAMD CCOLAMD LDL CHOLMOD SuiteSparse *~
