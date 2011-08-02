TRIANGLE=triangle.zip
TRIANGLE_URL=http://www.netlib.org/voronoi/$(TRIANGLE)

SUITESPARSE=SuiteSparse-3.6.1.tar.gz
SUITESPARSE_URL=http://www.cise.ufl.edu/research/sparse/SuiteSparse/$(SUITESPARSE)

default: prep triangle suitesparse

prep:
	mkdir -p lib

$(TRIANGLE): wget unzip
	wget -nc -nd $(TRIANGLE_URL)
#	rm -rf triangle
#	unzip -o -d triangle $(TRIANGLE)

triangle: prep $(TRIANGLE)
	echo "building triangle for ", ${OSTYPE}
 
	if [ "$(OSTYPE)" == "msys" -o "$(MSYSTEM)" == "MINGW32" ] ; then \
		sed -i -e 's/-DLINUX/-DCPU86/g' triangle/makefile ; \
	fi
	sed -i -e 's/CC = cc/CC = gcc -fPIC/g' triangle/makefile; \
	
	make -C triangle/ trilibrary
	rm -f lib/libtriangle.a
	ar cqs lib/libtriangle.a triangle/triangle.o

$(SUITESPARSE): wget unzip
	wget -nc -nd $(SUITESPARSE_URL)
	tar -xzvf $(SUITESPARSE)

suitesparse: prep $(SUITESPARSE)

	echo "building suitesparse for ", $(OSTYPE)

	if [ $(OSTYPE) == "msys" ]; then \
		sed -i -e 's/CC = cc/CC = gcc /g' SuiteSparse/UFconfig/UFconfig.mk; \
	else \
		sed -i -e 's/CC = cc/CC = gcc -fPIC/g' SuiteSparse/UFconfig/UFconfig.mk; \
	fi

	sed -i -e 's/CHOLMOD_CONFIG =/CHOLMOD_CONFIG = -DNPARTITION/g' SuiteSparse/UFconfig/UFconfig.mk; 

	for i in AMD CAMD COLAMD CCOLAMD LDL; do \
		$(MAKE) -C SuiteSparse/$$i clean; \
		$(MAKE) -C SuiteSparse/$$i; \
		cp SuiteSparse/$$i/Lib/*.a lib/; \
	done

	$(MAKE) -C SuiteSparse/CHOLMOD/Lib clean
	$(MAKE) -C SuiteSparse/CHOLMOD/Lib 
	cp SuiteSparse/CHOLMOD/Lib/*.a lib/; \

wget:
	if [ "$(OSTYPE)" == "msys" -o "$(MSYSTEM)" == "MINGW32" ] ; then \
		if ( wget --version ); then \
			echo "########### wget found: ok" ; \
		else \
			echo "########### Installing wget" ; \
			mingw-get.exe install msys-wget ; \
		fi; \
	fi
	
unzip:
	if [ "$(OSTYPE)" == "msys" -o "$(MSYSTEM)" == "MINGW32" ] ; then \
		if ( unzip -v ); then \
			echo "########### unzip found: ok" ; \
		else \
			echo "########### Installing unzip" ; \
			mingw-get.exe install msys-unzip ; \
		fi; \
	fi

clean:
	echo "Cleaning ..."
	rm -rf lib triangle AMD CAMD COLAMD CCOLAMD LDL CHOLMOD SuiteSparse *~

	
