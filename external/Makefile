TRIANGLE=triangle.zip
TRIANGLE_URL=http://www.netlib.org/voronoi/$(TRIANGLE)

SUITESPARSE=SuiteSparse-3.6.1.tar.gz
SUITESPARSE_URL=http://www.cise.ufl.edu/research/sparse/SuiteSparse/$(SUITESPARSE)

default: prep triangle suitesparse

prep:
	mkdir -p lib

$(TRIANGLE): wget unzip
	#wget -nc -nd $(TRIANGLE_URL)
	#unzip -o -d triangle $(TRIANGLE)

triangle: prep $(TRIANGLE)
	echo "building triangle for ", ${OSTYPE}
 
	if [ "$(OSTYPE)" == "msys" -o "$(MSYSTEM)" == "MINGW32" ] ; then \
		sed -i -e 's/-DLINUX/-DCPU86/g' triangle/makefile ; \
	fi
	sed -i -e 's/CC = cc/CC = gcc -fPIC/g' triangle/makefile; \
	
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

	