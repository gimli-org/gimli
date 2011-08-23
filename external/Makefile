TRIANGLE=triangle.zip
TRIANGLE_URL=http://www.netlib.org/voronoi/$(TRIANGLE)

BLAS=lib/blas.dll
BLAS_SRC=blas.tgz
BLAS_URL=http://www.netlib.org/blas/$(BLAS_SRC)

LAPACK_SRC=lapack.tgz
LAPACK_URL=http://www.netlib.org/lapack/$(LAPACK_SRC)

SUITESPARSE=SuiteSparse-3.6.1.tar.gz
SUITESPARSE_URL=http://www.cise.ufl.edu/research/sparse/SuiteSparse/$(SUITESPARSE)


#default: prep triangle blas lapack suitesparse
default: prep triangle suitesparse

prep:
	mkdir -p lib

$(TRIANGLE): wget unzip
	wget -nc -nd $(TRIANGLE_URL)
	rm -rf triangle
	unzip -o -d triangle $(TRIANGLE)

triangle: prep $(TRIANGLE)
	echo "building triangle for ", ${OSTYPE}
 
	if [ "$(OSTYPE)" == "msys" -o "$(MSYSTEM)" == "MINGW32" ] ; then \
		sed -i -e 's/-DLINUX/-DCPU86/g' triangle/makefile ; \
	fi
	sed -i -e 's/CC = cc/CC = gcc -fPIC/g' triangle/makefile; \
	
	make -C triangle/ trilibrary
	rm -f lib/libtriangle.a
	ar cqs lib/libtriangle.a triangle/triangle.o

$(LAPACK_SRC): wget
	wget -nc -nd $(LAPACK_URL)
#	rm -rf lapack-*
#	tar -xzvf $(LAPACK_SRC)

BLASFILES = $(wildcard lapack-*/BLAS/SRC/*.f)
BLASOBJECTS = $(patsubst %.f,%.o,$(BLASFILES))

LAPACKFILES = $(wildcard lapack-*/src/*.f)
LAPACKOBJECTS = $(patsubst %.f,%.o,$(LAPACKFILES))

.f.o:
	gfortran -O3 -c -o $(patsubst %.f,%.o, $<) $<	
	
$(BLAS): prep $(LAPACK_SRC) $(BLASOBJECTS)
	echo "building blas for ", $(OSTYPE)
	gfortran -o $(BLAS) --strip-all -shared $(BLASOBJECTS)

lapack: prep $(LAPACK_SRC) $(BLAS) $(LAPACKOBJECTS)
	echo "building lapack for ", $(OSTYPE)
	gfortran -o $(LAPACK) --strip-all -shared $(LAPACKOBJECTS)
	
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

	
