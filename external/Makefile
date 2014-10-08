TRIANGLE=triangle.zip
TRIANGLE_URL=http://www.netlib.org/voronoi/$(TRIANGLE)

BLAS=lib/libblas.dll

LAPACK=lib/liblapack.dll
LAPACK_VER=3.4.2
LAPACK_SRC=lapack-$(LAPACK_VER).tgz
LAPACK_URL=http://www.netlib.org/lapack/$(LAPACK_SRC)

SUITESPARSE=SuiteSparse-4.3.1.tar.gz
SUITESPARSE_URL=http://www.cise.ufl.edu/research/sparse/SuiteSparse/$(SUITESPARSE)

#LDFLAGS="-static-libgfortran -static-libgcc"
#default: prep triangle blas lapack suitesparse
default: prep triangle suitesparse

prep:
	mkdir -p lib
	mkdir -p include

$(TRIANGLE): wget unzip patch
	wget -nc -nd $(TRIANGLE_URL)
	rm -rf triangle
	unzip -o -d triangle $(TRIANGLE)

triangle: prep $(TRIANGLE)
	echo "building triangle for ", ${OSTYPE}, $(MSYSTEM)
 
	if [ "$(OSTYPE)" == "msys" -o "$(MSYSTEM)" == "MINGW32" ] ; then \
		sed -i -e 's/-DLINUX/-DCPU86/g' triangle/makefile ; \
	fi
	sed -i -e 's/CC = cc/CC = gcc -fPIC/g' triangle/makefile; \
	
	patch triangle/triangle.c -i patches/triangle-mingw-win64.patch
	
	make -C triangle/ trilibrary
	rm -f lib/libtriangle.a
	ar cqs lib/libtriangle.a triangle/triangle.o

$(LAPACK_SRC): wget cmake
	wget -nc -nd $(LAPACK_URL)
	rm -rf lapack-$(LAPACK_VER)
	rm -rf lapack-build
	tar -xzvf $(LAPACK_SRC)

lapack: prep $(LAPACK_SRC)
	echo "building lapack for ", $(OSTYPE)
	mkdir -p lapack-build
	pushd lapack-build; \
	if [ "$(OSTYPE)" == "msys" -o "$(MSYSTEM)" == "MINGW32" ]; then \
		cmake -D BUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release -G 'MSYS Makefiles' ../lapack-$(LAPACK_VER); \
		$(MAKE); \
		cp bin/liblapack.dll ../$(LAPACK); \
		cp bin/libblas.dll ../$(BLAS); \
	else \
		cmake -D BUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release ../lapack-$(LAPACK_VER); \
		$(MAKE); \
		cp lib/liblapack.so ../lib/; \
		cp lib/libblas.so ../lib/; \
	fi; \
	popd

	
$(SUITESPARSE): wget unzip
	wget -nc -nd $(SUITESPARSE_URL)
	tar -xzvf $(SUITESPARSE)
		
suitesparse: prep $(SUITESPARSE)

	echo "building suitesparse for ", $(OSTYPE)

	if [ "$(OSTYPE)" == "msys" -o "$(MSYSTEM)" == "MINGW32" ]; then \
		echo "LIB = -lm" >> SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk; \
		echo "CC = gcc" >> SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk; \
		echo "BLAS = -L$(PWD)/lib -lblas" >> SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk; \
	elif [ "$(OSTYPE)" == "darwin13" ]; then \
		echo "LIB = -lm" >> SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk; \
		echo "CC = gcc -fPIC" >> SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk; \
	else \
		echo "CC = gcc -fPIC" >> SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk; \
	fi
	echo "INSTALL_LIB = $(PWD)/lib" >> SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk;
	echo "INSTALL_INCLUDE = $(PWD)/include" >> SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk;
	
# 	echo "BLAS=-L$(PWD)/lib -lblas" >> SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk;
# 	echo "LAPACK=-L$(PWD)/lib -llapack" >> SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk;

	$(MAKE) -C SuiteSparse/ library
	$(MAKE) -C SuiteSparse/ install

wget:
	if [ "$(OSTYPE)" == "msys" -o "$(MSYSTEM)" == "MINGW32" ] ; then \
		if ( wget --version ); then \
			echo "########### wget found: ok" ; \
		else \
			echo "########### Installing wget" ; \
			mingw-get.exe install msys-wget ; \
		fi; \
	fi
	
patch:
	if [ "$(OSTYPE)" == "msys" -o "$(MSYSTEM)" == "MINGW32" ] ; then \
		if ( patch --version ); then \
			echo "########### patch found: ok" ; \
		else \
			echo "########### Installing patch" ; \
			mingw-get.exe install msys-patch ; \
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

cmake:
	if [ "$(OSTYPE)" == "msys" -o "$(MSYSTEM)" == "MINGW32" ] ; then \
		if ( cmake --version ); then \
			minor=$(cmake --version | cut -d ' ' -f 3 | cut -d '.' -f 2); \
			minor=(cmake --version); \
			echo $(minor); \
			age=`cmake --version | cut -d ' ' -f 3 | cut -d'.' -f 3`; \
			echo "########### cmake found $(minor) $(age): ok" ; \
		else \
			echo "########### need cmake with version 2.8.7 or higher, please install from:" ; \
			echo "http://cmake.org/cmake/resources/software.html"; \
		fi; \
	fi
	
clean:
	echo "Cleaning ..."
	rm -rf lib include triangle AMD CAMD COLAMD CCOLAMD LDL CHOLMOD SuiteSparse lapack-* *~
	echo "done"
	
