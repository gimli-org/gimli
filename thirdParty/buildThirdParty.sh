#!/usr/bin/env bash

BOOST_VERSION_DEFAULT=1.57.0
BOOST_URL=http://sourceforge.net/projects/boost/files/boost/

LAPACK_VERSION=3.4.2
LAPACK_URL=http://www.netlib.org/lapack/

SUITESPARSE_VERSION=4.4.1
SUITESPARSE_URL=http://faculty.cse.tamu.edu/davis/SuiteSparse/

TRIANGLE_URL=http://www.netlib.org/voronoi/
GCCXML_URL=git://github.com/gccxml/gccxml.git

PYGCCXML_URL=https://github.com/gccxml/pygccxml
PYPLUSPLUS_URL=https://bitbucket.org/ompl/pyplusplus

CPPUNIT_URL=http://svn.code.sf.net/p/cppunit/code/trunk

GIMLI_URL=http://svn.code.sf.net/p/libgimli/code/trunk
BERT_URL=http://geo27.geo.tu-berlin.de/svn/bert

checkTOOLSET(){
	if [ "$TOOLSET" == "none" ]; then
		echo "No TOOLSET set .. using default gcc"
		SYSTEM=UNIX
		SetGCC_TOOLSET
	fi

	needPYTHON

    SRC_DIR=$PREFIX/src
    mkdir -p $SRC_DIR
    
    if [ -z "$TARGETNAME" ]; then
        TARGETNAME=-$TOOLSET-$ADDRESSMODEL
    fi

	BUILD_DIR=$PREFIX/build$TARGETNAME
	mkdir -p $BUILD_DIR

	DIST_DIR=$PREFIX/dist$TARGETNAME
	mkdir -p $DIST_DIR
		
    BUILDSCRIPT_HOME=${0%/*}
    
	echo "-----------------------------------------------------------"
	echo "SYSTEM: $SYSTEM"
	echo "USING TOOLSET=$TOOLSET and CMAKE_GENERATOR=$CMAKE_GENERATOR, PARRALELBUILD: $PARALLEL_BUILD"
	echo "ADDRESSMODEL=$ADDRESSMODEL"
	echo "PYTHON=$PYTHONVERSION"
	echo "SRCPATH=$SRC_DIR, BUILD=$BUILD_DIR, DISTRIBUTION=$DIST_DIR"
	echo "-----------------------------------------------------------"
}
SetMSVC_TOOLSET(){
	TOOLSET=vc10
	SYSTEM=WIN
	export PATH=/C/Program\ Files\ \(x86\)/Microsoft\ Visual\ Studio\ 10.0/VC/bin/:$PATH
	export PATH=/C/Program\ Files\ \(x86\)/Microsoft\ Visual\ Studio\ 10.0/Common7/IDE/:$PATH
	CMAKE_GENERATOR='Visual Studio 10'
	COMPILER='msvc-10.0'
	ADDRESSMODEL=32
    CPUCOUNT=1
}
SetGCC_TOOLSET(){
	needGCC
	TOOLSET=gcc-`gcc -dumpversion`
	B2TOOLSET=''
    MAKE=make
	
	if [ "$OSTYPE" == "msys" -o "$MSYSTEM" == "MINGW32" ]; then
		#CMAKE_GENERATOR='MSYS Makefiles'
        CMAKE_GENERATOR='MinGW Makefiles'
		SYSTEM=WIN
		B2TOOLSET=mingw
        MAKE=mingw32-make 
	elif [ "$OSTYPE" == "darwin13" ]; then
		CMAKE_GENERATOR='Unix Makefiles'
		SYSTEM=MAC
	else
		CMAKE_GENERATOR='Unix Makefiles'
		SYSTEM=UNIX
	fi
	 
    CMAKE_MAKE=$MAKE
	COMPILER='gcc'
	GCCVER=`gcc -dumpmachine`-`gcc -dumpversion`
	GCCARCH=`gcc -dumpmachine`
	CPUCOUNT=$PARALLEL_BUILD
	[ -f /proc/cpuinfo ] && CPUCOUNT=`cat /proc/cpuinfo | awk '/^processor/{print $3}' | tail -1`
    [ "$CPUCOUNT" == 0 ] && CPUCOUNT=$PARALLEL_BUILD

	if [ "$GCCARCH" == "mingw32" -o "$GCCARCH" == "i686" ]; then
		ADDRESSMODEL=32
	else
		ADDRESSMODEL=64
	fi
}
needWGET(){
    if ( which wget ); then 
        echo "########### wget found: ok" ; 
    else 
        echo "########### Installing wget" ; 
        mingw-get.exe install msys-wget ; 
    fi; 
}
needTAR(){
	if ( which tar ); then 
        echo "########### tar found: ok" ; 
    else 
        echo "########### Install tar" ; 
		  exit
        #mingw-get.exe install msys-wget ; 
    fi; 
}
needZIP(){
	if ( which unzip ); then 
        echo "########### unzip found: ok" ; 
    else 
        echo "########### Install unzip" ;
		  exit
        #mingw-get.exe install msys-wget ; 
    fi; 
}
needSED(){
	if ( which sed ); then 
        echo "########### sed found: ok" ; 
    else 
        echo "########### Install sed" ;
		  exit
        #mingw-get.exe install msys-wget ; 
    fi; 
}
getWITH_WGET(){
    needWGET
	needTAR
	needZIP
	_URL_=$1
	_SRC_=$2
	_PAC_=$3
    echo "wget $_URL_/$_PAC_"
    
	if [ ! -d $_SRC_ ]; then
        echo "Copying sources into $_SRC_"
        pushd $SRC_DIR
            wget -nc -nd $_URL_/$_PAC_
			if [ "${_PAC_##*.}" = "zip" ]; then
				unzip -o -d $_SRC_ $_PAC_ 
			else
				tar -xzvf $_PAC_
			fi
        popd
    else 
        echo "skipping .. sourcetree allready exist."
    fi

}
    
needSVN(){
    SVN="svn"
}
getWITH_SVN(){
	needSVN
    _URL_=$1
    _SRC_=$2
    _BRANCH_=$3
	
	echo "----SVN--$_URL_ -> $_SRC_ : $_BRANCH_----------------"
	echo "--------------------------------------------------"
		
	if ( [ -d $_SRC_ ] ); then 
		pushd $_SRC_
            "$SVN" up
		popd
	else
        pushd $SRC_DIR
            "$SVN" co $_URL_ $_SRC_
        popd
	fi
}
needHG(){
    echo "Toolcheck for mecurial here"
    HG="hg"
}
getWITH_HG(){
	needHG
    _URL_=$1
    _SRC_=$2
    _BRANCH_=$3
	
	echo "----HG--$_URL_ -> $_SRC_ : $_BRANCH_----------------"
	echo "--------------------------------------------------"
		
	if ( [ -d $_SRC_ ] ); then 
		pushd $_SRC_
            "$HG" pull
		popd
	else
        pushd $SRC_DIR
            "$HG" clone $_URL_ $_SRC_
        popd
	fi
}
needGIT(){
    if ( (git --version) ); then
        GIT="git"
    echo "... found $GIT_EXE, good"
    elif ( (/c/Program\ Files\ \(x86\)/Git/bin/git --version) );then
        GIT="/c/Program Files (x86)/Git/bin/git"
        echo "... found $GITEXE, good"
    else
        echo "need git client"
        echo "get one from http://msysgit.github.io/"
        echo "if already .. ensure git installation directory is in your PATH"
        exit
    fi
}
getWITH_GIT(){
	needGIT
    _URL_=$1
    _SRC_=$2
    _BRANCH_=$3
	
	echo "----GIT--$_URL_ -> $_SRC_ : $_BRANCH_----------------"
	echo "--------------------------------------------------"
		
	if ( [ -d $_SRC_ ] ); then 
		pushd $_SRC_
			#"$GIT" pull && git update $_BRANCH_
            "$GIT" pull
		popd
	else
        pushd $SRC_DIR
            "$GIT" clone $_URL_ $_SRC_
            #pushd $_SRC_
            #    "$GITEXE" pull && git update $_BRANCH_
            #popd
        popd
	fi
}
needGCC(){
    echo "looking for gcc ..."
    if ( (gcc --version) );then
        echo "... found, good"
        HAVEGCC=1
    else
        echo "need a working gcc installation"
        exit
    fi
}
needPYTHON(){
    echo ""
    echo "looking for python ..."
    if ( (python --version) );then
		HAVEPYTHON=1
        echo "... found, good"
    else
        echo "need python2.7 installation"
        echo "get one from http://www.python.org/"
        echo "if already .. ensure python27 installation directory is in your PATH"
        exit
    fi
	PYTHONVERSION=`python -c 'import sys; print(sys.version)'`
	PYTHONMAJOR=`python -c 'import sys; print(sys.version_info.major)'`
	PYTHONMINOR=`python -c 'import sys; print(sys.version_info.minor)'`
	#echo $PYTHONVERSION $PYTHONMAJOR
	
	PYTHON_HOME=`which python`
	PYTHON_HOME=${PYTHON_HOME%/*}
	
    if [ $SYSTEM == "win" ]; then
        PYTHONBASE=python$PYTHONMAJOR$PYTHONMINOR
        PYTHONDLL=$PYTHON_HOME/$PYTHONBASE.dll
        PYTHONLIB=$PYTHON_HOME/libs/lib$PYTHONBASE.a
        echo "looking for $PYTHONLIB"

        if [ ! -f $PYTHONLIB ]; then
            echo "creating"
            gendef.exe $PYTHONDLL
            dlltool.exe --dllname $PYTHONDLL --def $PYTHONBASE.def --output-lib $PYTHONLIB
            rm $PYTHONBASE.def
        else
            echo "found, ok"
        fi
    fi
	
}
needCMAKE(){
    echo ""
    echo "looking for cmake ..."
    if ( (cmake --version) );then
        echo "... found, good"
    else
        echo "need cmake"
        echo "get one from http://www.cmake.org/cmake/resources/software.html"
        echo "if already .. ensure cmake installation directory is in your PATH"
        exit
    fi
}
cmakeBuild(){
	needCMAKE
	_SRC_=$1
	_BUILD_=$2
	_DIST_=$3
	_EXTRA_=$4
	
	echo "SRC" $_SRC_
	mkBuildDIR $_BUILD_
    pushd $_BUILD_
		echo "SRC" $_SRC_
		cmake $_SRC_ -G "$CMAKE_GENERATOR" \
             -DCMAKE_MAKE_PROGRAM=$CMAKE_MAKE \
			-DCMAKE_INSTALL_PREFIX=$_DIST_ $_EXTRA_
            
        
        #make -j$PARALLEL_BUILD install
		cmake --build . --config release --target install -- -j$PARALLEL_BUILD 
    popd
}
mkBuildDIR(){
	_BUILDMB_=$1
	_SRCMB_=$2
    _FORCEMB_=$3
	if [ -n "$CLEAN" -o -n "$_FORCEMB_" ]; then
        echo " Cleaning $_BUILDMB_"
        rm -rf $_BUILDMB_
    fi
	mkdir -p $_BUILDMB_
	
	if [ -n "$_SRCMB_" ]; then # if nonzero
		echo "copy $_SRCMB_ to $_BUILDMB_"
		if [ "$(ls -A $_BUILDMB_)" ]; then
			echo " $_BUILDMB_ not empty ... do nothing"
		else
			[ -d $_SRCMB_ ] && cp -rf $_SRCMB_/* $_BUILDMB_
		fi
	fi
	
	
}
prepBOOST(){
    BOOST_VER=boost_${BOOST_VERSION//./_}
	BOOST_SRC=$SRC_DIR/$BOOST_VER
	BOOST_DIST_NAME=$BOOST_VER-$TOOLSET-$ADDRESSMODEL-py$PYTHONMAJOR$PYTHONMINOR
	BOOST_DIST=$DIST_DIR/$BOOST_DIST_NAME
	
	export BOOST_ROOT=$BOOST_DIST
	
	BOOST_ROOT_WIN=${BOOST_ROOT/\/c\//C:\/}
	BOOST_ROOT_WIN=${BOOST_ROOT_WIN/\/d\//D:\/}
	BOOST_ROOT_WIN=${BOOST_ROOT_WIN/\/e\//E:\/}
}
buildBOOST(){
	checkTOOLSET
	prepBOOST

	getWITH_WGET $BOOST_URL/$BOOST_VERSION $BOOST_SRC $BOOST_VER'.tar.gz'
	
	pushd $BOOST_SRC
        
		if [ "$SYSTEM" == "WIN" ]; then
			if [ ! -f ./b2.exe ]; then
			
				cmd /c "bootstrap.bat $B2TOOLSET" # try this first .. works for 54 with mingw
				
				if [ ! -f ./b2.exe ]; then
					./bootstrap.sh --with-toolset=$B2TOOLSET # only mingw does not work either
				fi
				#sed -e s/gcc/mingw/ project-config.jam > project-config.jam
			fi
            B2="./b2.exe"
		elif [ "$SYSTEM" == "UNIX" ]; then
			sh bootstrap.sh
            B2="./b2"
		fi
		
		[ $HAVEPYTHON -eq 1 ] && WITHPYTHON='--with-python'
		
		"$B2" toolset=$COMPILER variant=release link=static,shared threading=multi address-model=$ADDRESSMODEL install \
        -j $PARALLEL_BUILD \
		--prefix=$BOOST_DIST \
		--layout=tagged \
		$WITHPYTHON \
		--with-system \
		--with-thread \
		--with-date_time \
		--with-chrono \
		--with-regex \
		--with-filesystem \
		--with-atomic 
	popd
	echo $BOOST_DIST_NAME > $DIST_DIR/.boost.dist
}

prepGCCXML(){
	GCCXML_VER=gccxml
	GCCXML_SRC=$SRC_DIR/$GCCXML_VER
	GCCXML_BUILD=$BUILD_DIR/$GCCXML_VER
	GCCXML_DIST=$DIST_DIR
}
buildGCCXML(){
	checkTOOLSET
	prepGCCXML
    needCMAKE
    
 	getWITH_GIT $GCCXML_URL $GCCXML_SRC
	
	EXTRA=''
	if [ "$OSTYPE" == "darwin13" ]; then
			# on APPLE will wrong insert -no-cpp-precomp, avoid it with -DCMAKE_C_COMPILER_ID=GNUAPPLE \
		EXTRA='-DCMAKE_C_COMPILER_ID=GNUAPPLE'
	fi

    if [ "$SYSTEM" == "WIN" ]; then
        # on windows system gccxml with 64bit seems to be broken so we build with default 32 bit compiler
        cp /etc/fstab /etc/fstab_back
        OLDPATH=$PATH
        umount /mingw
        mount c:\\MINGW/ /mingw
        export PATH=/c/MINGW/bin:$PATH
        cmakeBuild $GCCXML_SRC $GCCXML_BUILD $GCCXML_DIST $EXTRA
        cp /etc/fstab_back /etc/fstab
        export PATH=$OLDPATH
    else
        cmakeBuild $GCCXML_SRC $GCCXML_BUILD $GCCXML_DIST $EXTRA
    fi
	 
}
prepPYGCCXML(){
	PYGCCXML_VER=pygccxml
	PYGCCXML_SRC=$SRC_DIR/$PYGCCXML_VER
	PYGCCXML_BUILD=$BUILD_DIR/$PYGCCXML_VER
	PYGCCXML_DIST=$DIST_DIR/$PYGCCXML_VER
    
    PYPLUSPLUS_VER=pyplusplus
    PYPLUSPLUS_SRC=$SRC_DIR/$PYPLUSPLUS_VER
	PYPLUSPLUS_BUILD=$BUILD_DIR/$PYPLUSPLUS_VER
	PYPLUSPLUS_DIST=$DIST_DIR/$PYPLUSPLUS_VER
        
    PYGCCXML_DIST_WIN=${PYGCCXML_DIST/\/c\//C:\\/}
	PYGCCXML_DIST_WIN=${PYGCCXML_DIST_WIN/\/d\//D:\\/}
    PYGCCXML_DIST_WIN=${PYGCCXML_DIST_WIN/\/e\//E:\\/}
}

buildPYGCCXML(){
	checkTOOLSET
	prepPYGCCXML

	getWITH_GIT $PYGCCXML_URL $PYGCCXML_SRC
	getWITH_HG $PYPLUSPLUS_URL $PYPLUSPLUS_SRC
	
	mkBuildDIR $PYGCCXML_BUILD $PYGCCXML_SRC
    pushd $PYGCCXML_BUILD
		python setup.py build
		echo "copy build->dist"
		cp -rf $PYGCCXML_BUILD/build/lib*/pygccxml $PYGCCXML_DIST
		#export PYTHONPATH=$PYTHONPATH:$PYGCCXML_DIST/Lib/site_packages/
		#python setup.py install --prefix=$PYGCCXML_DIST_WIN
    popd
	
	mkBuildDIR $PYPLUSPLUS_BUILD $PYPLUSPLUS_SRC
	pushd $PYPLUSPLUS_BUILD
		python setup.py build
		echo "copy build->dist"
		cp -rf $PYPLUSPLUS_BUILD/build/lib*/pyplusplus $PYPLUSPLUS_DIST
		#export PYTHONPATH=$PYTHONPATH:$PYGCCXML_DIST/Lib/site_packages/
		#python setup.py install --prefix=$PYGCCXML_DIST_WIN
    popd
}

prepLAPACK(){
	LAPACK_VER=lapack-$LAPACK_VERSION
	LAPACK_SRC=$SRC_DIR/$LAPACK_VER
	LAPACK_BUILD=$BUILD_DIR/$LAPACK_VER
	LAPACK_DIST=$DIST_DIR
}
buildLAPACK(){
	echo "############### LAPACK ##########################"
	checkTOOLSET
	prepLAPACK
	getWITH_WGET $LAPACK_URL $LAPACK_SRC $LAPACK_VER.tgz
	EXTRA="-DBUILD_SHARED_LIBS=ON"
	cmakeBuild $LAPACK_SRC $LAPACK_BUILD $LAPACK_DIST $EXTRA
}
prepTRIANGLE(){
	TRIANGLE_VER=triangle
	TRIANGLE_SRC=$SRC_DIR/$TRIANGLE_VER
	TRIANGLE_BUILD=$BUILD_DIR/$TRIANGLE_VER
	TRIANGLE_DIST=$DIST_DIR/lib
}
buildTRIANGLE(){
	echo "############### TRIANGLE ##########################"
	checkTOOLSET
	prepTRIANGLE
	getWITH_WGET $TRIANGLE_URL $TRIANGLE_SRC $TRIANGLE_VER.zip
	needSED
	rm -rf $TRIANGLE_BUILD
	mkBuildDIR $TRIANGLE_BUILD $TRIANGLE_SRC
	
	pushd $TRIANGLE_BUILD
		if [ "$SYSTEM" == "WIN" ]; then 
			sed -i -e 's/-DLINUX/-DCPU86/g' makefile ;
			patch triangle.c -i $BUILDSCRIPT_HOME/patches/triangle-mingw-win64.patch
		fi
	
		if [ "$ADDRESSMODEL" == "64" ]; then
			sed -i -e 's/CC = cc/CC = gcc -fPIC/g' makefile;
		else
			sed -i -e 's/CC = cc/CC = gcc/g' makefile;
		fi
	
		make trilibrary
		mkdir -p $TRIANGLE_DIST
		ar cqs $TRIANGLE_DIST/libtriangle.a triangle.o
        mkdir -p $DIST_DIR/include
        cp triangle.h $DIST_DIR/include
	popd
}
prepSUITESPARSE(){
	SUITESPARSE_VER=SuiteSparse-$SUITESPARSE_VERSION
	SUITESPARSE_SRC=$SRC_DIR/$SUITESPARSE_VER
	SUITESPARSE_BUILD=$BUILD_DIR/$SUITESPARSE_VER
	SUITESPARSE_DIST=$DIST_DIR
}
buildSUITESPARSE(){
	checkTOOLSET
	prepLAPACK
	prepSUITESPARSE
	
	getWITH_WGET $SUITESPARSE_URL $SUITESPARSE_SRC $SUITESPARSE_VER.tar.gz
	[ -d $SRC_DIR/SuiteSparse ] && mv $SRC_DIR/SuiteSparse $SUITESPARSE_SRC
	
	mkBuildDIR $SUITESPARSE_BUILD $SUITESPARSE_SRC 1
	
	pushd $SUITESPARSE_BUILD
		if [ "$SYSTEM" == "WIN" ]; then 
			patch -p1 < $SRC_DIR/patches/SuiteSparse-4.4.1.patch
			echo "LIB = -lm" >> SuiteSparse_config/SuiteSparse_config.mk
			echo "CC = gcc" >> SuiteSparse_config/SuiteSparse_config.mk
			echo "BLAS = -L$TRIANGLE_DIST -lblas" >> SuiteSparse_config/SuiteSparse_config.mk
		elif [ "$OSTYPE" == "darwin13" ]; then
			echo "LIB = -lm" >> SuiteSparse_config/SuiteSparse_config.mk
			echo "CC = gcc -fPIC" >> SuiteSparse_config/SuiteSparse_config.mk
		else 
			echo "CC = gcc " >> SuiteSparse_config/SuiteSparse_config.mk
		fi
	
		mkdir -p $SUITESPARSE_DIST/lib
		mkdir -p $SUITESPARSE_DIST/include
		echo "INSTALL_LIB = $SUITESPARSE_DIST/lib" >> SuiteSparse_config/SuiteSparse_config.mk;
		echo "INSTALL_INCLUDE = $SUITESPARSE_DIST/include" >> SuiteSparse_config/SuiteSparse_config.mk;

		MODULE='.'
		if [ -n "$1" ]; then
			MODULES=$1
		fi
		echo "Installing $MODULES"
		pushd $MODULE
            		"$MAKE" -j$PARALLEL_BUILD library
		        "$MAKE" install
		popd 
	popd
}

prepCPPUNIT(){
    CPPUNIT_VER=cppunit
	CPPUNIT_SRC=$SRC_DIR/$CPPUNIT_VER
	CPPUNIT_BUILD=$BUILD_DIR/$CPPUNIT_VER
	CPPUNIT_DIST=$DIST_DIR
}
buildCPPUNIT(){
    checkTOOLSET
    prepCPPUNIT
    getWITH_SVN $CPPUNIT_URL $CPPUNIT_SRC 
	mkBuildDIR $CPPUNIT_BUILD $CPPUNIT_SRC
    
    pushd $CPPUNIT_SRC
        pushd cppunit
            ./autogen.sh 
            ./configure --prefix=$CPPUNIT_DIST
            make -j $PARALLEL_BUILD
            make install
            
        popd
    popd
    
}
# prepGIMLI(){
	# GIMLI_VER=gimli
	# GIMLI_SRC=$SRC_DIR/$GIMLI_VER
	# GIMLI_BUILD=$BUILD_DIR/$GIMLI_VER
	# GIMLI_DIST=$DIST_DIR/$GIMLI_VER-py$PYTHONMAJOR$PYTHONMINOR
# }

# buildGIMLI(){
	# checkTOOLSET
	# prepGIMLI
	# prepBOOST
	# prepLAPACK
    # prepGCCXML
	# prepSUITESPARSE

	# getWITH_SVN $GIMLI_URL $GIMLI_SRC 
	# mkBuildDIR $GIMLI_BUILD

	# pushd $GIMLI_BUILD
		# echo "cmake $GIMLI_SRC -G "$CMAKE_GENERATOR" \
			# -DCMAKE_INSTALL_PREFIX=$GIMLI_DIST \
			# -DBOOST_ROOT=$BOOST_ROOT \
			# -DEXTERNAL_DIR=$DIST_DIR \
			# -DGCCXML_EXECUTABLE=$GCCXML_DIST/bin/gccxml.exe"
		# cmake $GIMLI_SRC -G "$CMAKE_GENERATOR" \
            # -DCMAKE_MAKE_PROGRAM=$CMAKE_MAKE \
			# -DCMAKE_INSTALL_PREFIX=$GIMLI_DIST \
			# -DBOOST_ROOT=$BOOST_ROOT \
			# -DEXTERNAL_DIR=$DIST_DIR \
			# -DGCCXML_EXECUTABLE=$GCCXML_DIST/bin/gccxml.exe	
		# "$MAKE" -j$PARALLEL_BUILD 
		
        # OLDPATH=$PATH
		# export PATH=/c/mingw/bin:$PATH
		# export PYTHONPATH=$PYTHONPATH:$DIST_DIR
		# "$MAKE" pygimli J=$PARALLEL_BUILD
        # export PATH=$OLDPATH
		
		# #	-DLAPACK_LIBRARIES=C:/Users/CarstenRuecker/src/gimli/trunk/external/lib/liblapack.dll
	# popd
# }
# prepBERT(){
	# BERT_VER=bert
	# BERT_SRC=$SRC_DIR/$BERT_VER
	# BERT_BUILD=$BUILD_DIR/$BERT_VER
	# BERT_DIST=$DIST_DIR/$BERT_VER-py$PYTHONMAJOR$PYTHONMINOR
# }
# buildBERT(){
	# checkTOOLSET
	# prepBERT
	# prepGIMLI
	# prepBOOST
    # prepGCCXML
	
	# getWITH_SVN $BERT_URL $BERT_SRC 
	# mkBuildDIR $BERT_BUILD

	# pushd $BERT_BUILD
		# cmake $BERT_SRC/trunk -G "$CMAKE_GENERATOR" \
			# -DCMAKE_INSTALL_PREFIX=$BERT_DIST \
            # -DCMAKE_MAKE_PROGRAM=$CMAKE_MAKE \
			# -DGIMLI_SRC=$GIMLI_SRC \
			# -DGIMLI_BUILD=$GIMLI_BUILD \
			# -DBOOST_ROOT=$BOOST_ROOT \
			# -DEXTERNAL_DIR=$DIST_DIR \
			# -DGCCXML_EXECUTABLE=$GCCXML_DIST/bin/gccxml.exe
		
		# "$MAKE" -j$PARALLEL_BUILD
		# export PATH=$PATH:/c/mingw/bin
		# export PYTHONPATH=$PYTHONPATH:$DIST_DIR
		# "$MAKE" pybert J=$PARALLEL_BUILD

	# popd
# }

slotAll(){
	buildBOOST
	buildLAPACK
	buildTRIANGLE
	buildSUITESPARSE
	buildGCCXML
	buildPYGCCXML
}

showHelp(){
	echo "boost | lapack | triangle | suitesparse | gccxml | pygccxml | all"
}

# script starts here 
TOOLSET=none

if [ -n "$BOOST_VERSION" ]; then
	BOOST_VERSION=$BOOST_VERSION
else
	BOOST_VERSION=$BOOST_VERSION_DEFAULT
fi

if [ -n "$PREFIX" ]; then
    PREFIX=`readlink -m $PREFIX`
else 
    PREFIX=`pwd`
fi

if [ -z "$PARALLEL_BUILD" ]; then
    PARALLEL_BUILD=1
fi
echo "Installing at " $PREFIX


CMAKE_BUILD_TYPE=Release

for arg in $@
do
	echo $arg
    case $arg in
    msvc) 
        SetMSVC_TOOLSET;;
	mingw) 
        SetMINGW_TOOLSET;;
    all) 
        slotAll;;
    help)
        showHelp
		exit;;
	boost)
		buildBOOST;;
	lapack)
		buildLAPACK;;
	triangle)
		buildTRIANGLE;;
	suitesparse)
		buildSUITESPARSE;;
    umfpack)
        buildSUITESPARSE UMFPACK;;
	gccxml)
		buildGCCXML;;
	pygccxml)
		buildPYGCCXML;;
    cppunit)
		buildCPPUNIT;;
	
    *) 
        echo "Don't know what to do."
        showHelp;;
    esac
done
