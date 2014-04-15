#!/usr/bin/env bash

BOOST_VERSION_DEFAULT=1.53.0

checkTOOLSET(){
	ADRESSMODEL=32
	if [ "$TOOLSET" == "none" ]; then
		echo "No TOOLSET set .. using default gcc"
		SYSTEM=UNIX
		SetGCC_TOOLSET
	fi

	SRC_DIR=$PREFIX
	DIST_DIR=$PREFIX
	
	echo "-----------------------------------------------------------"
	echo "SYSTEM: $SYSTEM"
	echo "USING TOOLSET=$TOOLSET and CMAKE_GENERATOR=$CMAKE_GENERATOR"
	echo "ADRESSMODEL=$ADRESSMODEL"
	echo "-----------------------------------------------------------"
}

SetMSVC_TOOLSET(){
	TOOLSET=vc10
	SYSTEM=WIN
	export PATH=/C/Program\ Files\ \(x86\)/Microsoft\ Visual\ Studio\ 10.0/VC/bin/:$PATH
	export PATH=/C/Program\ Files\ \(x86\)/Microsoft\ Visual\ Studio\ 10.0/Common7/IDE/:$PATH
	CMAKE_GENERATOR='Visual Studio 10'
	COMPILER='msvc-10.0'
	ADRESSMODEL=32
    CPUCOUNT=1
}
SetGCC_TOOLSET(){
	TOOLSET=gcc
	
	if [ "$OSTYPE" == "msys" -o "$MSYSTEM" == "MINGW32" ]; then
		CMAKE_GENERATOR='MSYS Makefiles'
		SYSTEM=WIN
	else
		CMAKE_GENERATOR='Unix Makefiles'
	fi
	
	COMPILER='gcc'
	GCCVER=`gcc -dumpmachine`-`gcc -dumpversion`
	GCCARCH=`gcc -dumpmachine`
    CPUCOUNT=`cat /proc/cpuinfo | awk '/^processor/{print $3}' | tail -1`
	
    if [ "$CPUCOUNT"==1 ]; then
        CPUCOUNT=1
    fi

	if [ "$GCCARCH" == "mingw32" ]; then
		ADRESSMODEL=32
	else
		ADRESSMODEL=64
	fi
}
needWGET(){
    if ( wget --version ); then 
        echo "########### wget found: ok" ; 
    else 
        echo "########### Installing wget" ; 
        mingw-get.exe install msys-wget ; 
    fi; 
}
needSVN(){
    SVNEXE="svn"
}
needHG(){
    echo "Toolcheck for mecurial here"
    HG="hg"
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
needGCC(){
    echo "looking for gcc ..."
    if ( (gcc --version) );then
        echo "... found, good"
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

prepBOOST(){
    BOOST_VER=boost_${BOOST_VERSION//./_}
	BOOST_SRC=$SRC_DIR/$BOOST_VER
	BOOST_DIST=$DIST_DIR/$BOOST_VER-$TOOLSET-$ADRESSMODEL
	
	export BOOST_ROOT=$BOOST_DIST
	
	BOOST_ROOT_WIN=${BOOST_ROOT/\/c\//C:\/}
	BOOST_ROOT_WIN=${BOOST_ROOT_WIN/\/d\//D:\/}
	BOOST_ROOT_WIN=${BOOST_ROOT_WIN/\/e\//E:\/}
}
buildBOOST(){
	checkTOOLSET
	prepBOOST
    needWGET
	needPYTHON
    
    echo "--------------------------------------------------"
	echo "-------$BOOST_VER : $BOOST_SRC to $BOOST_DIST"
	echo "--------------------------------------------------"
	
	if [ ! -d $BOOST_SRC ]; then
        pushd $SRC_DIR
            wget -nc -nd http://sourceforge.net/projects/boost/files/boost/$BOOST_VERSION/$BOOST_VER'.tar.gz'
            tar -xzvf $BOOST_VER'.tar.gz'
        popd
    fi
    
	pushd $BOOST_SRC

        
		if [ "$SYSTEM" == "UNIX" ]; then
			sh bootstrap.sh
            B2="./b2"
		else
			if [ ! -f ./b2.exe ]; then
				cmd /c "bootstrap.bat "
			fi
            B2="./b2.exe"
		fi

		[ $HAVEPYTHON -eq 1 ] && WITHPYTHON='--with-python'
		
		"$B2" toolset=$COMPILER variant=release link=static,shared threading=multi address-model=$ADRESSMODEL install \
        -j $CPUCOUNT \
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
}

########## build scripts starting here
buildGCCXML(){
	checkTOOLSET
    needGIT
    needCMAKE
    
    SRCDIR=$PREFIX/gccxml
    BUILDDIR=$PREFIX/gccxml-build
    DISTDIR=$PREFIX/gccxml-bin

    echo "#########################################"
    echo "install gccxml from $SRCDIR into $DISTDIR"
    echo "#########################################"

    pushd $PREFIX
		if ( [ -d $SRCDIR ] ); then 
			pushd $SRCDIR
				"$GIT" pull
			popd
		else
			"$GIT" clone https://github.com/gccxml/gccxml.git $SRCDIR
		fi

		mkdir -p $BUILDDIR
		
		pushd $BUILDDIR
            if [ $OSTYPE = "msys" ]; then
                cmake $SRCDIR  -G 'MSYS Makefiles' \
                    -DCMAKE_INSTALL_PREFIX=$DISTDIR 
            else
                cmake $SRCDIR \
                    -DCMAKE_INSTALL_PREFIX=$DISTDIR
            fi
			
			make -j4
			make install
		popd
	popd
}
fixGCCXML(){ # check if obsolete
	GCCXML_CFG=$GCCXML_BIN_ROOT/share/gccxml-0.9/gccxml_config
	pushd $PREFIX
		rm -rf *.gch
		echo "#include <string>" > test.h
		("$GCCXML_BIN_ROOT/bin/gccxml" --debug test.h > .test.log)

		
		if [ $? -gt 0 ]; then
			echo "gccxml test fail"
			USER_FLAGS=''
			#-isystemc:mingw_w64bin
			for i in `grep "isystemc:" .test.log | sed -e 's/isystemc:mingw/isystemc:\/mingw/' | sed -e 's/bin../\/bin\/../' | tr -s '"' '\ '`; do
			#for i in `grep "isystemc:" .test.log | sed -e 's/isystemc:mingwbin../isystemc:\/mingw\/bin\/../' | tr -s '"' '\ '`; do
				echo $i
				USER_FLAGS=$USER_FLAGS' '$i
			done
			echo -e 'GCCXML_USER_FLAGS="'$USER_FLAGS'"' >> "$GCCXML_CFG"
			echo "I will now try to fix this gccxml installation ... "
			echo "You may rerun this test. If this error occur please contact me. Carsten."
		else
			echo "gccxml seems to work"
		fi
		#rm -rf test.h .test.log
    popd
}
buildPYGCCXML(){
    needHG
    needPYTHON

    SRCDIR=$PREFIX/pygccxml

    echo "#########################################"
    echo "install pygccxml at $SRCDIR"
    echo "#########################################"

    pushd $PREFIX
        if ( [ -d $SRCDIR ] ); then 
            pushd $SRCDIR
                "$HG" pull
            popd
        else
            "$HG" clone https://bitbucket.org/ompl/pygccxml $SRCDIR
        fi

        pushd $SRCDIR
            python setup.py build
            #python setup.py install
        popd
    popd
}
buildPYPLUSPLUS(){
    needHG
    needPYTHON
    
    SRCDIR=$PREFIX/pyplusplus

    echo "#########################################"
    echo "install pyplusplus at $SRCDIR"
    echo "#########################################"

    pushd $PREFIX
        if ( [ -d $SRCDIR ] ); then 
            pushd $SRCDIR
                "$HG" pull
            popd
        else
            "$HG" clone https://bitbucket.org/ompl/pyplusplus $SRCDIR
        fi

        pushd $SRCDIR
            python setup.py build
            #python setup.py install
        popd
    popd
}

showHelp(){
    echo "Compilation helper script."
    echo "--------------------------"
    echo "Install a TOOL on the current path by calling"
    echo "TOOLS are so far: boost, gccxml, pygccxml, pyplusplus"
	echo "bash $0 TOOL"
	echo ""
    echo "You can specify an installation path by setting the PREFIX variable i.e.:"
    echo "PREFIX=/path/where/to/install/ bash $0 TOOL"
    echo ""
    echo "You can specify boost version setting the BOOST_VERSION (default is $BOOST_VERSION) variable i.e.:"
    echo "BOOST_VERSION=$BOOST_VERSION bash $0 boost"
    exit
}

slotAll(){
    echo "building all"
    buildGCCXML
    buildPYGCCXML
    buildPYPLUSPLUS
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

echo "Installing at " $PREFIX

#shift 1
for arg in $@
do
    echo $arg

    case $arg in
    msvc) 
        SetMSVC_TOOLSET;;
    mingw) 
        SetGCC_TOOLSET;;
    all) 
        slotAll;;
    help)
        showHelp
        exit;;
	boost)
        buildBOOST;;
    gccxml)
        buildGCCXML;;
    pygccxml)
        buildPYGCCXML;;
    py++)
        buildPYPLUSPLUS;;
    pyplusplus)
        buildPYPLUSPLUS;;
    *) 
        echo "Don't know what to do"
        showHelp
    esac
done
