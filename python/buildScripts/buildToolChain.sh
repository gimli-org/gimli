#!/usr/bin/env bash

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

########## build scripts starting here
buildGCCXML(){
    needGIT
    needCMAKE
    needGCC

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
    echo "sh $0 TOOL"
    echo "You can specify an installation path by setting the PREFIX variable i.e.:"
    echo "PREFIX=/path/where/to/install/ sh $0 TOOL"
    echo "TOOLS ar so far: all, gccxml, pygccxml, pyplusplus"
    exit
}

slotAll(){
    echo "building all"
    buildGCCXML
    buildPYGCCXML
    buildPYPLUSPLUS
}

# script starts here 

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
        SetMINGW_TOOLSET;;
    all) 
        slotAll;;
    help)
        showHelp
        exit;;
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
