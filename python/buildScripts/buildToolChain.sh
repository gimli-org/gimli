#!/usr/bin/env bash

if [ $# -eq 0 ]; then
	prefix=`pwd`
else 
	prefix=`readlink -m $1`
fi

echo "Installing at " $prefix

GCCXML_BIN_ROOT=$prefix/gccxml-bin

echo "looking for gcc ..."
if ( (gcc --version) );then
	echo "... found, good"
else
	echo "need a working gcc installation"
	exit
fi
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
echo ""
echo "looking for svn ..."
if ( (svn --version --quiet) );then
	echo "... found, good"
else
	echo "need svn client"
	echo "get one from http://www.sliksvn.com/en/download"
	echo "if already .. ensure svn installation directory is in your PATH"
	exit
fi
echo ""
echo "Installing sources at" $prefix

echo ""
echo "looking for git ..."
if ( (git --version) );then
	GITEXE="git"
    echo "... found $GIT_EXE, good"
elif ( (/c/Program\ Files\ \(x86\)/Git/bin/git --version) );then
	GITEXE="/c/Program Files (x86)/Git/bin/git"
	echo "... found $GITEXE, good"
else
    echo "need git client"
    echo "get one from http://msysgit.github.io/"
    echo "if already .. ensure git installation directory is in your PATH"
    exit
fi
echo ""
echo "Installing sources at" $prefix

installGCCXML(){
    echo "install gccxml"
    pushd $prefix
		if ( [ -d gccxml ] ); then 
			pushd gccxml
				"$GITEXE" pull
			popd
		else
			"$GITEXE" clone --depth 1 https://github.com/gccxml/gccxml.git gccxml
		fi
		rm -rf gccxml-build $GCCXML_BIN_ROOT
		mkdir -p gccxml-build
		mkdir -p $GCCXML_BIN_ROOT
		pushd gccxml-build
            if [ $OSTYPE = "msys" ]; then
                cmake -D CMAKE_INSTALL_PREFIX=$GCCXML_BIN_ROOT ../gccxml -G 'MSYS Makefiles' 
            else
                cmake -D CMAKE_INSTALL_PREFIX=$GCCXML_BIN_ROOT ../gccxml
            fi
			
			make -j4
			make install
		popd
	popd
}

fixGCCXML(){
	GCCXML_CFG=$GCCXML_BIN_ROOT/share/gccxml-0.9/gccxml_config
	pushd $prefix
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

WORKING_PYGCC_REV=1856

installPYGCCXML(){
    echo "install pygccxml"
    pushd $prefix
        echo "getting sources ..."
        svn checkout svn://svn.code.sf.net/p/pygccxml/svn/pygccxml_dev pygccxml
        #svn co https://pygccxml.svn.sourceforge.net/svnroot/pygccxml/pygccxml_dev -r $WORKING_PYGCC_REV pygccxml
        pushd pygccxml
            python setup.py build
            #python setup.py install
        popd
    popd
}

installPYPLUSPLUS(){
    echo "install pyplusplus"
    pushd $prefix
        echo "getting sources ..."
        svn checkout svn://svn.code.sf.net/p/pygccxml/svn/pyplusplus_dev pyplusplus
        #svn co https://pygccxml.svn.sourceforge.net/svnroot/pygccxml/pyplusplus_dev -r $WORKING_PYGCC_REV pyplusplus
        pushd pyplusplus
            python setup.py build
            #python setup.py install
        popd
    popd
}

installGCCXML
installPYGCCXML
installPYPLUSPLUS
#fixGCCXML
#fixGCCXML
