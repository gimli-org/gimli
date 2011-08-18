#!/usr/bin/env bash

prefix=`pwd`

if ( (gcc --version) );then
	echo "looking for gcc ... found, good"
	echo "... found, good"
else
	echo "need a working gcc installation"
	exit
fi

if ( (python --version) );then
	echo "looking for python ... found, good"
else
	echo "need python2.6 installation"
	echo "get one from http://www.python.org/"
	echo "if allready ensure python26 installation directory is in your PATH"
	exit
fi

if ( (cmake --version) );then
	echo "looking for cmake ... found, good"
else
	echo "need cmake"
	echo "get one from http://www.cmake.org/cmake/resources/software.html"
	echo "if allready ensure cmake installation directory is in your PATH"
	exit
fi

if ( (svn --version --quiet) );then
	echo "lookcing for svn ... found, good"
else
	echo "need svn client"
	echo "get one from http://www.sliksvn.com/en/download"
	echo "if allready ensure svn installation directory is in your PATH"
	exit
fi

echo "Installing sources at" $prefix

installGCCXML(){
    echo "install gccxml"
    oldpwd=`pwd`
    cd $prefix
    cvs -d :pserver:anoncvs@www.gccxml.org:/cvsroot/GCC_XML co gccxml/
    cd gccxml
    cmake ./ -G 'MSYS Makefiles'
    make
    make install

    cd $oldpwd
}

fixGCCXML(){
    oldpwd=`pwd`
    cd $prefix
    echo "#include <string>" > test.h
    if [ !$(/c/programme/gccxml/bin/gccxml --debug test.h >.test.log) ]; then
        echo "gccxml test fail"
        COMPILER=`grep "GCCXML_COMPILER" /c/programme/gccxml/share/gccxml-0.9/gccxml_config | cut -f2 -d'=' | cut -f2 -d':'`
        
        echo $COMPILER
        echo ${COMPILER%'"'}
        grep "isystemc:" .test.log
    else
        echo "gccxml seems to work"
    fi
    
    cd $oldpwd
}

installPYGCCXML(){
    echo "install pygccxml"
    oldpwd=`pwd`
    cd $prefix
	
	svn co https://pygccxml.svn.sourceforge.net/svnroot/pygccxml/pygccxml_dev -r 1842 pygccxml
    cd pygccxml
    python setup.py install
    cd $oldpwd   
}

installPYPLUSPLUS(){
    echo "install pyplusplus"
    oldpwd=`pwd`
    cd $prefix
	svn co https://pygccxml.svn.sourceforge.net/svnroot/pygccxml/pyplusplus_dev -r 1842 pyplusplus
    cd pyplusplus
    python setup.py install
    cd $oldpwd   
}

#installGCCXML
#installPYGCCXML
#installPYPLUSPLUS

fixGCCXML
