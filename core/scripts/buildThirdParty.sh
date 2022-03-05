#!/usr/bin/env bash

BOOST_VERSION_DEFAULT=1.68.0
#since 63 libboost_numpy
#since 64 python build broken

BOOST_URL=http://sourceforge.net/projects/boost/files/boost/

LAPACK_VERSION=3.4.2
LAPACK_URL=http://www.netlib.org/lapack/

#SUITESPARSE_VERSION=5.2.0
#SUITESPARSE_URL=http://faculty.cse.tamu.edu/davis/SuiteSparse/

SUITESPARSE_URL=https://github.com/DrTimothyAldenDavis/SuiteSparse
SUITESPARSE_REV=v5.8.1

TRIANGLE_URL=http://www.netlib.org/voronoi/

CASTXML_URL=https://github.com/CastXML/CastXML.git
#CASTXML_REV=d5934bd08651dbda95a65ccadcc5f39637d7bc59 #current functional
#CASTXML_REV=9d7a46d639ce921b8ddd36ecaa23c567d003294a #last functional

# Check for updates https://data.kitware.com/#search/results?query=castxml&mode=text
CASTXML_BIN_LINUX=https://data.kitware.com/api/v1/file/5b6c5b4d8d777f06857c323b/download
#CASTXML_BIN_LINUX=https://data.kitware.com/api/v1/item/57b5de948d777f10f2696371/download # seems broken
#CASTXML_BIN_LINUX=https://data.kitware.com/api/v1/file/57b5dea08d777f10f2696379/download # old ok
CASTXML_BIN_MAC=https://data.kitware.com/api/v1/file/57b5de9f8d777f10f2696378/download
CASTXML_BIN_WIN=https://data.kitware.com/api/v1/file/5b68bfc28d777f06857c1f44/download

PYGCCXML_URL=https://github.com/gccxml/pygccxml
PYGCCXML_REV=84be3367bf43cb494512f343068cb23704a47460 # for py3.8

# old bitbucked project not working anymore and moved to github
# PYPLUSPLUS_URL=https://bitbucket.org/ompl/pyplusplus
# PYPLUSPLUS_REV=1e30641 # tag 1.8.3 for py3.8
PYPLUSPLUS_URL=https://github.com/ompl/pyplusplus
PYPLUSPLUS_REV=d4811c8 # tag 1.8.3 for py3.8

CPPUNIT_URL=http://svn.code.sf.net/p/cppunit/code/trunk

checkTOOLSET(){
    if [ "$TOOLSET" == "none" ]; then
        echo "No TOOLSET set .. using default gcc"
        SYSTEM=UNIX
        SetGCC_TOOLSET
    elif [ "$TOOLSET" == "clang" ]; then
        SetCLANG_TOOLSET
    fi

    needPYTHON

    SRC_DIR=$GIMLI_PREFIX/src
    mkdir -p $SRC_DIR

    if [ -z "$TARGETNAME" ]; then
        TARGETNAME=-$TOOLSET-$ADDRESSMODEL
    fi

    BUILD_DIR=$GIMLI_PREFIX/build$TARGETNAME
    mkdir -p $BUILD_DIR

    DIST_DIR=$GIMLI_PREFIX/dist$TARGETNAME
    mkdir -p $DIST_DIR

    BUILDSCRIPT_HOME=${0%/*}

    echo "##################################################################################################"
    echo "SYSTEM: $SYSTEM"
    echo "USING TOOLSET=$TOOLSET and CMAKE_GENERATOR=$CMAKE_GENERATOR, PARRALELBUILD: $PARALLEL_BUILD"
    echo "ADDRESSMODEL = $ADDRESSMODEL"
    echo "PYTHON       = $PYTHONEXE"
    echo "Version      = $PYTHONVERSION"
    echo "SRC_DIR      = $SRC_DIR"
    echo "BUILD_DIR    = $BUILD_DIR"
    echo "DIST_DIR     = $DIST_DIR"
    echo "##################################################################################################"
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
        CMAKE_GENERATOR='Unix Makefiles'
        #CMAKE_GENERATOR='MinGW Makefiles'
        SYSTEM=WIN
        B2TOOLSET=mingw
        B2TOOLSET=gcc
        #MAKE=mingw32-make
    elif [ $(uname) == "Darwin" ]; then
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

    if [ "$GCCARCH" == "mingw32" -o "$GCCARCH" == "i686" -o "$GCCARCH" == "i686-pc-msys" -o "$GCCARCH" == "i686-w64-mingw32" ]; then
        ADDRESSMODEL=32
    else
        ADDRESSMODEL=64
    fi
}

getCLANG_NAME() {
    # this is probably not longer used
    if command -v clang-3.6 2>/dev/null; then
        CLANG="clang-3.6"
        CLANGPP="clang++-3.6"
    elif command -v clang-3.6.1 2>/dev/null; then
        CLANG="clang-3.6.1"
        CLANGPP="clang++-3.6.1"
    else
        CLANG="clang"
        CLANGPP="clang++"
    fi

    echo "Setting clang to: " $CLANG
}

SetCLANG_TOOLSET(){
    #needGCC
    TOOLSET=clang-`clang -dumpversion`
    B2TOOLSET='clang'
    MAKE=make

    if [ "$OSTYPE" == "msys" -o "$MSYSTEM" == "MINGW32" ]; then
        CMAKE_GENERATOR='Unix Makefiles'
        SYSTEM=WIN
    elif [ $(uname) == "Darwin" ]; then
        CMAKE_GENERATOR='Unix Makefiles'
        SYSTEM=MAC
    else
        CMAKE_GENERATOR='Unix Makefiles'
        SYSTEM=UNIX
    fi

    CMAKE_MAKE=$MAKE
    COMPILER='clang'
    GCCVER=`clang -dumpmachine`-`gcc -dumpversion`
    GCCARCH=`clang -dumpmachine`
    CPUCOUNT=$PARALLEL_BUILD
    [ -f /proc/cpuinfo ] && CPUCOUNT=`cat /proc/cpuinfo | awk '/^processor/{print $3}' | tail -1`
    [ "$CPUCOUNT" == 0 ] && CPUCOUNT=$PARALLEL_BUILD

    if [ "$GCCARCH" == "mingw32" -o "$GCCARCH" == "i686" -o "$GCCARCH" == "i686-pc-msys" -o "$GCCARCH" == "i686-w64-mingw32" ]; then
        ADDRESSMODEL=32
    else
        ADDRESSMODEL=64
    fi
}

getWITH_WGET(){
    _URL_=$1
    _SRC_=$2
    _PAC_=$3

    echo "** get with wget ** " $_URL_ $_SRC_ $_PAC_
    echo "wget -nc -nd $_URL_/$_PAC_"

    if [ -n "$CLEAN" ]; then
        rm -rf $_SRC_
        #rm -rf $_PAC_
    fi

    if [ ! -d $_SRC_ ]; then
        echo "Copying sources into $_SRC_"
        pushd $SRC_DIR
            if [[ `wget -S --spider $_URL_/$_PAC_  2>&1 | grep 'HTTP/1.1 200 OK'` ]]; then
                echo $_URL_/$_PAC_ " found. Downloading."
                wget -nc -nd $_URL_/$_PAC_
            else
                echo $_URL_/$_PAC_ " not found. Downloading index.html and copying to " $_PAC_
                wget -nc -nd $_URL_

                echo "cp download" $_PAC_
                cp download download.tar.gz
                cmake -E tar -xzf download.tar.gz
                mv $_PAC_ download.dir/
                echo "cp download.dir/$_PAC_ ."
                mv download.dir/$_PAC_ .
            fi
            
            if [ "${_PAC_##*.}" = "zip" ]; then
                mkdir -p $_SRC_
                pushd $_SRC_
                    cmake -E tar -xz ../$_PAC_
                popd
            else
                cmake -E tar -xzf $_PAC_
            fi
        popd
    else
        echo "Skipping .. source tree already exists. Use with CLEAN=1 if you want to force installation."
    fi

}

getWITH_SVN(){
    SVN="svn"
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

getWITH_HG(){
    HG="hg"
    _URL_=$1
    _SRC_=$2
    _BRANCH_=$3

    echo "----HG--$_URL_ -> $_SRC_ : $_BRANCH_----------------"
    echo "--------------------------------------------------"

    if ( [ -d $_SRC_ ] ); then
        pushd $_SRC_
#            "$HG" fetch
            $HG pull -u
            $HG up
        popd
    else
        pushd $SRC_DIR
            ${HG} --config ui.clonebundles=false clone ${_URL_} ${_SRC_}
        popd
    fi
    if [ -n $_BRANCH_ ]; then
        pushd $_SRC_
          echo $_SRC_ $_BRANCH_
          $HG checkout $_BRANCH_
          #"$HG" revert -r $_BRANCH_ --all
        popd
    fi
}

getWITH_GIT(){
    GIT="git"
    _URL_=$1
    _SRC_=$2
    _BRANCH_=$3

    echo "----GIT--$_URL_ -> $_SRC_ : $_BRANCH_----------------"
    echo "--------------------------------------------------"

    if ( [ -d $_SRC_ ] ); then
        pushd $_SRC_
            "$GIT" stash
            "$GIT" pull
         popd
    else
        pushd $SRC_DIR
            "$GIT" clone $_URL_ $_SRC_
        popd
    fi
    if [ -n $_BRANCH_ ]; then
        pushd $_SRC_
          echo $_SRC_ $_BRANCH_
          echo $GIT checkout --force $_BRANCH_ .
          "$GIT" checkout $_BRANCH_ .
        popd
    fi

}
needGCC(){
    HAVEGCC=1
}
needPYTHON(){

    echo "*** Need python "

    if [ "$SYSTEM" == "WIN" ]; then
        # don't check for python3 in win .. it will find msys python not anaconda
        if command -v python 2>/dev/null; then
            PYTHONEXE=python
        else
            echo "Cannot find python interpreter. Please consider install anaconda."
        fi
    else
        if command -v python3 2>/dev/null; then
            PYTHONEXE=python3
        elif command -v python 2>/dev/null; then
            PYTHONEXE=python
        else
            echo "cannot find python interpreter"
        fi
    fi

    HAVEPYTHON=1
    PYTHONVERSION=`"$PYTHONEXE" -c 'import sys; print(sys.version)'`
    PYTHONMAJOR=`"$PYTHONEXE" -c 'import sys; print(sys.version_info.major)'`
    PYTHONMINOR=`"$PYTHONEXE" -c 'import sys; print(sys.version_info.minor)'`
    echo "Python version found py" $PYTHONMAJOR $PYTHONMINOR $PYTHONVERSION 

    PYTHON_HOME=`which $PYTHONEXE`
    PYTHON_HOME=${PYTHON_HOME%/*}

    echo "PYTHON HOME: $PYTHON_HOME"
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
            echo "found, ok: $PYTHONLIB "
        fi
    fi
}

cmakeBuild(){
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
    BOOST_DIST_NAME=$BOOST_VER-$TOOLSET-$ADDRESSMODEL-'py'$PYTHONMAJOR$PYTHONMINOR
    BOOST_DIST=$DIST_DIR/$BOOST_DIST_NAME
    BOOST_BUILD=$BUILD_DIR/$BOOST_VER-'py'$PYTHONMAJOR$PYTHONMINOR
    export BOOST_ROOT=$BOOST_DIST

    BOOST_ROOT_WIN=${BOOST_ROOT/\/c\//C:\/}
    BOOST_ROOT_WIN=${BOOST_ROOT_WIN/\/d\//D:\/}
    BOOST_ROOT_WIN=${BOOST_ROOT_WIN/\/e\//E:\/}
}
buildBOOST(){
    echo "*** building boost ... "
    
    checkTOOLSET
    prepBOOST

    getWITH_WGET $BOOST_URL/$BOOST_VERSION $BOOST_SRC $BOOST_VER'.tar.gz'

    if [ ! -d $BOOST_BUILD ]; then
        echo "copying sourcetree into build: $BOOST_BUILD"
        cp -r $BOOST_SRC $BOOST_BUILD
    fi
    pushd $BOOST_BUILD
        echo "Try to build b2 for TOOLSET: $B2TOOLSET"

        if [ "$SYSTEM" == "WIN" ]; then
            if [ ! -f ./b2.exe ]; then
                echo "Try with cmd /c \"bootstrap.bat $B2TOOLSET\""
                cmd.exe /c "bootstrap.bat $B2TOOLSET" # try this first .. works for 54 with mingw

                if [ ! -f ./b2.exe ]; then
                    echo "Try with ./bootstrap.sh --with-toolset=$B2TOOLSET"
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

        #"$B2" toolset=$COMPILER --verbose-test test

        #quit
        if [ $SYSTEM == 'WIN' -a $ADDRESSMODEL == '64' ]; then
            EXTRADEFINES='define=BOOST_USE_WINDOWS_H define=MS_WIN64'
            echo "+++++++++++++++++++ :$EXTRADEFINES"
        fi
        echo "Build with python: $WITHPYTHON"

        "$B2" toolset=$COMPILER variant=release link=static,shared threading=multi address-model=$ADDRESSMODEL $EXTRADEFINES install \
        -j $PARALLEL_BUILD \
        -d 0 \
        -a \
        --prefix=$BOOST_DIST \
        --platform=msys \
        --layout=tagged \
        --debug-configuration \
        $WITHPYTHON 
    
    	#--with-system \
        #--with-thread 

        # --with-date_time \
        # --with-chrono \
        # --with-regex \
        # --with-filesystem \
        # --with-atomic
    popd
    echo $BOOST_DIST_NAME > $DIST_DIR/.boost-py$PYTHONMAJOR.dist
}

prepCASTXMLBIN(){
    CASTXML_VER=castxml
    CASTXML_SRC=$SRC_DIR/$CASTXML_VER
    CASTXML_BUILD=$BUILD_DIR/$CASTXML_VER
    CASTXML_DIST=$DIST_DIR
}

buildCASTXMLBIN(){
    checkTOOLSET
    prepCASTXMLBIN

    if [ "$SYSTEM" == "WIN" ]; then
        if [ -n "$CLEAN" ]; then
            rm -f $SRC_DIR/castxml-windows.zip
        fi
        getWITH_WGET $CASTXML_BIN_WIN $CASTXML_SRC castxml-windows.zip
        cp -r $CASTXML_SRC/castxml/* $CASTXML_DIST
        CASTXMLBIN=castxml.exe
    elif [ "$SYSTEM" == "MAC" ]; then
        getWITH_WGET $CASTXML_BIN_MAC $CASTXML_SRC castxml-macosx.tar.gz
        cp -r $CASTXML_SRC/* $CASTXML_DIST
        CASTXMLBIN=castxml
    else
        if [ -n "$CLEAN" ]; then
            rm -f $SRC_DIR/castxml-linux.tar.gz
        fi
        getWITH_WGET $CASTXML_BIN_LINUX $CASTXML_SRC castxml-linux.tar.gz
        cp -r $CASTXML_SRC/* $CASTXML_DIST
        CASTXMLBIN=castxml
    fi

    if "$CASTXML_DIST/bin/$CASTXMLBIN" --version; then
        echo "Binary castxml seems to work"
    else
        echo "Binary castxml does not seems to work. Removing binary installation"
        rm $CASTXML_DIST/bin/*
        rm -rf $CASTXML_DIST/share/castxml
    fi
}

prepCASTXML(){
    CASTXML_VER=castxmlSRC
    CASTXML_SRC=$SRC_DIR/$CASTXML_VER
    CASTXML_BUILD=$BUILD_DIR/$CASTXML_VER
    CASTXML_DIST=$DIST_DIR
}

buildCASTXML(){
    checkTOOLSET
    prepCASTXML

    getWITH_GIT $CASTXML_URL $CASTXML_SRC $CASTXML_REV

    if [ "$SYSTEM" == "WIN" ]; then

        mkBuildDIR $CASTXML_BUILD
        pushd $CASTXML_BUILD
            if [ $ADDRESSMODEL == '32' ]; then
                CC=/mingw32/bin/gcc CXX=/mingw32/bin/g++ cmake $CASTXML_SRC -G "$CMAKE_GENERATOR" \
                    -DCMAKE_MAKE_PROGRAM=$CMAKE_MAKE \
                    -DCMAKE_INSTALL_PREFIX=$CASTXML_DIST

                sed -i -e 's/\/mingw32/C:\/msys32\/mingw32/g' src/CMakeFiles/castxml.dir/linklibs.rsp
            else
                                #CC=/mingw64/bin/gcc CXX=/mingw64/bin/g++ cmake $CASTXML_SRC -G "$CMAKE_GENERATOR" \
                #    -DCMAKE_MAKE_PROGRAM=$CMAKE_MAKE \
                #    -DCMAKE_INSTALL_PREFIX=$CASTXML_DIST
                #sed -i -e 's/\/mingw64/C:\/msys64\/mingw64/g' src/CMakeFiles/castxml.dir/linklibs.rsp

                sed -i 's/if(DEFINED LLVM_BUILD_BINARY_DIR)/if(DEFINED LLVM_BUILD_BINARY_DIR_)/' $CASTXML_SRC/CMakeLists.txt

                #-DCLANG_RESOURCE_DIR=/mingw64/lib/clang/3.7.0/

                cmake $CASTXML_SRC -G "$CMAKE_GENERATOR" \
                    -DCMAKE_MAKE_PROGRAM=$CMAKE_MAKE \
                    -DCMAKE_INSTALL_PREFIX=$CASTXML_DIST
            fi

            #make -j$PARALLEL_BUILD install
            cmake --build . --config release --target install -- -j$PARALLEL_BUILD
        popd
    elif [ "$SYSTEM" == "MAC" ]; then
        CC=/usr/bin/gcc CXX=/usr/bin/g++ cmakeBuild $CASTXML_SRC $CASTXML_BUILD $CASTXML_DIST \
          "-DLLVM_DIR=/usr/local/Cellar/llvm/HEAD/share/llvm/cmake/"
    else
        # getCLANG_NAME not longer used
        #CC=$CLANG CXX=$CLANGPP cmakeBuild $CASTXML_SRC $CASTXML_BUILD $CASTXML_DIST
        cmakeBuild $CASTXML_SRC $CASTXML_BUILD $CASTXML_DIST

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

    getWITH_GIT $PYGCCXML_URL $PYGCCXML_SRC $PYGCCXML_REV
    getWITH_GIT $PYPLUSPLUS_URL $PYPLUSPLUS_SRC $PYPLUSPLUS_REV

    mkBuildDIR $PYGCCXML_BUILD $PYGCCXML_SRC 1
    pushd $PYGCCXML_BUILD
        "$PYTHONEXE" setup.py build
        echo "copy build->dist"
        if [ -n "$CLEAN" ]; then
            rm -rf $PYGCCXML_DIST
        fi

        cp -rf $PYGCCXML_BUILD/build/lib*/pygccxml $PYGCCXML_DIST
        #export PYTHONPATH=$PYTHONPATH:$PYGCCXML_DIST/Lib/site_packages/
        #python setup.py install --prefix=$PYGCCXML_DIST_WIN
    popd

    mkBuildDIR $PYPLUSPLUS_BUILD $PYPLUSPLUS_SRC 1
    pushd $PYPLUSPLUS_BUILD
        "$PYTHONEXE" setup.py build
        echo "copy build->dist"
        if [ -n "$CLEAN" ]; then
            rm -rf $PYPLUSPLUS_DIST
        fi

        pushd $PYPLUSPLUS_BUILD
            # fix slice bug
            sed -i -e 's/m_start = std::max/\/\/m_start = std::max/g' build/lib*/pyplusplus/code_repository/indexing_suite/slice_header.py
            sed -i -e 's/m_stop = std::max/\/\/m_stop = std::max/g' build/lib*/pyplusplus/code_repository/indexing_suite/slice_header.py
            #patch -p1 < $BUILDSCRIPT_HOME/patches/pyplusplus-slice-fix.patch
        popd

        cp -rf $PYPLUSPLUS_BUILD/build/lib*/pyplusplus $PYPLUSPLUS_DIST
        pushd $PYPLUSPLUS_DIST

        popd
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

    rm -rf $TRIANGLE_BUILD
    mkBuildDIR $TRIANGLE_BUILD $TRIANGLE_SRC

    pushd $TRIANGLE_BUILD
        if [ "$SYSTEM" == "WIN" ]; then
            sed -i -e 's/-DLINUX/-DCPU86/g' makefile ;
            patch triangle.c -i $BUILDSCRIPT_HOME/patches/triangle-mingw-win64.patch
        elif [ "$SYSTEM" == "MAC" ]; then
            sed -i -e 's/-DLINUX//g' makefile ;
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

    getWITH_GIT $SUITESPARSE_URL $SUITESPARSE_SRC $SUITESPARSE_REV
    #getWITH_WGET $SUITESPARSE_URL $SUITESPARSE_SRC $SUITESPARSE_VER.tar.gz
    [ -d $SRC_DIR/SuiteSparse ] && mv $SRC_DIR/SuiteSparse $SUITESPARSE_SRC

    mkBuildDIR $SUITESPARSE_BUILD $SUITESPARSE_SRC 1

    pushd $SUITESPARSE_BUILD
        if [ "$SYSTEM" == "WIN" ]; then
            #patch -p1 < $BUILDSCRIPT_HOME/patches/SuiteSparse-4.4.1.patch
            echo "LIB = -lm" >> SuiteSparse_config/SuiteSparse_config.mk
            echo "CC = gcc" >> SuiteSparse_config/SuiteSparse_config.mk
            echo "CFLAGS=-std=c90" >> SuiteSparse_config/SuiteSparse_config.mk
            echo "BLAS = -L$TRIANGLE_DIST -lblas" >> SuiteSparse_config/SuiteSparse_config.mk
        elif [ "$SYSTEM" == "MAC" ]; then
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
                    CFLAGS='-std=c90' "$MAKE" -j$PARALLEL_BUILD library
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

slotAll(){
    buildBOOST
    buildLAPACK
    buildTRIANGLE
    buildSUITESPARSE
    buildCASTXML
    buildPYGCCXML
}

showHelp(){
    echo "boost | lapack | triangle | suitesparse | castxml | castxmlbin | pygccxml | all"
}

# script starts here
if [ -z "$TOOLSET" ]; then
    TOOLSET=none
fi
echo "****** Third party build via bash script for: $@ ******"
echo "env: TOOLSET set to: " $TOOLSET

if [ -n "$BOOST_VERSION" ]; then
    BOOST_VERSION=$BOOST_VERSION
else
    BOOST_VERSION=$BOOST_VERSION_DEFAULT
fi

if [ -n "$GIMLI_PREFIX" ]; then
    GIMLI_PREFIX=`readlink -m $GIMLI_PREFIX`
else
    GIMLI_PREFIX=`pwd`
fi

if [ -z "$PARALLEL_BUILD" ]; then
    PARALLEL_BUILD=1
fi
echo "Installing at " $GIMLI_PREFIX

CMAKE_BUILD_TYPE=Release

for arg in $@
do
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
    castxml)
        buildCASTXML;;
    castxmlbin)
        buildCASTXMLBIN;;
    pygccxml)
        buildPYGCCXML;;
    cppunit)
        buildCPPUNIT;;

    *)
        echo "Don't know what to do."
        showHelp;;
    esac
done
