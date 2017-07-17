#!/usr/bin/env bash

ROOT=$PWD
[ -z $PARALLEL_BUILD ] && PARALLEL_BUILD=1
[ -z $PYTHON_MAJOR ] && PYTHON_MAJOR=3

CMAKE_GENERATOR='Unix Makefiles'

echo "CLEAN: ", $CLEAN

buildGIMLI(){
       
    echo "Building at:" $PWD
    mkdir -p build
    pushd build
        cmake -G "$CMAKE_GENERATOR" ../gimli -DBLAS_LIBRARIES=/mingw64/lib/libopenblas.a

        make -j$PARALLEL_BUILD && make pygimli J=$PARALLEL_BUILD
                 
    popd
}

buildGIMLI
