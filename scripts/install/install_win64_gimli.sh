#!/usr/bin/env bash

ROOT=$PWD
[ -z $PARALLEL_BUILD ] && PARALLEL_BUILD=1
[ -z $PYTHON_MAJOR ] && PYTHON_MAJOR=3

CMAKE_GENERATOR='MSYS Makefiles'

buildGIMLI(){
        
    mkdir -p build
    pushd build
        cmake -G "$CMAKE_GENERATOR" ../gimli -DBLAS_LIBRARIES=/mingw64/lib/libopenblas.a

        make -j$PARALELL && make pygimli J=$PARALELL 
        echo ""
        echo ""
        echo "============================================================================"
        echo "---  try some basic test: calling pygimli once------------------------------"
        export PYTHONPATH=$PYTHONPATH:$ROOT/gimli/python
        python -c 'import pygimli as pg; print("pygimli version:", pg.__version__)'
        echo "--- ------------------------------------------------------------------------"
        echo "export PYTHONPATH=\$PYTHONPATH:$ROOT/gimli/python" > $ROOT/.bash_hint_pygimli
        echo "export PATH=\$PATH:$ROOT/gimli/python/apps" >> $ROOT/.bash_hint_pygimli
           
    popd
}

buildGIMLI
