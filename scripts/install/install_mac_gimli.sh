#!/usr/bin/env bash

ROOT=$PWD
[ -z $PARALLEL_BUILD ] && PARALLEL_BUILD=1
[ -z $PYTHON_MAJOR ] && PYTHON_MAJOR=3

CMAKE_GENERATOR='Unix Makefiles'

PYTHONSPECS=''

buildGIMLI(){
    if [ -d "gimli" ]; then
        pushd gimli
            git pull
        popd
    else
        git clone https://github.com/gimli-org/gimli.git
    fi

    echo $(pwd)
    pushd gimli
        echo "switching to branch: " $BRANCH
        [ -n "$BRANCH" ] && git checkout $BRANCH
    popd
    

    chmod +x gimli/python/apps/*
    
    if [ $UPDATE_ONLY -eq 0 ] ; then
        rm -rf build/
    fi
    mkdir -p build

    pushd build
        cmake -G "$CMAKE_GENERATOR" ../gimli $PYTHONSPECS

        make -j$PARALLEL_BUILD && make pygimli J=$PARALLEL_BUILD
        echo ""
        echo ""
        echo "============================================================================"
        echo "------------------------  TEST pyGIMLi installation ------------------------"
        export PYTHONPATH=$ROOT/gimli/python:$PYTHONPATH
        python -c 'import pygimli as pg; print("pygimli version:", pg.__version__)'
        if [ -x "$(command -v pytest)" ]; then
            python -c 'import pygimli as pg; pg.test()'
        fi
        echo "--- ------------------------------------------------------------------------"
        echo "export PYTHONPATH=$ROOT/gimli/python:\$PYTHONPATH" > $ROOT/.bash_hint_pygimli
        echo "export PATH=$ROOT/gimli/python/apps:\$PATH" >> $ROOT/.bash_hint_pygimli

    popd
}


buildGIMLI
