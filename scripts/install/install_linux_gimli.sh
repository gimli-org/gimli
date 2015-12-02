#!/usr/bin/env bash

ROOT=$PWD

[ -z $PARALLEL_BUILD ] && PARALLEL_BUILD=1
[ -z $PYTHON_MAJOR ] && PYTHON_MAJOR=3

CMAKE_GENERATOR='Unix Makefiles'

case "$(grep "ID=" /etc/os-release)" in
    *"gentoo"*)
        echo "Gentoo system found"
        PYTHONSPECS=''
    ;;
    *"debian"*)
        echo "Debian system found"
        if [ $PYTHON_MAJOR -eq 3 ] ; then
            PYTHONSPECS='-DPYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libpython3.4m.so.1.0
                        -DBoost_PYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libboost_python-py34.so
                        -DPYTHON_EXECUTABLE=/usr/bin/python3'
        else
            PYTHONSPECS=''
        fi

    ;;
    *"arch"*)
        echo "Arch Linux system found"
        PYTHONSPECS='-DBoost_PYTHON_LIBRARY=/usr/lib64/libboost_python3.so'
    ;;
    *"ubuntu"*)
        echo "Ubuntu Linux system found"
        PYTHONSPECS='-DPYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libpython3.4m.so.1.0
                     -DBoost_PYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libboost_python-py34.so
                     -DPYTHON_EXECUTABLE=/usr/bin/python3'
    ;;
    *)
        echo $(grep "ID=" /etc/os-release) "system found: trying defaults"
        PYTHONSPECS=''
    ;;
esac

buildGIMLI(){
    if [ -d "gimli" ]; then
        pushd gimli
            git pull
        popd
    else
        git clone https://github.com/gimli-org/gimli.git
    fi

    [ -n "$BRANCH" ] && git checkout $BRANCH

    chmod +x gimli/python/apps/*

    rm -rf build/
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
