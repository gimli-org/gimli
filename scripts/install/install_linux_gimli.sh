#!/usr/bin/env bash

ROOT=$PWD
PARALELL=2

CMAKE_GENERATOR='Unix Makefiles'

case "$(grep "ID=" /etc/os-release)" in
    *"gentoo"*)
        echo "Gentoo system found"
        PYTHONSPECS=''
    ;;
    *"debian"*)
        echo "Debian system found"
        PYTHONSPECS='-DPYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libpython3.4m.so.1.0
                    -DBoost_PYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libboost_python-py34.so 
                    -DPYTHON_EXECUTABLE=/usr/bin/python3'
    ;;
    *"arch"*)
        echo "Arch Linux system found"
        PYTHONSPECS='-DBoost_PYTHON_LIBRARY=/usr/lib64/libboost_python3.so'
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

        make -j$PARALELL && make pygimli J=$PARALELL 
        echo ""
        echo ""
        echo "============================================================================"
        echo "---  try some basic test: calling pygimli once------------------------------"
        export PYTHONPATH=$ROOT/gimli/python:$PYTHONPATH
        python -c 'import pygimli as pg; print("pygimli version:", pg.__version__)'
        echo "--- ------------------------------------------------------------------------"
        echo "export PYTHONPATH=$ROOT/gimli/python:\$PYTHONPATH" > $ROOT/.bash_hint_pygimli
        echo "export PATH=$ROOT/gimli/python/apps:\$PATH" >> $ROOT/.bash_hint_pygimli
           
    popd
}


buildGIMLI
