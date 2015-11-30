#!/usr/bin/env bash

SRC=$pwd
PARALELL=2

CMAKE_GENERATOR='MSYS Makefiles'

buildGIMLI(){
    mkdir -p gimli
    
    pushd gimli
        if [ -d "gimli" ]; then
            pushd gimli
                git pull   
            popd  
        else
            git clone https://github.com/gimli-org/gimli.git
        fi
            
        chmod +x gimli/python/apps/*
        
        rm -rf build/
        mkdir -p build
        pushd build
            cmake -G "$CMAKE_GENERATOR" ../gimli -DBLAS_LIBRARIES=/mingw64/lib/libopenblas.a

            make -j$PARALELL && make pygimli J=$PARALELL 
            #python -c 'import pygimli as pg; pg.test()'

            echo "PYTHONPATH=$SRC/gimli/gimli/python" > .bash_hint_pygimli
            echo "PATH=\$PATH:$SRC/gimli/gimli/apps" >> .bash_hint_pygimli
        popd
    popd
}

buildGIMLI