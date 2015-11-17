#!/usr/bin/env bash

# annoyed from repeating your passwd:
# git config credential.helper store
#
#

SRC=$pwd
PARALELL=2

if [ "$OSTYPE" == "msys" -o "$MSYSTEM" == "MINGW32" ]; then
    CMAKE_GENERATOR='MSYS Makefiles'
else
    CMAKE_GENERATOR='Unix Makefiles'
fi

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
        mkdir build
        pushd build
            cmake -G "$CMAKE_GENERATOR" ../gimli
            make -j$PARALELL && make pygimli J=$PARALELL 
            #python -c 'import pygimli as pg; pg.test()'
        popd
    popd
}

buildBert(){
    mkdir -p bert
    
    pushd bert
        if [ -d "bert" ]; then
            pushd bert
                git pull   
            popd  
        else
            git clone https://gitlab.com/resistivity-net/bert.git
        fi

        chmod +x bert/python/apps/*
        
        rm -rf build/
        mkdir build
        pushd build
            cmake -G "$CMAKE_GENERATOR" ../bert
            make -j$PARALELL && make pybert J=$PARALELL
        popd
    popd
}

slotAll(){
    buildGIMLI
    buildBert
}

help(){
    echo "bert | gimli | all"
}

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
    bert)
        buildBert;;
    gimli)
        buildGIMLI;;
    
    *) 
        echo "Don't know what to do."
        help;;
    esac
done