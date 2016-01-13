#!/usr/bin/env bash

[ -z $GIMLI_ROOT ] && $(pwd)
[ -z $GIMLI_SOURCE ] && GIMLI_SOURCE=$GIMLI_ROOT/gimli
[ -z $GIMLI_BUILD ] && GIMLI_BUILD=$GIMLI_ROOT/build
[ -z $PARALLEL_BUILD ] && PARALLEL_BUILD=1
[ -z $PYTHON_MAJOR ] && PYTHON_MAJOR=3

CMAKE_GENERATOR='Unix Makefiles'

if [ -z $PYTHONSPECS ]; then

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
fi

buildGIMLI(){

    mkdir -p $GIMLI_BUILD

    pushd $GIMLI_BUILD

        [ -n $CASTXML ] && CASTXML="-DCASTER_EXECUTABLE=$CASTXML"

        cmake -G "$CMAKE_GENERATOR" $GIMLI_SOURCE $PYTHONSPECS $CASTXML

        make -j$PARALLEL_BUILD && make pygimli J=$PARALLEL_BUILD

    popd
}

buildGIMLI
