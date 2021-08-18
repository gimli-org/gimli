#!/usr/bin/env bash

slotCMAKE(){
    cmake -G 'Unix Makefiles' ../gimli -DOpenBLAS_INCLUDE_DIR=C:/Users/ruc-bfe/bin/mingw64/mingw64/include/OpenBLAS
}
slotAll(){
    slotCMAKE
    slotPG
}
slotPG(){
    slotGimli
    #mv ../gimli/pygimli/core/_pygimli_.pyd ../gimli/pygimli/core/_pygimli_.pyd.back
    make pygimli J=4
}

slotGimli(){
    mv ../gimli/pygimli/core/libgimli.dll ../gimli/pygimli/core/libgimli.dll.back
    make -j4
    cp ./bin/libgimli.dll ../gimli/pygimli/core/libgimli.dll
}

slotRePG(){
    mv ../gimli/pygimli/core/_pygimli_.pyd ../gimli/pygimli/core/_pygimli_.pyd.back
    make -j4 pg

}
for arg in $@
do
    echo $arg
    case $arg in
    all)
        slotAll;;
    repg)
	slotRePG;;
    pg)
	slotPG;;
    g)
        slotGimli;;
    *)
        echo "Don't know what to do."
    esac
done
