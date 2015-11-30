#!/usr/bin/env bash

if [ $# -lt 1 ]; then
    GIMLI_ROOT=$PWD/gimli
else
    GIMLI_ROOT=$1
fi
# echo $GIMLI_ROOT
# [ -z "$GIMLI_ROOT" ] && GIMLI_ROOT=$PWD/gimli

echo " "
echo "GIMLI_ROOT="$GIMLI_ROOT

SCRIPT_REPO='-Ls https://raw.githubusercontent.com/gimli-org/gimli/dev/scripts/install'
GET="curl" 

#SCRIPT_REPO=$PREFIX/gimli/scripts/install
#GET="cat"

PYTHON_MAJOR=3

echo "=========================================="
if [ $(uname -o) == "Msys" ]; then
    if [ $(uname -m) == "x86_64" ]; then
        SYSTEM='win64'
        echo "Determining system ... Msys WIN64 system found"
        # just run this for the initial call .. after a mingw-shell restart these setting is automatic
        export PATH=/mingw64/bin:$PATH
    else
        SYSTEM='win32'
        echo "Determining system ... Msys WIN32 system found"
        # just run this for the initial call .. after a mingw-shell restart these setting is automatic
        export PATH=/mingw32/bin:$PATH
    fi
    # check for winpython
    "$GET" $SCRIPT_REPO/install_$SYSTEM'_winpython.sh' | bash -s $PYTHON_MAJOR
    [ -f $GIMLI_ROOT/.bash_hint_python ] && source $GIMLI_ROOT/.bash_hint_python        
    
elif [ $(uname -o) == "GNU/Linux" ]; then
    SYSTEM='linux'
    echo "Determining system ... LINUX system found"
fi
echo "------------------------------------------"

echo "=========================================="
echo "Install system prerequisites for" $SYSTEM
echo "------------------------------------------"
"$GET" $SCRIPT_REPO/install_$SYSTEM'_prereqs.sh' | bash 

echo "=========================================="
echo "Install gimli for" $SYSTEM
echo "------------------------------------------"

mkdir -p $GIMLI_ROOT
pushd $GIMLI_ROOT
    "$GET" $SCRIPT_REPO/install_$SYSTEM'_gimli.sh' | bash
popd

echo ""
echo "=========================================="
echo "set the followng setting to use pygimli, either local per session or permanently in your $HOME/.bashrc"
echo "------------------------------------------"
echo ""
[ -f $GIMLI_ROOT./bash_hint_python ] && cat $GIMLI_ROOT/.bash_hint_python 
[ -f $GIMLI_ROOT/.bash_hint_pygimli ] && cat $GIMLI_ROOT/.bash_hint_pygimli
echo ""
echo "------------------------------------------"

#https://raw.githubusercontent.com/gimli-org/gimli/master/scripts/install/install**