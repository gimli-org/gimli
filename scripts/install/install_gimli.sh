#!/usr/bin/env bash

GIMLI_ROOT=$PWD/gimli
PARALLEL_BUILD=1

help(){
    echo "==============================================================="
    echo "Command line options:"
    echo "---------------------------------------------------------------"
    echo " add these option with out - at the end of the command line    "
    echo "---------------------------------------------------------------"
    echo "h"
    echo "      This help"
    echo "path=GIMLI_ROOT "
    echo "      Path for the gimli root dir. Default path=$GIMLI_ROOT"
    echo "j=PARALLEL_BUILD"
    echo "      Number of parallel compiler jobs. Default j=$PARALLEL_BUILD"
    exit
}

for i in "$@"; do
    case "$i" in
        *gimli_root=*|*path=*)
            GIMLI_ROOT=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
        ;;
        *j=*)
            PARALLEL_BUILD=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
        ;;
        *h*|*help*)
            help
        ;;
    esac
done

echo "Installing at "$GIMLI_ROOT
echo "Parallelize with j="$PARALLEL_BUILD

#SCRIPT_REPO='-Ls -x http://wwwproxy:8080 https://raw.githubusercontent.com/gimli-org/gimli/dev/scripts/install'
SCRIPT_REPO='-Ls https://raw.githubusercontent.com/gimli-org/gimli/dev/scripts/install'
GET="curl" 

#SCRIPT_REPO=$GIMLI_ROOT/gimli/scripts/install
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