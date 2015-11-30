#!/bin/env bash

SCRIPT_REPO=https://raw.githubusercontent.com/gimli-org/gimli/dev/scripts/install

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
    curl -s $SCRIPT_REPO/install_$SYSTEM'_winpython.sh' | bash
    [ -f .bash_hint_python ] && source .bash_hint_python        
    
elif [ $(uname -o) == "GNU/Linux" ]; then
    SYSTEM='linux'
    echo "Determining system ... LINUX system found"
fi
echo "------------------------------------------"

echo "=========================================="
echo "Install system prerequisites for" $SYSTEM
echo "------------------------------------------"
curl -s $SCRIPT_REPO/install_$SYSTEM'_prereqs.sh' | bash

echo "=========================================="
echo "Install gimli for" $SYSTEM
echo "------------------------------------------"
curl -s $SCRIPT_REPO/install_$SYSTEM'_gimli.sh' | bash

echo ""
echo "=========================================="
echo "set the followng setting to use pygimli, either local per session or permanently in your $HOME/.bashrc"
echo "------------------------------------------"
echo ""
[ -f .bash_hint_python ] && cat .bash_hint_python 
[ -f .bash_hint_pygimli ] && cat .bash_hint_pygimli
echo ""
echo "------------------------------------------"

#https://raw.githubusercontent.com/gimli-org/gimli/master/scripts/install/install**