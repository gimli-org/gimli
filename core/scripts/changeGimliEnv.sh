GIMLIDIR=$1
python -c 'import pygimli as pg; print(pg.version())'

echo "Switching to " $1

export PATH=$GIMLIDIR/gimli/apps:$PATH
export LD_LIBRARY_PATH=$GIMLIDIR/build/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$GIMLIDIR/gimli:$PYTHONPATH

python -c 'import pygimli as pg; print(pg.version())'

# GPATH=$1
# LD=$GPATH/build/lib
# PY=$GPATH/gimli/python
# echo "LD:" $LD
# echo "PY:" $PY

# export LD_LIBRARY_PATH=$LD:$LD_LIBRARY_PATH
export PYTHONPATH=$PYTHONPATH:$HOME/src
