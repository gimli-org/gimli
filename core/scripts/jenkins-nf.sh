#!/bin/bash
if [[ $- == *i* ]]; then
    # run with . *.sh
    echo "Running on interactive mode. Aliases will work"
else
    # run with bash *.sh
    echo "Running on NON interactive mode. Aliases will not work. Setting then now."
    shopt -s expand_aliases ## else aliases will be ignored in this bash session
    alias python='python3'
fi

# This script is for continuous integration using Jenkins (http://jenkins-ci.org/)
# It is called for the newfea branch from the jenkins workspace, i.e., 
# bash -xe ./gimli/core/scripts/jenkins-nf.sh

#rm -rf /var/lib/jenkins/.cache/pygimli # We need a better cleaning process
GIMLIROOT=$pwd
GIMLISRC=$GIMLIROOT/gimli

echo "Starting automatic build #$BUILD_NUMBER on" `date`
start=$(date +"%s")

# Show last change to repo in build log
echo `git --git-dir $GIMLISRC/.git log -1 --pretty="Last change by %cn (%h): %B"`

# Show system information
lsb_release -d
uname -a
gcc --version
cmake --version
python --version
python -c "import numpy; print('Numpy:', numpy.__version__)"
python -c "import matplotlib; print('MPL:', matplotlib.__version__)"

# Check if core was changed
core_update=$(git --git-dir=$GIMLISRC/.git diff --name-only $GIT_PREVIOUS_COMMIT $GIT_COMMIT | grep -c core/src || true)

# Set this to 1 if you want clean build (also of dependencies)
export CLEAN=0
export GIMLI_NUM_THREADS=$((`nproc --all` - 4))

################
#  Main build  #
################

function build(){
    PYVERSION=$1
    
    pushd $GIMLIROOT
        PYPLATFORM=cp$PYVERSION-cp$PYVERSION
        BUILDENV=venv_$PYPLATFORM
        BUILDDIR=build_$PYPLATFORM
        $PYTHON -m venv $BUILDENV
        source $BUILDENV/bin/activate

        # build minimum
        python -m pip install --upgrade pip
        python -m pip install numpy pygccxml pyplusplus 
        # build wheel needs
        python -m pip install wheel auditwheel twine 

        python3 -m pip install -r $GIMLISRC/core/scripts/requirements.txt

        mkdir -p $BUILDDIR
        pushd $BUILDDIR

        if [ $core_update -ge 1 ]; then
            echo "# Core changes detected. #"
            rm -rf CMakeCache.txt

            cmake $GIMLISRC

            make -j $GIMLI_NUM_THREADS 
            make pygimli J=$GIMLI_NUM_THREADS
        
        else
            echo "# No core changes detected. #"
        fi

        make whlTest            

        popd        
    popd
}

#build

#############################
#  Testing & documentation  #
#############################

# # Setup for 3D visualizations
# # ------------------------------------------------
# export DISPLAY=:99.0
# export PYVISTA_OFF_SCREEN=true
# # ------------------------------------------------

# # Test pygimli
# export PYTHONPATH=`pwd`/../trunk:$PYTHONPATH
# OMP_THREAD_LIMIT=4 python -c "import pygimli; print(pygimli.Report()); pygimli.test(show=False, abort=True, htmlreport=\"build_tests.html\", devTests=True)"

# # Build documentation
# make clean-gallery
# make doc # = doxygen, sphinxapi, sphinxpdf, sphinxhtml
# end=$(date +"%s")
# echo "Ending automatic build #$BUILD_NUMBER".
# diff=$(($end-$start))
# echo "$(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."

# # If this script fails, a mail with the log file will be send to mail@pygimli.org.
# # If it succeeds, the documentation will be uploaded to pygimli.org.
# # In any case, a badge icon of the current status and the log file will be uploaded.
