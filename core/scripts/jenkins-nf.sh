#!/bin/bash
# This script is for continuous integration using Jenkins (http://jenkins-ci.org/)
# It is called for the newfea branch from the jenkins workspace, i.e., 
# bash -xe ./gimli/core/scripts/jenkins-nf.sh
if [[ $- == *i* ]]; then
    # run with . *.sh
    echo "Running on interactive mode. Aliases will work"
else
    # run with bash *.sh
    echo "Running on NON interactive mode. Aliases will not work. Setting then now."
    shopt -s expand_aliases ## else aliases will be ignored in this bash session
    alias python='python3'
fi


#rm -rf /var/lib/jenkins/.cache/pygimli # We need a better cleaning process
echo "Starting automatic build #$BUILD_NUMBER on" `date` 

echo "JOB_NAME=$JOB_NAME"
echo "JOB_BASE_NAME=$JOB_BASE_NAME"
echo "JENKINS_HOME=$JENKINS_HOME"
echo "WORKSPACE=$WORKSPACE"
echo "BUILD_TAG=$BUILD_TAG"

#rm -rf /var/lib/jenkins/.cache/pygimli # We need a better cleaning process
AGENTS_ROOT=$WORKSPACE/../
GIMLIROOT=$WORKSPACE
GIMLISRC=$GIMLIROOT/gimli

echo "GIMLIROOT=$GIMLIROOT"
echo "GIMLISRC=$GIMLISRC"


echo "Starting automatic build #$BUILD_NUMBER on" `date` "ROOT: $GIMLIROOT"
start=$(date +"%s")

# Show system information
lsb_release -d
uname -a
python --version

PYVERSION=$1
    
if [ -z $PYVERSION ]; then
    PYVERSION=$(python -c 'import sys; print(f"{sys.version_info.major}{sys.version_info.minor}")')
    echo "building for python: $PYVERSION"
else
    echo "building for python: $PYVERSION (forced by setting PYVERSION)"
fi

# Show last change to repo in build log
echo `git --git-dir $GIMLISRC/.git log -1 --pretty="Last change by %cn (%h): %B"`

# Check if core was changed
core_update=$(git --git-dir=$GIMLISRC/.git diff --name-only $GIT_PREVIOUS_COMMIT $GIT_COMMIT | grep -c core/src || true)

# Set this to 1 if you want clean build (also of dependencies)
export CLEAN=0
export GIMLI_NUM_THREADS=$((`nproc --all` - 4))

function build(){

    PYVERSION=$1
    PYPLATFORM=cp$PYVERSION-cp$PYVERSION

    GIMLIBLD=$GIMLIROOT/build_$PYPLATFORM
    BUILDENV=$GIMLIROOT/build_$PYPLATFORM'_venv'
    TESTVENV=$GIMLIROOT/test_$PYPLATFORM'_venv'

    pushd $GIMLIROOT
        
        python -m venv $BUILDENV
        source $BUILDENV/bin/activate

        echo "Updating pip ..."
        python -m pip install --upgrade pip
        echo "Installing or update build venv ..."
        python -m pip install --quiet -r $GIMLISRC/core/scripts/requirements.txt

        mkdir -p $GIMLIBLD
        pushd $GIMLIBLD

            if [ $core_update -ge 1 ]; then
                echo "# Core changes detected. #"
                rm -rf CMakeCache.txt

            else
                echo "# No core changes detected. #"
            fi

            cmake $GIMLISRC -DTESTVENV=$TESTVENV
            make -j $GIMLI_NUM_THREADS 
            make pygimli J=$GIMLI_NUM_THREADS
            make whlTest

        popd        
    popd
    
}


build $PYVERSION
#test
