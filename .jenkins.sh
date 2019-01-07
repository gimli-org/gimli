# This script is for continuous integration using Jenkins (http://jenkins-ci.org/)
# It is called from the parent directory, i.e. bash -xe trunk/.jenkins.sh

echo "Starting automatic build #$BUILD_NUMBER on" `date`
start=$(date +"%s")

# Show last change to repo in build log
echo `git --git-dir trunk/.git log -1 --pretty="Last change by %cn (%h): %B"`

# link python to python3
ln -sf /usr/bin/python3 python
export PATH=`pwd`:$PATH

# Show system information
lsb_release -d
uname -a
gcc --version
cmake --version
python --version
python -c "import numpy; print(numpy.__version__)"
python -c "import matplotlib; print(matplotlib.__version__)"

# Check if core was changed
core_update=$(git --git-dir=trunk/.git diff-tree -r $GIT_COMMIT | grep -c src || true)

# Set this to 1 if you want clean build (also of dependencies)
export CLEAN=0

export GIMLI_NUM_THREADS=4

################
#  Main build  #
################

# just do this if something is wrong with the thirdparty sources
#rm -rf thirdParty/src
#rm -rf build # Uncomment for clean build (expensive, but necessary sometimes)
#rm -f build/build_tests.html # remove old test report
#rm -f build/CMakeCache.txt # clean old cache

mkdir -p build
cd build

if [ ! -f build/CMakeCache.txt ]; then
    # Always rebuild core when Cmake cache does not exist
    core_update=2
fi

if [[ $core_update -ge 1 ]]; then
  echo "# Core changes detected. #"
  cmake ../trunk \
     -DPYVERSION=3 \
     -DPYTHON_EXECUTABLE=/usr/bin/python3 \
     -DPYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libpython3.5m.so \
     -DBoost_PYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libboost_python-py35.so
  make -j 8 gimli
  make pygimli J=4
else
  echo "# No core changes detected. #"
fi

#############################
#  Testing & documentation  #
#############################

# Test pygimli
export PYTHONPATH=`pwd`/../trunk/python:$PYTHONPATH

OMP_THREAD_LIMIT=1 python -c "import pygimli; pygimli.test(show=False, abort=True, htmlreport=\"build_tests.html\")"

# Build documentation
export PUBLISH="True" # for correct PATH settings in sidebar gallery
export PATH=`pwd`/../trunk/python/apps:$PATH
chmod +x ../trunk/python/apps/*

make clean-gallery
make doc # = doxygen, sphinxapi, sphinxpdf, sphinxhtml
end=$(date +"%s")
echo "Ending automatic build #$BUILD_NUMBER".
diff=$(($end-$start))
echo "$(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."

# If this script fails, a mail with the log file will be send to mail@pygimli.org.
# If it succeeds, the documentation will be uploaded to pygimli.org.
# In any case, a badge icon of the current status and the log file will be uploaded.
