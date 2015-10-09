# This script is for continuous integration using Jenkins (http://jenkins-ci.org/)
# It is called from the parent directory, i.e. bash -xe trunk/.jenkins.sh

echo "Starting automatic build #$BUILD_NUMBER on" `date`
start=$(date +"%s")

# Show last change to repo in build log
echo `git --git-dir trunk/.git log -1 --pretty="Last change by %cn (%h): %B"`

# Show system information
lsb_release -d
uname -a
gcc --version
cmake --version
python3 --version
python3 -c "import numpy; print(numpy.__version__)"
python3 -c "import matplotlib; print(matplotlib.__version__)"

################
#  Main build  #
################

# just do this if something is wrong with the thirdparty sources
#rm -rf thirdParty/src
#rm -rf build # Uncomment for clean build (expensive, but necessary sometimes)
rm -f build/build_tests.html # remove old test report
rm -f build/CMakeCache.txt # clean old cache

mkdir -p build
cd build
cmake ../trunk \
    -DCASTER=castxml \
    -DPYVERSION=3 \
    -DPYTHON_EXECUTABLE=/usr/bin/python3 \
    -DPYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libpython3.4m.so \
    -DBoost_PYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libboost_python-py34.so

make -j 16 gimli
make pygimli J=12

#############################
#  Testing & documentation  #
#############################

# Test gimli
make check
./bin/gimliUnitTest

# Test pygimli
export PYTHONPATH=`pwd`/../trunk/python:$PYTHONPATH
python3 << END
import pygimli as pg
print(pg.__version__)
pg.test(onlydoctests=False, htmlreport='build_tests.html')
END

# Build documentation
export PATH=/opt/texbin:$PATH # for building pdf
export PUBLISH="True" # for correct PATH settings in sidebar gallery
export PATH=`pwd`/../trunk/python/apps:$PATH
chmod +x ../trunk/python/apps/*
make doc # = doxygen, sphinxapi, sphinxpdf, sphinxhtml
end=$(date +"%s")
echo "Ending automatic build #$BUILD_NUMBER".
diff=$(($end-$start))
echo "$(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."

# If this script fails, a mail with the log file will be send to mail@pygimli.org.
# If it succeeds, the documentation will be uploaded to pygimli.org.
# In any case, a badge icon of the current status and the log file will be uploaded.
