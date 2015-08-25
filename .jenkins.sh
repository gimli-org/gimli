# This script is for continuous integration using Jenkins (http://jenkins-ci.org/)
# It is called from the parent directory, i.e. bash -xe trunk/.jenkins.sh

echo "Starting automatic build #$BUILD_NUMBER on" `date`
start=$(date +"%s")
uname -a
gcc --version
python --version
cmake --version
python -c "import numpy; print(numpy.__version__)"

# just do this if something is wrong with the thirdparty sources
# rm -rf thirdParty/src

# Main build
#rm -rf build # Uncomment for clean build (expensive, but necessary sometimes)
mkdir -p build
cd build
cmake ../trunk
make -j 16 gimli
make pygimli J=12

# Test gimli
make check
./bin/gimliUnitTest

# Test pygimli
export PYTHONPATH=`pwd`/../trunk/python:$PYTHONPATH
python -c "import pygimli; print(pygimli.__version__)"
python -c "import pygimli; pygimli.test(onlydoctests=False, coverage=True)"

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
