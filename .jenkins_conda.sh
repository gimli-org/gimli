# This script is for continuous integration using Jenkins (http://jenkins-ci.org/)
# It is called from the parent directory, i.e. bash -xe trunk/.jenkins.sh

echo "Starting automatic build #$BUILD_NUMBER on" `date`
start=$(date +"%s")

# Show last change to repo in build log
echo `git --git-dir trunk/.git log -1 --pretty="Last change by %cn (%h): %B"`

# Show system information
lsb_release -d
uname -a

conda env update
conda activate pgdev

which python
OMP_THREAD_LIMIT=1 python -c "import pygimli; pygimli.test(show=False, abort=True, htmlreport=\"build_tests.html\")"

# Build documentation

# Setup for 3D visualizations
# ------------------------------------------------
export DISPLAY=:99.0
export PYVISTA_OFF_SCREEN=True
# ------------------------------------------------

cd doc
make clean
make html

end=$(date +"%s")
echo "Ending automatic build #$BUILD_NUMBER".
diff=$(($end-$start))
echo "$(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."

# If this script fails, a mail with the log file will be send to mail@pygimli.org.
