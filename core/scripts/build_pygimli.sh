#!/usr/bin/env sh
# This script checks out pygimli from git and builds it. Afterwards, Python
# source and binary dist tar.gz files are created, as well as a .whl and a .deb
# file

# Refer to https://www.pygimli.org/compilation.html#useful-cmake-settings
# for additional cmake settings

git clone https://github.com/gimli-org/gimli.git

mkdir build
cd build
cmake ../gimli
make -j 2 gimli
make -j 2 apps
make -j 2 pygimli

cd ../gimli
# source and binary tar.gz files
python3 setup.py sdist
python3 setup.py bdist
# build a wheel that can be installed using pip
pip wheel .

outdir="/HOST/pg_artifacts"
test -d "${outdir}" || mkdir "${outdir}"
cp -r dist "${outdir}/"
cp *.whl "${outdir}"

# build a quick and dirty .deb file
export PYBUILD_DISABLE=test
python3 setup.py --command-packages=stdeb.command bdist_deb
cp -r deb_dist "${outdir}/"

# make sure we can read/write/delete the output data on the host
chmod -R 0777 "${outdir}"
