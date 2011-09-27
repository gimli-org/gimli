BOOST_VERSION=1_47_0

if [ $# -eq 0 ]; then
	prefix=`pwd`
else 
	prefix=`readlink -m $1`
fi

echo "Installing boost at " $prefix

BOOST_SRC=$prefix
#BOOST_SRC=/c/Home

if ( (python --version) );then
	PYTHON_ROOT=`python -c "import sys; print sys.prefix.replace(\"\\\\\\\\\",\"/\")"`
	echo "... found at $PYTHON_ROOT, good"
else
	echo "need > python2.6 installation"
	echo "get one from http://www.python.org/"
	echo "if allready done, ensure python26 installation directory is in your PATH"
	exit
fi

BOOST_SRC_DIR=$BOOST_SRC/boost_$BOOST_VERSION
GCCVER=mingw-`gcc -v 2>&1 | tail -n1 | cut -d' ' -f3`

pushd $BOOST_SRC_DIR

	DISTDIR=$BOOST_SRC_DIR/boost_$BOOST_VERSION-$GCCVER

	echo calling from $OLDDIR
	echo Installing at $DISTDIR

	if [ ! -f ./bjam.exe ]; then
		./bootstrap.sh --with-toolset=mingw 
	fi

	# if you experience bjam complains something like founding no python
	# edit ./tools/build/v2/tools/python.jam:486 (comment out line to disable quotation adding)
	# edit ./tools/build/v2/tools/python-config.jam:12 (add 2.7 2.6 2.5) but not necessary

	./bootstrap.sh --prefix=$DISTDIR --with-bjam=./bjam.exe --with-toolset=gcc \
		--with-python-root=$PYTHON_ROOT --with-libraries=python,system,thread,regex
	
	./b2 --prefix=$DISTDIR --layout=tagged --build-type=complete --variant=release install

	cp -r $DISTDIR ../boost
popd
