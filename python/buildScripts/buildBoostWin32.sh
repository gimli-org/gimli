BOOST_VERSION=1_49_0

if [ $# -eq 0 ]; then
	prefix=`pwd`
else 
	prefix=`readlink -m $1`
fi

echo "Installing boost at " $prefix

BOOST_SRC=$prefix

if ( (python --version) );then
	PYTHON_ROOT=`python -c "import sys; print sys.prefix.replace(\"\\\\\\\\\",\"/\")"`
	echo "... found at $PYTHON_ROOT, good"
else
	echo "need > python2.7 installation"
	echo "get one from http://www.python.org/"
	echo "if allready done, ensure python27 installation directory is in your PATH"
	exit
fi

# fixing python2.7 win64 installation -- there is a file missing
#if [ ! -f $PYTHON_ROOT/libs/libpython27.a ]; then
#	echo missing $PYTHON_ROOT/libs/libpython27.a
#	echo "try creating them"
	# http://wiki.cython.org/InstallingOnWindows?action=AttachFile&do=get&target=python27.def
	# pexports python24.dll > python24.def 
	# mingw-get.exe install pexports
	# pexports python27.dll > python27.def 
	# dlltool -D python27.dll -d python27.def -l libpython27.a
#	dlltool --dllname python27.dll --def libgimli/trunk/external/patches/python27.def --output-lib $PYTHON_ROOT/libs/libpython27.a
#fi

BOOST_SRC_DIR=$BOOST_SRC/boost_$BOOST_VERSION
GCCVER=mingw-`gcc -v 2>&1 | tail -n1 | cut -d' ' -f2`
if [ "$GCCVER" = "mingw-version" ]; then
    GCCVER=mingw-`gcc -v 2>&1 | tail -n1 | cut -d' ' -f3`
fi

arch=`python -c 'import platform; print platform.architecture()[0]'`
if [ "$arch" == "64bit" ]; then
	ADRESSMODEL=64
else
	ADRESSMODEL=32
fi

if [ ! -d $BOOST_SRC_DIR ]; then
     wget -nc -nd http://sourceforge.net/projects/boost/files/boost/1.48.0/boost_1_48_0.tar.gz
     tar -xzvf boost_1_48_0.tar.gz 
fi

pushd $BOOST_SRC_DIR
	
    # ./tools/build/v2/tools/python.jam:486 (comment out line to disable quotation adding)
    sed -i 's/python-cmd = \\"\$(python-cmd)\\" ;/# python-cmd = \\"\$(python-cmd)\\" ;/' ./tools/build/v2/tools/python.jam

    # is this still necessary ??
    # edit ./tools/build/v2/tools/python-config.jam:12 (add 2.7 2.6 2.5) but not necessary
    
    DISTDIR=$BOOST_SRC_DIR/boost_$BOOST_VERSION-$GCCVER-$ADRESSMODEL

	echo Calling from $OLDDIR
	echo Installing at $DISTDIR

	if [ ! -f ./bjam.exe ]; then
		./bootstrap.sh --with-toolset=mingw 
	fi

	LDFLAGS='-static-libgcc -static-libstdc++' ./bootstrap.sh --prefix=$DISTDIR --with-bjam=./bjam.exe --with-toolset=gcc \
		--with-python-root=$PYTHON_ROOT --with-libraries=python,system,thread,regex
		
	LDFLAGS='-static-libgcc -static-libstdc++' ./b2 install -d+2 --prefix=$DISTDIR --layout=tagged \
			address-model=$ADRESSMODEL variant=release link=shared \
			threading=multi

    echo "copying into new boost dir", ../boost
	mkdir -p ../boost
	cp -r $DISTDIR/include ../boost
	cp -r $DISTDIR/lib ../boost
popd
