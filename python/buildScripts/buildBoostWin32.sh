BOOST_VERSION=1.53.0

if [ $# -eq 0 ]; then
	prefix=`pwd`
else 
	prefix=`readlink -m $1`
fi

echo "Installing boost at " $prefix

BOOST_NAME=boost_${BOOST_VERSION//./_}
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

if ( wget --version ); then 
	echo "########### wget found: ok" ; 
else 
	echo "########### Installing wget" ; 
	mingw-get.exe install msys-wget ; 
fi; 

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

BOOST_SRC_DIR=$BOOST_SRC/$BOOST_NAME
GCCVER=`gcc -dumpmachine`-`gcc -dumpversion`
GCCARCH=`gcc -dumpmachine`

if [ "$GCCARCH" == "mingw32" ]; then
	ADRESSMODEL=32
else
	ADRESSMODEL=64
fi

echo 'Found gcc ' $GCCVER 
WITHPYTHON='--with-python'
arch=`python -c 'import platform; print platform.architecture()[0]'`

if [ "$arch" == "64bit" ]; then
	if [ "$ADRESSMODEL" == "32" ]; then
		echo "32bit GCC but 64bit python .. omitting boostpython."
		WITHPYTHON=''
	fi
else
	if [ "$ADRESSMODEL" == "64" ]; then
		echo "64bit GCC but 32bit python .. omitting boostpython."
		WITHPYTHON=''
	fi
fi

if [ ! -d $BOOST_SRC_DIR ]; then
     wget -nc -nd http://sourceforge.net/projects/boost/files/boost/$BOOST_VERSION/$BOOST_NAME'.tar.gz'
     tar -xzvf $BOOST_NAME'.tar.gz'
fi

pushd $BOOST_SRC_DIR
	
    # ./tools/build/v2/tools/python.jam:486 (comment out line to disable quotation adding)
    sed -i 's/python-cmd = \\"\$(python-cmd)\\" ;/# python-cmd = \\"\$(python-cmd)\\" ;/' ./tools/build/v2/tools/python.jam

    # is this still necessary ??
    # edit ./tools/build/v2/tools/python-config.jam:12 (add 2.7 2.6 2.5) but not necessary
    
    DISTDIR=$BOOST_SRC_DIR/build/boost_$BOOST_VERSION-$GCCVER-$ADRESSMODEL

	echo Calling from $OLDDIR
	echo Installing at $DISTDIR

	if [ ! -f ./b2.exe ]; then
		 cmd /c "bootstrap.bat mingw --with-python-root=$PYTHON_ROOT "
	fi

	./b2.exe toolset=gcc variant=release link=shared threading=multi address-model=$ADRESSMODEL install \
		--prefix=$DISTDIR \
		--layout=tagged \
		$WITHPYTHON \
		--with-system \
		--with-thread \
		--with-date_time \
		--with-chrono \
		--with-regex 
		#libraries=python,system,thread,regex 
		
	
	#sed -e 's/ using mingw / using gcc /' project-config.jam > tmp
	#mv tmp project-config.jam
	
	#exit
	#LDFLAGS='-static-libgcc -static-libstdc++' 
	#./b2 toolset=gcc install -d+2 --prefix=$DISTDIR --layout=tagged \
	#		address-model=$ADRESSMODEL variant=release link=shared \
	#		threading=multi

    echo "copying into new boost dir", ../boost
	mkdir -p ../boost
	cp -r $DISTDIR/include ../boost
	cp -r $DISTDIR/lib ../boost
popd
