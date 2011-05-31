rm -rf autom4te.cache

#autoreconf -i -s -f &&
#echo "- libtoolize"; libtoolize --force --copy &&
echo "- autoreconf"; autoreconf --force --install &&
./configure


# Fire up autotools
#echo "- libtoolize"; libtoolize --force --copy 
#echo "- aclocal."; aclocal #-I Scripts/m4 $ACLOCAL_FLAGS 
#echo "- autoconf." ; autoconf 
#echo "- autoheader."; autoheader 
#echo "- automake."; automake --include-deps --add-missing --foreign --copy 

