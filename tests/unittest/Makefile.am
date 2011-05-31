bin_PROGRAMS=\
	gimliUnitTest

TESTS = unittest.sh

#compile only when make check is called
check_PROGRAMS=\
	gimliUnitTest

INCLUDES = \
        -I $(top_srcdir)/src
 
gimliUnitTest_SOURCES = \
	unittest.sh \
	unitTest.cpp \
	testGIMLiMisc.h 

# 	testOptionMap.h \
# 	testShape.h \
# 	testVector.h \
# 	testGeometry.h \
# 	testFEM.h 

gimliUnitTest_LDADD = ../../src/.libs/libgimli.la -lcppunit
