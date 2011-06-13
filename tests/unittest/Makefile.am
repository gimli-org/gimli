TESTS = ./unittest.sh

#compile only when make check is called
check_PROGRAMS=\
	gimliUnitTest

INCLUDES = \
	-I $(top_srcdir)/src \
	$(CPPUNIT_CFLAGS)
 
gimliUnitTest_SOURCES = \
	unittest.sh \
	unitTest.cpp \
	testGIMLiMisc.h \
	testDataContainer.h \
	testVector.h \
	testShape.h \
	testGeometry.h \
	testFEM.h \
	testExternals.h 

# 	testOptionMap.h 

if HAVE_UNITTEST
#gimliUnitTest_LDADD = $(top_srcdir)/src/.libs/libgimli.la $(CPPUNIT_LIBS)
gimliUnitTest_LDADD = ../../src/.libs/libgimli.la $(CPPUNIT_LIBS)
else
gimliUnitTest_LDADD = 
endif

