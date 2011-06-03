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
	testGIMLiMisc.h 

# 	testOptionMap.h \
# 	testShape.h \
# 	testVector.h \
# 	testGeometry.h \
# 	testFEM.h 

if HAVE_UNITTEST
gimliUnitTest_LDADD = $(top_srcdir)/src/.libs/libgimli.la $(CPPUNIT_LIBS)
else
gimliUnitTest_LDADD = 
endif

