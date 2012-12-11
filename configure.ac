# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.61])

#define( [SUMMARY], "GIMLi summary" )
#AC_DEFINE( SUMMARY )

define([SVNVERSION], esyscmd([sh -c "svnversion -n"]))
AC_DEFINE(SVN_REVISION, "SVNVERSION", [SVN Revision])

#http://sourceforge.net/projects/libgimli/

#AC_INIT(package, version, bug-report-address)
#     information on the package
#     checks for programs
#     checks for libraries
#     checks for header files
#     checks for types
#     checks for structures
#     checks for compiler characteristics
#     checks for library functions
#     checks for system services
#     AC_CONFIG_FILES([file...])
#     AC_OUTPUT

AC_INIT( [gimli], 0.9.0.SVNVERSION, [carsten@resistivity.net] )

GIMLI_MAJOR_VERSION=0
GIMLI_MINOR_VERSION=9
GIMLI_MICRO_VERSION=0
GIMLI_INTERFACE_AGE=0
GIMLI_BINARY_AGE=0

AC_SUBST(GIMLI_MAJOR_VERSION)
AC_SUBST(GIMLI_MINOR_VERSION)
AC_SUBST(GIMLI_MICRO_VERSION)
AC_SUBST(GIMLI_INTERFACE_AGE)
AC_SUBST(GIMLI_BINARY_AGE)

GIMLI_VERSION=$GIMLI_MAJOR_VERSION.$GIMLI_MINOR_VERSION.$GIMLI_MICRO_VERSION.$SVNVERSION
#AC_CONFIG_SRCDIR([src/foo.c])

AM_INIT_AUTOMAKE
AM_SILENT_RULES([yes])

################################################################################
# Init magic libtool
################################################################################

LT_RELEASE=$GIMLI_MAJOR_VERSION.$GIMLI_MINOR_VERSION
LT_CURRENT=`expr $GIMLI_MICRO_VERSION - $GIMLI_INTERFACE_AGE`
LT_REVISION=$GIMLI_INTERFACE_AGE
LT_AGE=`expr $GIMLI_BINARY_AGE - $GIMLI_INTERFACE_AGE`
AC_SUBST(LT_RELEASE)
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)

LT_PREREQ([2.2])
LT_INIT([dlopen,win32-dll])

AC_SUBST([CR_ARCH],$host_cpu)
AC_SUBST([CR_CPU],$host_cpu)

################################################################################
# Check for rpmbuild
################################################################################

#AC_PATH_PROGS([RPMBUILD],[rpmbuild rpmb rpm],[none],[$PATH:/usr/lib/rpm])
AC_PATH_PROGS([RPMBUILD],[rpmbuild rpm],[none],[$PATH:/usr/lib/rpm])
if $RPMBUILD -bs 2>&1 | grep 'Usage' >/dev/null 2>/dev/null; then
    :
else
    RPMBUILD=none
fi
AM_CONDITIONAL(HAVE_RPMBUILD,[test $RPMBUILD != none])

################################################################################


#BB_ENABLE_DOXYGEN

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([Makefile.am]) 
#AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADER([config.h])

AC_PREFIX_DEFAULT(`pwd`)

################################################################################
# Checks for programs.
################################################################################

AC_PROG_CXX
AC_LANG([C++])

# check why this doesn't work
#AS_COMPILER_FLAGS( GLOBAL_CFLAGS, "-Wall -ansi -pedantic" )

AC_SUBST([CXXFLAGS], ["-O2 -g -Wall -ansi -pedantic -pipe"])
#AC_SUBST([CXXFLAGS], ["-O2 -g -Wall -ansi -pedantic -Wl,-noinhibit-exec"])

#AC_SUBST([CXXFLAGS], ["-O2 -g -Wall -fPIC -ansi -pedantic"])
#CXXFLAGS="$CXXFLAGS -O2 -g -Wall -fPIC -ansi -pedantic #-D svn-\`svnversion -n\` "
#CXXFLAGS="$CXXFLAGS -O2 -g -Wall -fPIC -ansi -pedantic -D ${PACKAGE_VERSION}"
#AC_SUBST([CPPFLAGS], ["${CPPFLAGS} -Wall"])


#
# Checks for header files.
#
AC_FUNC_ALLOCA
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_HEADER_STDBOOL

AC_CHECK_HEADER([sys/timeb.h])
AC_CHECK_HEADER([proc/readproc.h])
if test x"$ac_cv_header_proc_readproc_h" = "xyes"; then
    AC_DEFINE( [HAVE_PROC_READPROC], [1], [Define and we can use memwatch])  
fi


#
# Checks for libraries.
#
AC_CHECK_LIB( m, sqrt )
AC_CHECK_LIB( proc, look_up_our_self )

#
# CHECK for boost/bind.hpp
#
BOOST_REQUIRE([1.41],[AC_MSG_ERROR([Sorry, boost is requried.
                    You can add your boost installation path (absolute) to the configure call, e.g.,
                    BOOST_ROOT=/c/home/boost/ ./autogen.sh
					or
					./configure --with-boost=/c/home/boost])])
BOOST_FIND_HEADER([boost/bind.hpp])
BOOST_FIND_HEADER([boost/thread.hpp])
BOOST_FIND_HEADER([boost/interprocess/managed_shared_memory.hpp],
                  [AC_MSG_WARN("*******************************************"
                               "does not found or non-working boost/interprocess/managed_shared_memory.hpp"
                               "*******************************************")] )
AC_SUBST([CPPFLAGS], ["${CPPFLAGS} ${BOOST_CPPFLAGS}"])
AC_CHECK_LIB( [boost_thread], [main])


#
# Checks for typedefs, structures, and compiler characteristics.
#

AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

#
# Checks for library functions.
#
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([strerror])

#	
# Checks for packages.
#
#PKG_CHECK_MODULES( CHECK, [ check >= 0.8.0 ] ) 

AM_PATH_CPPUNIT(,[have_unittest="yes"],
  AC_MSG_WARN([cppunit not found; cannot run unit tests!])
  [have_unittest="no"])

AM_CONDITIONAL(HAVE_UNITTEST, test x"$have_unittest" = "xyes")
if test x"$have_unittest" = "xyes"; then
	AC_DEFINE( [HAVE_UNITTEST], [1], [Define and we can use unittest])  
fi

dnl ############################################################################################
dnl # define arguments for configure command line
dnl ############################################################################################

#AC_ARG_VAR( VAR_NAME , [ Description ] )
#AC_ARG_WITH( tool, AC_HELP_STRING([--with-tools], [check for development tools not necessary for a simple compile (default is NO)]),
#                ac_cv_use_devtools=yes, ac_cv_use_devtools=no)

dnl ############################################################################################
dnl # Checks for global external numeric libraries on the current system defaults (gentoo
dnl ############################################################################################

#AC_SEARCH_LIBS( function, search-libs, [action-if-found], [action-if-not-found], [other-libraries])

AC_CHECK_HEADERS( [ldl.h] )
AC_CHECK_LIB( [ldl],     [main])

AC_CHECK_HEADERS( [amd.h] )
AC_CHECK_LIB( [amd],     [main])

AC_CHECK_HEADERS( [cholmod.h] )
AC_CHECK_LIB( [cholmod], [main])

# for mkl only
#AC_CHECK_LIB( gfortran, [ _gfortran_free ], [], [], [-pthread] )
#AC_CHECK_LIB( pthread, [ _pthread_atfork ]  )
#AC_SUBST([CXXFLAGS], [$CXXFLAGS "-pthread"])

# for mkl only
#AC_CHECK_LIB( pthread, [_pthread_atfork], [], [], [ -lgfortran ]  )
#AC_CHECK_LIB( gfortran, [_gfortran_free] )

AC_MSG_RESULT()

dnl ############################################################################################
dnl # Checks for local external libraries and haeders
dnl ############################################################################################

GIMLI_PATH=`pwd`
GIMLI_EXTERNAL_LIBS_PATH=$GIMLI_PATH/external/lib
#LDFLAGS_OLD="$LDFLAGS"
LDFLAGS="$LDFLAGS -L$GIMLI_EXTERNAL_LIBS_PATH"

AC_CHECK_LIB( [triangle], [triangulate], [], [], [-lm] )
#$LDFLAGS="$LDFLAGS_OLD"
#AC_DEFINE( [HAVE_LIBTRIANGLE], [1], [Define and we can use unittest])  
#LDFLAGS="$LDFLAGS -L$GIMLI_EXTERNAL_LIBS_PATH -Wl,-ltriangle,-lcrtdll,-lmingwex,-lkernel32,-lgcc"

if test "$ac_cv_lib_ldl_main" == "no"; then
    AC_CHECK_LIB( [ldl],          [ldl_symbolic])
fi

if test "$ac_cv_lib_amd_main" == "no"; then
    AC_CHECK_LIB( [amd], [amd_order])
fi

if test "$ac_cv_lib_cholmod_main" == "no"; then
    AC_CHECK_LIB( [cholmod], [cholmod_start])
fi

if test "$ac_cv_lib_cholmod_main" == "yes" || test "$ac_cv_lib_cholmod_cholmod_start" == "yes"; then
    AC_SEARCH_LIBS( colamd, colamd )
    AC_SEARCH_LIBS( dgemm_, blas )
    AC_SEARCH_LIBS( dpotrf_, lapack )
fi

#AC_CHECK_LIB( [blas],   [main])
#AC_CHECK_LIB( [lapack], [main])
#AC_CHECK_LIB( [colamd], [colamd_printf] )

AC_CONFIG_FILES([Makefile
                src/Makefile
                tests/Makefile
                tests/unittest/Makefile
                apps/Makefile
                apps/fit/Makefile
                apps/em1d/Makefile
                apps/traveltime/Makefile
                rpm/Makefile
                rpm/specfile
                ] )

#                apps/bert/Makefile
#                apps/jointinv/Makefile
#                apps/misc/Makefile


AC_OUTPUT

AC_MSG_RESULT(
-------------------------------------------------

${PACKAGE_NAME} Version ${PACKAGE_VERSION} )

AC_MSG_RESULT( );

#if test "$ac_cv_lib_boost_thread_main" == "yes"; then
#	AC_MSG_RESULT(  boost::thread           [FOUND])
#else
#	AC_MSG_RESULT(  boost::thread           [Not FOUND])
#	AC_MSG_RESULT(  		(recommended for use with multicore cpu ert-sensitivity calculation) )
#fi

if test "$ac_cv_lib_ldl_main" == "yes"; then
    AC_MSG_RESULT(  ldl (system)            [FOUND])
elif test "$ac_cv_lib_ldl_ldl_symbolic" == "yes"; then
    AC_MSG_RESULT(  ldl (local)             [FOUND])
else
    AC_MSG_RESULT(  ldl                     [Not FOUND])
    AC_MSG_RESULT(           (basic linear solver))
fi

if test "$ac_cv_lib_amd_main" == "yes"; then
    AC_MSG_RESULT(  amd (system)            [FOUND])
elif test "$ac_cv_lib_amd_amd_order" == "yes"; then
	AC_MSG_RESULT(  amd (local)             [FOUND])
else
	AC_MSG_RESULT(  amd                     [Not FOUND])
	AC_MSG_RESULT(			    (needed for ldl))
fi

if test "$ac_cv_lib_cholmod_main" == "yes"; then
    AC_MSG_RESULT(  cholmod (system)        [FOUND])
elif test "$ac_cv_lib_cholmod_cholmod_start" == "yes"; then
    AC_MSG_RESULT(  cholmod (local)         [FOUND])
else
    AC_MSG_RESULT(  cholmod                 [Not FOUND])
    AC_MSG_RESULT(           (recommended for faster linear solving) )
fi

AC_MSG_RESULT( );

if test "$ac_cv_lib_triangle_triangulate" == "yes"; then
	AC_MSG_RESULT(  triangle		[FOUND])
else
	AC_MSG_RESULT(  triangle		[Not FOUND])
	AC_MSG_RESULT(			 (recommended for inbuild meshing) )
fi

AC_MSG_RESULT()

AC_MSG_RESULT( 

GIMLi configuration complete for ${host}

Source directory:       ${srcdir} 
Installation directory: ${prefix}
C++ compiler:           ${CXX}-${CXXVERSION}
OSTYPE:                 ${OSTYPE}
CXXFLAGS:               ${CXXFLAGS}
CPPFLAGS:               ${CPPFLAGS}
Libraries:              ${LDFLAGS} ${LIBS}  )

AC_MSG_RESULT()

if test x"$have_unittest" = "xyes"; then
	AC_MSG_RESULT( cppunit 		[unittest INSTALLED])
else
	AC_MSG_RESULT( cppunit 		[unittest Not INSTALLED])
fi

AC_MSG_RESULT()
AC_MSG_RESULT([Configuration successful! Now type 'make' to build the library, tests, applications and examples])
AC_MSG_RESULT()
