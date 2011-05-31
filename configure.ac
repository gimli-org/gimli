# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.61)

define([SVNVERSION], esyscmd([sh -c "svnversion -n"]))
AC_DEFINE(SVN_REVISION, "SVNVERSION", [SVN Revision])

AC_INIT( [gimli], [0.7.0-SVNVERSION], [carsten@resistivity.net] )
AM_INIT_AUTOMAKE( gimli, 0.7.0 )

# init magic libtool
# -all-static
LT_PREREQ([2.2])
LT_INIT([dlopen,win32-dll])

AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADER([config.h])

AC_PREFIX_DEFAULT(`pwd`)

# Checks for programs.
AC_PROG_CXX

AC_SUBST([CXXFLAGS], ["-O2 -g -Wall -fPIC -ansi -pedantic"])

#CXXFLAGS="$CXXFLAGS -O2 -g -Wall -fPIC -ansi -pedantic #-D svn-\`svnversion -n\` "
#CXXFLAGS="$CXXFLAGS -O2 -g -Wall -fPIC -ansi -pedantic -D ${PACKAGE_VERSION}"
#AC_SUBST([CPPFLAGS], ["${CPPFLAGS} -Wall"])

# Checks for system.
#AC_MINGW32

# Checks for libraries.
AC_CHECK_LIB(m, sqrt )

# Checks for header files.
AC_CHECK_HEADERS([sys/timeb.h])
AC_FUNC_ALLOCA
AC_HEADER_STDC
AC_HEADER_SYS_WAIT

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([strerror])

#AC_CHECK_LIB( [boost_thread], [main])

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
#AC_CHECK_LIB( pthread, [ _pthread_atfork], [], [], [ -lgfortran ]  )
#AC_CHECK_LIB( gfortran, [_gfortran_free] )


AC_MSG_RESULT()

dnl ############################################################################################
dnl # Checks for local external libraries and haeders
dnl ############################################################################################

GIMLI_PATH=`pwd`
GIMLI_EXTERNAL_LIBS_PATH=$GIMLI_PATH/external/lib
LDFLAGS="$LDFLAGS -L$GIMLI_EXTERNAL_LIBS_PATH"

AC_CHECK_LIB( [triangle], [triangulate])

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
                ] )
#                tests/Makefile
#                tests/unittest/Makefile
#                apps/Makefile
#                apps/bert/Makefile
#                apps/em1d/Makefile
#                apps/fit/Makefile
#                apps/jointinv/Makefile
#                apps/misc/Makefile
#                apps/traveltime/Makefile

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
Libraries:              ${LDFLAGS} ${LIBS}  )

AC_MSG_RESULT()
AC_MSG_RESULT([Configuration successful! Now type 'make' to build the library, tests, applications and examples])
AC_MSG_RESULT()
