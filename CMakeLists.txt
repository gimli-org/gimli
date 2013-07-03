# main cmake configuration file

cmake_minimum_required(VERSION 2.8.6)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Sets the configuration to build (Release, Debug, etc...)")
endif()

project(libgimli)

set(LIBGIMLI_VERSION_MAJOR 0)
set(LIBGIMLI_VERSION_MINOR 9)
set(LIBGIMLI_VERSION_PATCH 1)
set(LIBGIMLI_VERSION ${LIBGIMLI_VERSION_MAJOR}.${LIBGIMLI_VERSION_MINOR}.${LIBGIMLI_VERSION_PATCH}) 

set(PACKAGE_NAME \"${PROJECT_NAME}\")
set(PACKAGE_VERSION \"${LIBGIMLI_VERSION}\")
set(PACKAGE_BUGREPORT \"carsten@resistivity.net\")
set(PACKAGE_AUTHORS \"carsten@resistivity.net thomas@resistivity.net\")

################################################################################
# Check for libs and other packages we might use.
################################################################################

include(CheckIncludeFileCXX)

IF(WIN32)
	set(Boost_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/../../boost/include)
endif(WIN32)

find_package(Boost 1.48.0 COMPONENTS python system thread REQUIRED
			)

if (Boost_INCLUDE_DIR)
    message(STATUS "boost include: ${Boost_INCLUDE_DIR}")
    include_directories(${Boost_INCLUDE_DIR})
	set (CMAKE_REQUIRED_INCLUDES ${Boost_INCLUDE_DIR})
		check_include_file_cxx(boost/interprocess/managed_shared_memory.hpp HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP
						   HINT ${Boost_INCLUDE_DIR})
		check_include_file_cxx(boost/bind.hpp BOOST_BIND_FOUND
						   HINT ${Boost_INCLUDE_DIR})
		set(BOOST_BIND_FOUND ON)
	if (HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP)
		message( "HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP ${HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP} ${Boost_INCLUDE_DIR}")
		option (USE_IPC "Use support for interprocess communication" ON) 
	else()
		message( "HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP 0 ${HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP} ${Boost_INCLUDE_DIR}")
		set(USE_IPC  OFF) 
		set(HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP 0 )
	endif()
endif (Boost_INCLUDE_DIR)

find_package(Triangle)
find_package(BLAS)
find_package(LAPACK)
find_package(CHOLMOD)

IF(WIN32)
	set(CppUnit_DIR C:/MinGW/msys/1.0/local/include)
endif(WIN32)
	
find_package(CppUnit)
if (CPPUNIT_FOUND)
    option (LIBGIMLI_BUILD_TESTS "Build unittests" ON) 
endif (CPPUNIT_FOUND)


find_package(Doxygen)

if (DOXYGEN_FOUND)
    add_custom_target(html
        "${DOXYGEN_EXECUTABLE}"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen"
        COMMENT "Generating documentation" VERBATIM)
endif()

IF(WIN32)
	set(READPROC_FOUND FALSE)
ELSE(WIN32)
    find_package(readproc)
ENDIF(WIN32)


configure_file (
  "${PROJECT_SOURCE_DIR}/config.cmake.h.in"
  "${PROJECT_BINARY_DIR}/config.cmake.h"
  )
add_definitions(-DHAVE_CONFIG_CMAKE_H)
# 
# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config1.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)


set(LIBRARY_INSTALL_DIR lib)
set(INCLUDE_INSTALL_DIR include/gimli/)

################################################################################
# Define the configurable options
################################################################################
include(CMakeDependentOption)
cmake_dependent_option( GIMLI_BUILD_SHARED_LIBS_WITH_STATIC_DEPENDENCIES 
"For WIN32 and APPLE where the dependency pack is used, prefer the static dependency libs over the shared/dynamic ones.  NOTE: On Windows you also need to be mindful of which C/C++ runtime setting has been used to compile the various
components - they must all match or you will get crashes, heap corruption and/or
other issues." FALSE "WIN32 OR APPLE" FALSE)

################################################################################
# Add main header locations (for everything we build)
################################################################################
include_directories( ${libgimli_SOURCE_DIR}/src)
include_directories("${PROJECT_BINARY_DIR}")

################################################################################
# descend into subdirs
################################################################################

if (LIBGIMLI_BUILD_TESTS)
     enable_testing()
endif(LIBGIMLI_BUILD_TESTS)

add_subdirectory(src)
add_subdirectory(tests EXCLUDE_FROM_ALL)
add_subdirectory(apps)

################################################################################
# CPack (mostly for source tarballs)
################################################################################
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "GIMLi")
set(CPACK_PACKAGE_VENDOR ${PACKAGE_AUTHORS})
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/doc/gimli.rst")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "${GIMLI_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${GIMLI_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${GIMLI_VERSION_PATCH}")

set(CPACK_SOURCE_GENERATOR TBZ2 ZIP)

set(CPACK_SOURCE_PACKAGE_FILE_NAME "gimlib-${GIMLI_VERSION}" CACHE INTERNAL "tarball basename")

set(CPACK_SOURCE_IGNORE_FILES
# repository history should not be in source tarballs
"\\\\.svn.*"
# most likely not needed in the tarball
"promo/"
"cppcheck-output"
"perform-cppcheck"

# the rest is modelled after .hgignore
"build/"
"~$"

"Thumbs.db"
"\\\\.directory"

"\\\\.kdev4"
"\\\\.settings"
"\\\\.project"
"\\\\.cproject"
"\\\\.pydevproject"
)

include(CPack)