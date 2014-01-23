# main cmake configuration file

cmake_minimum_required(VERSION 2.8.7)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
include(GIMLImacros)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
add_custom_target(pygimli)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Sets the configuration to build (Release, Debug, etc...)")
endif()

project(libgimli)

if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANGXX)
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -pipe -g -ansi -pedantic -Wall -Wno-long-long -Wno-unused-result -Wno-unused-value -Wno-strict-aliasing -Wno-unused-local-typedefs")
    
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-Wl,--no-undefined -static-libgcc -Wl,-O2 -Wl,--as-needed -Wl,--sort-common -s")
endif ()
 
if(MSVC)
	# MSVC complain a lot of possible unsecure std function
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif(MSVC)
 
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

if (BOOST_ROOT)
    get_filename_component(BOOST_ROOT "${BOOST_ROOT}" ABSOLUTE)
	set(Boost_INCLUDE_DIR ${BOOST_ROOT}/include)
else()
	IF(WIN32)
		set(Boost_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/../../boost/include)
	endif(WIN32)
endif()

if (MSVC)
	find_package(Boost 1.46.0 COMPONENTS system thread)
	
	if (NOT Boost_THREAD_FOUND)
		set(Boost_THREAD_FOUND ON)
	endif (NOT Boost_THREAD_FOUND)
	link_directories(${BOOST_LIBRARYDIR})	
	
else(MSVC)
	find_package(Boost 1.46.0 COMPONENTS system thread REQUIRED)
	link_directories(${BOOST_LIBRARYDIR})
endif(MSVC)

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
		message( STATUS "HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP ${HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP} ${Boost_INCLUDE_DIR}")
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

find_package(CppUnit)
if (CPPUNIT_FOUND)
    option (LIBGIMLI_BUILD_TESTS "Build unittests" ON) 
endif (CPPUNIT_FOUND)

find_package(Doxygen)

IF(WIN32)
	set(READPROC_FOUND FALSE)
ELSE(WIN32)
    find_package(readproc)
ENDIF(WIN32)

################################################################################
# Check for python stuff
################################################################################

set(Python_ADDITIONAL_VERSIONS 3.3 3.2 2.7)
find_package(PythonInterp)
if (WIN32)
	get_filename_component(PYTHONPATH "${PYTHON_EXECUTABLE}" PATH )
	set (PYTHON_LIBRARY ${PYTHONPATH}/libs/python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.lib)
	set (PYTHON_INCLUDE_DIR ${PYTHONPATH}/include)
else()
	set(Python_ADDITIONAL_VERSIONS ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
	set(PythonLibs_FIND_VERSION ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
endif()
find_package(PythonLibs)
find_package(Boost 1.46.0 COMPONENTS python )

find_package(gccxml)
find_package(pygccxml)

if(PYTHON_EXECUTABLE AND PYTHONLIBS_FOUND AND GCCXML_EXECUTABLE 
	AND PYGCCXML_FOUND AND PYPLUSPLUS_FOUND
	AND Boost_SYSTEM_LIBRARIES)
    set (PYGIMLI 1)
    message (STATUS "pygimli can be build.")
endif ()

configure_file ("${PROJECT_SOURCE_DIR}/config.cmake.h.in" "${PROJECT_BINARY_DIR}/config.cmake.h" )
configure_file( doc/doxygen/doxyfile.in doc/doxygen/doxyfile )

add_definitions(-DHAVE_CONFIG_CMAKE_H)

if (DOXYGEN_FOUND)
    add_custom_target(html
        "${DOXYGEN_EXECUTABLE}"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen"
        COMMENT "Generating documentation" VERBATIM)
endif()

# 
# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config1.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)

################################################################################
# Define the configurable options
################################################################################
include(CMakeDependentOption)
cmake_dependent_option( GIMLI_BUILD_SHARED_LIBS_WITH_STATIC_DEPENDENCIES 
"For WIN32 and APPLE where the dependency pack is used, prefer the static dependency libs over the shared/dynamic ones.  NOTE: On Windows you also need to be mindful of which C/C++ runtime setting has been used to compile the various
components - they must all match or you will get crashes, heap corruption and/or
other issues." FALSE "WIN32 OR APPLE" FALSE)

# cmake_dependent_option(GIMLI_BUILD_PYTHON_MODULES "Specifies whether to build the Python extension module(s)" "PYTHONINTERP_FOUND;PYTHONLIBS_FOUND;Boost_PYTHON_FOUND")

################################################################################
# 
################################################################################
set (LIBRARY_INSTALL_DIR lib)
set (INCLUDE_INSTALL_DIR include/gimli/)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY_MINSIZEREL ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
if (WIN32)
    set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
else()
    set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
endif()
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY_MINSIZEREL ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# set up install sub-directories
if (CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS "${CMAKE_INSTALL_PREFIX}/lib64")
    set( GIMLI_LIB_INSTALL_DIR lib64 )
elseif (CMAKE_SIZEOF_VOID_P EQUAL 4 AND EXISTS "${CMAKE_INSTALL_PREFIX}/lib32")
    set( GIMLI_LIB_INSTALL_DIR lib32 )
else()
    set( GIMLI_LIB_INSTALL_DIR lib )
endif()

set( GIMLI_VER_INSTALL_SUBDIR "/${CMAKE_PROJECT_NAME}-${GIMLI_VERSION_MAJOR}" )
set( GIMLI_DATA_INSTALL_DIR "share${GIMLI_VER_INSTALL_SUBDIR}" )
set( GIMLI_DOC_INSTALL_DIR "share/doc${GIMLI_VER_INSTALL_SUBDIR}" )

set( GIMLI_PKGCONFIG_INSTALL_DIR "${GIMLI_LIB_INSTALL_DIR}/pkgconfig" )
set( GIMLI_MODULE_INSTALL_DIR "${GIMLI_LIB_INSTALL_DIR}/${CMAKE_PROJECT_NAME}-${GIMLI_VERSION_MAJOR}.${GIMLI_VERSION_MINOR}" )
set( GIMLI_SAMPLE_INSTALL_DIR "${GIMLI_MODULE_INSTALL_DIR}" ) # TODO: put into /samples subdir!
set( GIMLI_INCLUDE_INSTALL_DIR "include${GIMLI_VER_INSTALL_SUBDIR}" )
set( GIMLI_PYTHON_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/python/pygimli)

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

if (PYGIMLI)
    add_subdirectory(python EXCLUDE_FROM_ALL)
endif (PYGIMLI)

#what is this good for?
export(TARGETS gimli
       FILE "${PROJECT_BINARY_DIR}/Gimli.cmake")

################################################################################
# CPack (mostly for source tarballs)
################################################################################
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "GIMLi")
set(CPACK_PACKAGE_VENDOR ${PACKAGE_AUTHORS})
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.rst")
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
"cppcheck-output"
"perform-cppcheck"

# the rest is modeled after .hgignore
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
