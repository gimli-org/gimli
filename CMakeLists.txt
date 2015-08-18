# main cmake configuration file

cmake_minimum_required(VERSION 2.8.8)

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR} )
    message( FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt." )
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
include(GIMLImacros)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
add_custom_target(pygimli)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Sets the configuration to build (Release, Debug, etc...)")
endif()

project(libgimli)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -Os -pipe -ansi -Wall -Wno-long-long -Wno-unused-result")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-unused-value -Wno-strict-aliasing -Wno-unused-local-typedefs")
#     set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -rdynamic")
 
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-Os -Wl,--no-undefined -static-libgcc -Wl,-O2 -Wl,--as-needed -Wl,--sort-common")
#     set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS} -lbfd")

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") 
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-overloaded-virtual")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Qunused-arguments")
    endif()

    if(WIN32)
    endif(WIN32)

    #-stdlib=libc++
elseif(MSVC)
	# MSVC complain a lot of possible unsecure std function
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()
 
set(LIBGIMLI_VERSION_MAJOR 0)
set(LIBGIMLI_VERSION_MINOR 9)
set(LIBGIMLI_VERSION_PATCH 7)
set(LIBGIMLI_VERSION ${LIBGIMLI_VERSION_MAJOR}.${LIBGIMLI_VERSION_MINOR}.${LIBGIMLI_VERSION_PATCH}) 

set(PACKAGE_NAME \"${PROJECT_NAME}\")
set(PACKAGE_VERSION \"${LIBGIMLI_VERSION}\")
set(PACKAGE_BUGREPORT \"carsten@resistivity.net\")
set(PACKAGE_AUTHORS \"carsten@resistivity.net thomas@resistivity.net\")

################################################################################
# Check for libs and other packages we might use.
################################################################################

if(NOT THIRDPARTY_DIR)
    if (NOT ADDRESSMODEL)
        if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
            message(STATUS "Target is 64 bits")
            set (ADDRESSMODEL "64")
        else()
            message(STATUS "Target is 32 bits")
            set (ADDRESSMODEL "32")
        endif()
    endif()

    set(TARGETNAME "-${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}-${ADDRESSMODEL}")
    set(ENV{TARGETNAME} ${TARGETNAME})
    set(THIRDPARTY_DIR ${PROJECT_SOURCE_DIR}/../thirdParty/)
    
    get_filename_component(EXTERNAL_DIR "${THIRDPARTY_DIR}/dist${TARGETNAME}" ABSOLUTE)
    get_filename_component(EXTERNAL_BUILD_DIR "${THIRDPARTY_DIR}/build${TARGETNAME}" ABSOLUTE)
    message(STATUS "ThirdParty set to: ${THIRDPARTY_DIR}")
    message(STATUS "External set to: ${EXTERNAL_DIR}")
endif()

if(NOT EXTERNAL_DIR)
    set(EXTERNAL_DIR ${PROJECT_SOURCE_DIR}/external/)
    message(STATUS "External set to: ${EXTERNAL_DIR}")
endif()

include(CheckIncludeFileCXX)

if (NOT Boost_INCLUDE_DIR)
    if (BOOST_ROOT)
        get_filename_component(BOOST_ROOT "${BOOST_ROOT}" ABSOLUTE)
        set(Boost_INCLUDE_DIR ${BOOST_ROOT}/include)
    else()
        find_or_build_package_check(Boost boost Boost_INCLUDE_DIR)
        # we should think about cmake superbuild with ExternalProject_Add as the replacement for find_or_build_package (next generation)
        if (NOT Boost_INCLUDE_DIR)
            FILE(READ "${EXTERNAL_DIR}/.boost.dist" BOOST_DIST_NAME)
            STRING(REGEX REPLACE "\n" "" BOOST_DIST_NAME "${BOOST_DIST_NAME}")

            get_filename_component(BOOST_ROOT "${EXTERNAL_DIR}/${BOOST_DIST_NAME}" ABSOLUTE)

            message(STATUS "BOOST_ROOT from ${EXTERNAL_DIR}/.boost.dist: ${BOOST_ROOT}")

            set(BOOST_ROOT ${BOOST_ROOT}  CACHE FILEPATH "Boost root distribution")
	
            set(Boost_INCLUDE_DIR ${BOOST_ROOT}/include)
            get_filename_component(Boost_INCLUDE_DIR "${Boost_INCLUDE_DIR}" ABSOLUTE)
            message(STATUS "Boost_INCLUDE_DIR, ${Boost_INCLUDE_DIR}")	
        endif()
        #exit()
    endif(BOOST_ROOT)
    message(STATUS "BOOST_ROOT, ${BOOST_ROOT}")
endif()

message(STATUS "Boost_INCLUDE_DIR, ${Boost_INCLUDE_DIR}")

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

mark_as_advanced(Boost_INCLUDE_DIR BOOST_ROOT)

# Find packages that CHOLMOD depends on
set(CMAKE_LIBRARY_PATH ${EXTERNAL_DIR}/lib $ENV{EXTERNAL_DIR}/lib ${CMAKE_LIBRARY_PATH})

if (J) # dummy to avoid error msg
endif()

find_or_build_package(Triangle triangle)
find_or_build_package(BLAS lapack)
find_or_build_package(LAPACK lapack)
find_or_build_package(CHOLMOD suitesparse)
find_or_build_package(UMFPACK umfpack)

find_package(CppUnit)
if (CPPUNIT_FOUND)
    option (LIBGIMLI_BUILD_TESTS "Build unittests" ON) 
endif (CPPUNIT_FOUND)

IF(WIN32)
	set(READPROC_FOUND FALSE)
ELSE(WIN32)
    find_package(readproc)
ENDIF(WIN32)

################################################################################
# Check for python stuff
################################################################################

if (PYVERSION)
    set(Python_ADDITIONAL_VERSIONS ${PYVERSION})
    #set(Python_ADDITIONAL_VERSIONS 3.4 3.3 3.2 2.7)
endif()

find_package(PythonInterp) # find version belongs to the active python command

if (WIN32)
	get_filename_component(PYTHONPATH "${PYTHON_EXECUTABLE}" PATH )
	set (PYTHON_LIBRARY ${PYTHONPATH}/libs/libpython${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.a)
	#set (PYTHON_LIBRARY ${PYTHONPATH}/libs/python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.lib)
	set (PYTHON_INCLUDE_DIR ${PYTHONPATH}/include)
else()
	set(Python_ADDITIONAL_VERSIONS ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
	set(PythonLibs_FIND_VERSION ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
endif()

find_package(PythonLibs)

find_python_module(numpy)

if (Boost_PYTHON_LIBRARY)
    set (Boost_PYTHON_LIBRARIES ${Boost_PYTHON_LIBRARY})
else()
    find_package(Boost 1.46.0 COMPONENTS python)
endif()

if (NOT CASTER)
    set(CASTER "gccxml")

    if (CMAKE_COMPILER_IS_GNUCC)
        string(REGEX MATCHALL "[0-9]+" GCC_VERSION_COMPONENTS ${CMAKE_CXX_COMPILER_VERSION})
        list(GET GCC_VERSION_COMPONENTS 0 GCC_MAJOR)
        list(GET GCC_VERSION_COMPONENTS 1 GCC_MINOR)

        if (${GCC_MAJOR} EQUAL 5)
            message(STATUS "Found gcc version above 5: choosing castxml")
            set(CASTER "castxml")
        endif()
            
    endif()
endif()

if (CASTER STREQUAL "gccxml")
    find_or_build_package(gccxml gccxml)
    set(CASTER_EXECUTABLE ${GCCXML_EXECUTABLE})
else()
    find_program(CLANG clang REQUIRED)
    find_program(CLANGXX clang++ REQUIRED)
    find_package(LLVM REQUIRED)

    find_or_build_package(castxml castxml)
    set(CASTER_EXECUTABLE ${CASTXML_EXECUTABLE})
endif()

find_or_build_package(pygccxml pygccxml)
find_or_build_package(pyplusplus pygccxml)

if (PYTHON_EXECUTABLE AND PYTHONLIBS_FOUND 
    AND (GCCXML_FOUND OR CASTXML_FOUND)
    AND PYGCCXML_FOUND AND PYPLUSPLUS_FOUND
    AND Boost_PYTHON_LIBRARIES AND numpy_FOUND)
    set (PYGIMLI 1)
    message (STATUS "pygimli can be build.")
else()
    message(STATUS "PYTHON_EXECUTABLE: ${PYTHON_EXECUTABLE}" )
    message(STATUS "PYTHONLIBS_FOUND: ${PYTHONLIBS_FOUND}" )
    message(STATUS "GCCXML_EXECUTABLE: ${GCCXML_EXECUTABLE}" )
    message(STATUS "CASTXML_EXECUTABLE: ${CASTXML_EXECUTABLE}" )
    message(STATUS "PYPLUSPLUS_FOUND: ${PYPLUSPLUS_FOUND}" )
    message(STATUS "Boost_PYTHON_LIBRARIES: ${Boost_PYTHON_LIBRARIES}" )
    message(STATUS "numpy_FOUND: ${numpy_FOUND}" )

    message (STATUS "pygimli cannot be build due to some missing packages. 
Check above for python, gccxml|castxml, pygccxml, pyplusplus, boost_python, numpy")
    
endif ()

message(STATUS "**********************************************************************")
message(STATUS "************************* Dependencies found *************************")
message(STATUS "**********************************************************************")
message(STATUS "Triangle_LIBRARIES: ${Triangle_LIBRARIES}")
message(STATUS "UMFPACK_LIBRARIES: ${UMFPACK_LIBRARIES}")
message(STATUS "CHOLMOD_LIBRARIES: ${CHOLMOD_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE: ${PYTHON_EXECUTABLE}" )
message(STATUS "PYTHON_LIBRARY: ${PYTHON_LIBRARY}" )
message(STATUS "Boost_PYTHON_LIBRARIES: ${Boost_PYTHON_LIBRARIES}" )
message(STATUS "XML Caster: ${CASTER_EXECUTABLE}")
message(STATUS "PYGCCXML: ${PYGCCXML_PATH}" )
message(STATUS "PYPLUSPLUS: ${PYPLUSPLUS_PATH}" )
message(STATUS "PY_NUMPY: ${PY_NUMPY}" )
message(STATUS "**********************************************************************")

configure_file("${PROJECT_SOURCE_DIR}/config.cmake.h.in" "${PROJECT_BINARY_DIR}/config.cmake.h" )
add_definitions(-DHAVE_CONFIG_CMAKE_H)

find_package(Doxygen)

find_package(Sphinx 1.2.0)
if (SPHINX_FOUND)
    add_subdirectory(doc EXCLUDE_FROM_ALL)
    file(GLOB_RECURSE DOCFILES RELATIVE "${CMAKE_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}" "*.rst")
    file(GLOB_RECURSE PNGFILES RELATIVE "${CMAKE_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}" "*.png")
    list(APPEND DOCFILES ${PNGFILES})
    file(GLOB_RECURSE BIBFILES RELATIVE "${CMAKE_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}" "*.bib")
    list(APPEND DOCFILES ${BIBFILES})
    
    set(RSTFILES ${DOCFILES} CACHE INTERNAL "RST source file that we need to copy")

    #foreach(file ${DOCFILES})
        #message(${file})
    #endforeach()

endif()


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
#set( GIMLI_PYTHON_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/python/pygimli)

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

if (PYTHON_EXECUTABLE AND PYTHONLIBS_FOUND 
    AND (GCCXML_FOUND OR CASTXML_FOUND)
    AND PYGCCXML_FOUND AND PYPLUSPLUS_FOUND
    AND Boost_PYTHON_LIBRARIES AND numpy_FOUND)
    set (PYGIMLI 1)
    message (STATUS "pygimli can be build.")
else()
    message(STATUS "PYTHON_EXECUTABLE: ${PYTHON_EXECUTABLE}" )
    message(STATUS "PYTHONLIBS_FOUND: ${PYTHONLIBS_FOUND}" )
    message(STATUS "XML Caster: ${CASTER_EXECUTABLE}")
    message(STATUS "PYPLUSPLUS_FOUND: ${PYPLUSPLUS_FOUND}" )
    message(STATUS "Boost_PYTHON_LIBRARIES: ${Boost_PYTHON_LIBRARIES}" )
    message(STATUS "numpy_FOUND: ${numpy_FOUND}" )

    message (STATUS "pygimli cannot be build due to some missing packages. 
Check above for python, gccxml|castxml, pygccxml, pyplusplus, boost_python, numpy")
    
endif ()

if (PYGIMLI)
#     set( PYGIMLI_SOURCE_DIR ${CMAKE_SOURCE_DIR}/python/pygimli )
#     set( PYGIMLI_BUILD_DIR ${CMAKE_BINARY_DIR}/package/pygimli )
# 
#     file(GLOB_RECURSE PYFILES RELATIVE "${PYGIMLI_SOURCE_DIR}" "python/pygimli/*.py")
#     set(PYGIMLI_FILES ${PYFILES} CACHE INTERNAL "pygimli source files that we need to copy")

    add_subdirectory(python EXCLUDE_FROM_ALL)
endif (PYGIMLI)

#CR what is this good for?
# export(TARGETS gimli
#        FILE "${PROJECT_BINARY_DIR}/Gimli.cmake")

