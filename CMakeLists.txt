cmake_minimum_required(VERSION 3.12) 

if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed.
            Please make a new directory (called a build directory)
            and run CMake from there. You may need to remove CMakeCache.txt." )
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/core/cmake")
include(GIMLImacros)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
add_custom_target(pygimli)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Sets the configuration to build (Release, Debug, etc...)")
endif()

project(libgimli)

include(CheckIncludeFileCXX)

message(STATUS "System: " ${CMAKE_SYSTEM_NAME} "; " ${CMAKE_SYSTEM} "; " ${CMAKE_SYSTEM_VERSION})

message(STATUS "UNIX:   " ${UNIX} " TRUE on all UNIX-like OS's, including Apple OS X and CygWin")
message(STATUS "WIN32:  " ${WIN32} " TRUE on Windows. Prior to 2.8.4 this included CygWin")
message(STATUS "APPLE:  " ${APPLE} " TRUE on Apple systems. Note this does not imply the
system is Mac OS X, only that APPLE is")
message(STATUS "MINGW:  " ${MINGW} " TRUE when using the MinGW compiler in Windows")
message(STATUS "MSYS:   " ${MSYS} " TRUE when using the MSYS developer environment in Windows")
message(STATUS "CYGWIN: " ${CYGWIN} " TRUE on Windows when using the CygWin version of cmake")

# Check if conda package is created
if(DEFINED ENV{CONDA_BUILD})
    message(STATUS "Conda package is being created.")
    if (APPLE)
	    set(CMAKE_MACOSX_RPATH 1)
	    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    endif()
    set(CONDA_BUILD TRUE)
    set(Boost_INCLUDE_DIR "${CMAKE_PREFIX_PATH}/include")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-class-memaccess") # new since gcc.8 maybe need to be checked, not on apple
else()
    set(CONDA_BUILD FALSE)
endif()

if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANGXX)

    add_definitions(-std=c++11)
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-Os")

    # run with -DCMAKE_BUILD_TYPE=Debug
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -p -g -O2 -std=c++11 -ansi -pedantic -fno-omit-frame-pointer -Wall")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wno-class-memaccess") # new since gcc.8 maybe need to be checked

    if (WIN32 OR CONDA_BUILD)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-deprecated-declarations -Wno-class-memaccess")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-class-memaccess") # new since gcc.8 maybe need to be checked
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-parentheses") # from boost since gcc.8 maybe need to be checked
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-attributes") # from gcc-8.2 + boost
    endif()

    if (WIN32 OR CONDA_BUILD)
        set(CMAKE_CXX_FLAGS_RELEASE
            "${CMAKE_CXX_FLAGS_RELEASE} -Wno-deprecated-declarations")

        # from boost since gcc.8 maybe need to be checked
        set(CMAKE_CXX_FLAGS_RELEASE
            "${CMAKE_CXX_FLAGS_RELEASE} -Wno-parentheses")
        # from gcc-8.2 + boost
        set(CMAKE_CXX_FLAGS_RELEASE
            "${CMAKE_CXX_FLAGS_RELEASE} -Wno-attributes")
    endif()

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS_RELEASE
            "${CMAKE_CXX_FLAGS_RELEASE} -Wno-overloaded-virtual")
        set(CMAKE_CXX_FLAGS_RELEASE
            "${CMAKE_CXX_FLAGS_RELEASE} -Qunused-arguments")

        if (APPLE)
            ## needet?
            # If it's libc++ and you're on <= 10.8, you need to compile with clang++ -stdlib=libc++. If it's libstdc++ and you're on 10.9 or later, you need to compile with clang++ -stdlib=libstdc++.
            # add_compile_options(-stdlib=libstdc++)
            # add_compile_options(-stdlib=libc++)
        endif()
    else() # if gcc
    endif()

    ### debug settings
    set(CMAKE_CXX_FLAGS_DEBUG
    "${CMAKE_CXX_FLAGS_DEBUG} -p -g -ansi -pedantic -fno-omit-frame-pointer -Wall")

    if (NOT WIN32 AND ASAN)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
    endif()
elseif(MSVC)
	# MSVC complain a lot of possible unsecure std function
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

find_package(Git)
if (GIT_FOUND)
    message(STATUS, ${GIT_FOUND})
    execute_process(
        COMMAND
            ${GIT_EXECUTABLE} describe --tags
        WORKING_DIRECTORY
            ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE
            GIMLI_GITVERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_STRIP_TRAILING_WHITESPACE
        )
    message(STATUS, "setting version from git description: ${GIMLI_GITVERSION}")
    set(LIBGIMLI_VERSION ${GIMLI_GITVERSION})
else()
    set(LIBGIMLI_VERSION "untagged")
endif()

set(PACKAGE_NAME  \"${PROJECT_NAME}\" )
set(PACKAGE_VERSION  \"${LIBGIMLI_VERSION}\" )
set(PACKAGE_BUGREPORT  \"carsten@pygimli.org\")
set(PACKAGE_AUTHORS  \"carsten@pygimli.org thomas@pygimli.org florian@pygimli.org\")

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

find_package(Threads REQUIRED)

################################################################################
# Check for python stuff
################################################################################
if (NOT PYVERSION)
    set(PYVERSION 3)
endif()

if (NOT PYTHON_EXECUTABLE AND NOT ${CMAKE_VERSION} VERSION_LESS "3.14.0")
    message(STATUS "Searching for python:" ${PYVERSION})
    find_package(Python ${PYVERSION} EXACT COMPONENTS Interpreter NumPy Development REQUIRED)

    # message(STATUS "Python_Interpreter_FOUND: ${Python_Interpreter_FOUND} ${Python_EXECUTABLE} version: ${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}.${Python_VERSION_PATCH}")
    # message(STATUS "Python_INCLUDE_DIRS: ${Python_INCLUDE_DIRS}")
    # message(STATUS "Python_LIBRARIES: ${Python_LIBRARIES} ${Python_LIBRARY_DIRS}")
    # message(STATUS "Python_NumPy: ${Python_NumPy_FOUND}, ${Python_NumPy_INCLUDE_DIRS}, ${Python_NumPy_VERSION}")
    set (PYTHON_FOUND TRUE CACHE BOOL "Is python available")
else()
    
    message(STATUS "Manually searching for python (pre cmake 3.14): ${PYVERSION}")
    find_package(PythonInterp ${PYVERSION} REQUIRED)
    find_package(PythonLibs)
    
    set(Python_EXECUTABLE ${PYTHON_EXECUTABLE} CACHE FILEPATH "Python interpreter")
    find_python_module(numpy)
    #set(Python_Libraries ${PYTHON_LIBRARY})
    
    message(STATUS "Python_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIR}")
    message(STATUS "Python_LIBRARIES: ${PYTHON_LIBRARY}")
    message(STATUS "Python_NumPy: ${numpy_FOUND}, ${PY_NUMPY}")
    
    set(Python_VERSION_MAJOR ${PYTHON_VERSION_MAJOR})
    set(Python_VERSION_MINOR ${PYTHON_VERSION_MINOR})
    set(Python_LIBRARIES ${PYTHON_LIBRARY} CACHE FILEPATH "python library")
    set(Python_INCLUDE_DIRS ${PYTHON_INCLUDE_DIR} CACHE FILEPATH "python includes")
    set(Python_NumPy_FOUND ${numpy_FOUND} CACHE BOOL "is numpy installed")
    set(Python_NumPy_INCLUDE_DIRS  ${PY_NUMPY}/core/include CACHE FILEPATH "Numpy includes")

    set (PYTHON_FOUND TRUE CACHE BOOL "Is python available")
endif()

function(find_boost_from_distname_file)
    set (BOOST_DIST_NAME_FILE "${EXTERNAL_DIR}/.boost-py${Python_VERSION_MAJOR}.dist")
    message(STATUS "locking for ${BOOST_DIST_NAME_FILE}")

    if (EXISTS ${BOOST_DIST_NAME_FILE})
        message(STATUS "reading ${BOOST_DIST_NAME_FILE}")
        FILE(READ "${EXTERNAL_DIR}/.boost-py${Python_VERSION_MAJOR}.dist" BOOST_DIST_NAME)
        STRING(REGEX REPLACE "\n" "" BOOST_DIST_NAME "${BOOST_DIST_NAME}")
        get_filename_component(BOOST_ROOT "${EXTERNAL_DIR}/${BOOST_DIST_NAME}" ABSOLUTE)
    else()
        message(STATUS "there is no ${BOOST_DIST_NAME_FILE}")
    endif()

    if (EXISTS ${BOOST_ROOT})
        message(STATUS "BOOST_ROOT from ${EXTERNAL_DIR}/.boost-py${Python_VERSION_MAJOR}.dist: ${BOOST_ROOT}")

        set(BOOST_ROOT ${BOOST_ROOT}  CACHE FILEPATH "Boost root distribution")
        set(Boost_INCLUDE_DIR ${BOOST_ROOT}/include CACHE FILEPATH "Boost inlude path")
        get_filename_component(Boost_INCLUDE_DIR "${Boost_INCLUDE_DIR}" ABSOLUTE)
        message(STATUS "found Boost_INCLUDE_DIR: ${Boost_INCLUDE_DIR}")
    else()
        message(STATUS "BOOST_ROOT does not exist: ${BOOST_ROOT}")
    endif()
endfunction(find_boost_from_distname_file)

if (NOT USE_BOOST_THREAD)
    set (USE_BOOST_THREAD FALSE)
endif()

if (1)
    message(STATUS "Trying to guess boost installation:")
    if (NOT Boost_INCLUDE_DIR)
        message(STATUS "no Boost_INCLUDE_DIR")
        if (BOOST_ROOT)
            message(STATUS "Found boost root: ${BOOST_ROOT}")
            get_filename_component(BOOST_ROOT "${BOOST_ROOT}" ABSOLUTE)
            set(Boost_INCLUDE_DIR ${BOOST_ROOT}/include CACHE FILEPATH "Boost inlude path")
            message(STATUS "set Boost_INCLUDE_DIR: ${Boost_INCLUDE_DIR}")
        else()
            if (NOT WIN32)
                find_package(Boost COMPONENTS python)
            endif()
            # we should think about cmake superbuild with ExternalProject_Add as the replacement for find_or_build_package (next generation)
            if (NOT Boost_INCLUDE_DIR)
                find_boost_from_distname_file()
            else()
                find_or_build_package_check(Boost boost Boost_INCLUDE_DIR False)
                find_boost_from_distname_file()
            endif()

        endif(BOOST_ROOT)
    endif(NOT Boost_INCLUDE_DIR)

    message(STATUS "BOOST_ROOT: ${BOOST_ROOT}")
    message(STATUS "Boost_INCLUDE_DIR: ${Boost_INCLUDE_DIR}")

    if (USE_BOOST_THREAD)
        find_package(Boost 1.46.0 COMPONENTS system thread REQUIRED)
        link_directories(${BOOST_LIBRARYDIR})
    endif()

    mark_as_advanced(Boost_INCLUDE_DIR BOOST_ROOT)
endif()

if (NOT USE_BOOST_THREAD)
    set(Boost_THREAD_FOUND OFF)
    set(BOOST_BIND_FOUND OFF)
    set(USE_IPC OFF)
    set(THREADS "std")
else()
    set(Boost_THREAD_FOUND ON)
    set(BOOST_BIND_FOUND ON)
    set(THREADS "boost")
endif()

# Find packages that CHOLMOD depends on
set(CMAKE_LIBRARY_PATH ${EXTERNAL_DIR}/lib $ENV{EXTERNAL_DIR}/lib ${CMAKE_LIBRARY_PATH})

if (J) # dummy to avoid error msg
endif()

set (UMFPACK_FOUND FALSE)

find_or_build_package(Triangle triangle LOCAL)
#find_package(Triangle)
set(BLA_VENDOR OpenBLAS)
find_or_build_package(BLAS lapack)
find_or_build_package(LAPACK lapack)
find_or_build_package(CHOLMOD suitesparse)
find_package(UMFPACK)

################################################################################
# Manuall check for openblas
################################################################################

#if (";${BLAS_LIBRARIES};" MATCHES "openblas")
message(STATUS "Blas found:" ${BLAS_FOUND})
if (BLAS_FOUND)
    message(STATUS "openblas is used: ${BLAS_LIBRARIES}")

    if (MINGW)
        get_filename_component(MINGW_ROOT ${BLAS_LIBRARIES} PATH)
        get_filename_component(MINGW_ROOT ${MINGW_ROOT}/../ ABSOLUTE)
        message(STATUS "MINGW_ROOT: " ${MINGW_ROOT})

        FIND_PATH(OpenBLAS_INCLUDE_DIR
            NAMES 
                cblas.h
            PATHS
                ${MINGW_ROOT}/include/OpenBLAS/
        )
    endif()
    
    if (NOT OpenBLAS_INCLUDE_DIR) 
        message(STATUS "openblas header not found")
        FIND_PATH(OpenBLAS_INCLUDE_DIR
            NAMES 
                cblas.h
            PATHS
                /usr/include/openblas/
                /usr/include/OpenBLAS/
                /mingw64/include/OpenBLAS/
            HINTS ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}
            NO_DEFAULT_PATH
        )
    endif()

    if (NOT OpenBLAS_INCLUDE_DIR) 
        message(STATUS "no openblas cblas.h in prefered paths .. searching defaults")

        FIND_PATH(OpenBLAS_INCLUDE_DIR 
            NAMES 
                cblas.h
        )
    endif()
    message(STATUS "openblas cblas.h: ${OpenBLAS_INCLUDE_DIR}")
    
    if (OpenBLAS_INCLUDE_DIR)
    set(OPENBLAS_CBLAS_FOUND TRUE)
    else()
    set(OPENBLAS_CBLAS_FOUND FALSE)
    endif()
    
    set(OPENBLAS_FOUND TRUE)
   
else()
    set(OPENBLAS_FOUND FALSE)
    set(OPENBLAS_CBLAS_FOUND FALSE)
endif()

if (NOT AVOID_CPPUNIT)
    find_package(CppUnit)
    if (CPPUNIT_FOUND)
        option (LIBGIMLI_BUILD_TESTS "Build unittests" ON)
    endif (CPPUNIT_FOUND)
endif()

if (NOT AVOID_READPROC)
if(WIN32)
    set(READPROC_FOUND FALSE)
else(WIN32)
    find_package(readproc)
endif(WIN32)
else()
    set(READPROC_FOUND FALSE)
endif()

################################################################################
# Manuall check for boost python
################################################################################
function(find_boost_python_manual bp_version)
    message(STATUS " ... searching for libboost_python${bp_version}")
    find_package(Boost 1.68 COMPONENTS "python${bp_version}")

    if (${Boost_PYTHON${bp_version}_FOUND})
        set(Boost_PYTHON_FOUND ${Boost_PYTHON${bp_version}_FOUND} CACHE BOOL "boost python lib found")
        set(Boost_PYTHON_LIBRARY ${Boost_PYTHON${bp_version}_LIBRARY_RELEASE} CACHE FILEPATH "boost python lib(find_boost_python_manual)")
        message(STATUS " ... found: ${Boost_PYTHON_LIBRARY}")
    else()
        message(STATUS " ... not found.")
    endif()
endfunction()

if (NOT Boost_PYTHON_LIBRARY)
    message(STATUS "boost-python not found by cmake ... trying to determine boost-python now manually ...")
    if (WIN32)
        set(bp_version "python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}-mt")
        string(TOUPPER ${bp_version} BP_VERSION)

        find_package(Boost 1.68 COMPONENTS "${bp_version}")
        set(Boost_PYTHON_FOUND ${Boost_${bp_version}_FOUND} CACHE BOOL "boost python found")
        set(Boost_PYTHON_LIBRARY ${Boost_${BP_VERSION}_LIBRARY} CACHE FILEPATH "boost python lib")

        
        if (NOT Boost_PYTHON_FOUND)
            set(Boost_PYTHON_LIBRARY ${BOOST_ROOT}/lib/libboost_python${bp_version}.dll CACHE FILEPATH "Boost python lib")
            message(STATUS "looking for:" ${Boost_PYTHON_LIBRARY} "...")
            if (EXISTS ${Boost_PYTHON_LIBRARY})
                message(STATUS ".. found.")
                set(Boost_PYTHON_FOUND TRUE CACHE BOOL "boost python lib found")
            else()
            message(STATUS "Sorry but no boost-python found.")
            endif()
        endif()
    else()
        # the name for py3 boost-python library lacks probably py3 suffix, which is different for different OS
        #set (Boost_PYTHON_FOUND False)
        #set (Boost_PYTHON_LIBRARY "None")
        message(STATUS "***** ${Boost_PYTHON_FOUND}, ${Boost_PYTHON_LIBRARY}")

        if (NOT Boost_PYTHON_FOUND)
            #/usr/lib/x86_64-linux-gnu/libboost_python-py34.so (debian)
            set(bp_version "-py${Python_VERSION_MAJOR}${Python_VERSION_MINOR}")
            find_boost_python_manual(${bp_version})
        endif()
        if (NOT Boost_PYTHON_FOUND)
            # /usr/lib/libboost_python3.so.1.64.0 (arch style)
            set(bp_version "${Python_VERSION_MAJOR}")
            find_boost_python_manual(${bp_version})
        endif()
        if (NOT Boost_PYTHON_FOUND)
            # /usr/lib/libboost_python-3.4.so (gentoo <boost-1.70)
            set(bp_version "-${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}")
            find_boost_python_manual(${bp_version})
        endif()
        if (NOT Boost_PYTHON_FOUND)
            # /usr/lib/libboost_python37.so # gentoo >=boost-1.70)
            set(bp_version "${Python_VERSION_MAJOR}${Python_VERSION_MINOR}")
            find_boost_python_manual(${bp_version})
        endif()
        if (NOT Boost_PYTHON_FOUND)
            # /usr/lib/libboost_python3-py37.so # weird debian)
            set(bp_version "${Python_VERSION_MAJOR}-py${Python_VERSION_MAJOR}${Python_VERSION_MINOR}")
            find_boost_python_manual(${bp_version})
        endif()

    endif()

    if (Boost_PYTHON_FOUND)
        message(STATUS "Boost_PYTHON_FOUND:" ${Boost_PYTHON_FOUND})
        message(STATUS "Boost_PYTHON_LIBRARY:" ${Boost_PYTHON_LIBRARY})
    else()
        message(STATUS "Sorry but no boost-python found. Please search for "
                       "any libboost_python* on your system and "
                       "provide it to the authors so we can try to fix "
                       "this here.")
        exit()
    endif()
endif()

if (NOT CASTER)
    set(CASTER "castxml")

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
    find_or_build_package(gccxml gccxml LOCAL)
    set(CASTER_EXECUTABLE ${GCCXML_EXECUTABLE})
else()
    if (NOT CASTER_EXECUTABLE)
        find_program(CASTXML castxml)

        if (NOT CASTXML OR CASTXML_LOCAL OR CASTXML_LOCALSRC)
            find_program(CLANG_EXECUTABLE
                NAMES clang++-3.9 clang++-3.8 clang++-3.7 clang++-3.6 clang
            )
            if(CLANG_EXECUTABLE)
                message(STATUS "clang++ found : ${CLANG_EXECUTABLE}")
            else()
                message(FATAL_ERROR "Can't found program: clang")
            endif()

            find_or_build_package(castxml castxmlbin LOCAL)
            set(CASTER_EXECUTABLE ${CASTXML_EXECUTABLE})

            if (NOT CASTXML_EXECUTABLE OR CASTXML_LOCALSRC)
                find_or_build_package(castxml castxml LOCAL)
                set(CASTER_EXECUTABLE ${CASTXML_EXECUTABLE})
            endif()
        else()
            set(CASTER_EXECUTABLE ${CASTXML})
        endif()
    endif()
endif()

find_or_build_package(pygccxml pygccxml LOCAL)
find_or_build_package(pyplusplus pygccxml LOCAL)

if (CASTER_EXECUTABLE)
    set (CASTER_FOUND TRUE)
endif()

find_package(Doxygen)
find_package(Sphinx 1.5)

if (SPHINX_FOUND AND NOT SKIP_SPHINX)
    add_subdirectory(doc EXCLUDE_FROM_ALL)
    message(STATUS "Search recursive for *.rst from here: " ${CMAKE_SOURCE_DIR})
    file(GLOB_RECURSE DOCFILES RELATIVE "${CMAKE_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}" "*.rst")
    message(STATUS "DOCFILES: " ${DOCFILES})

    file(GLOB_RECURSE PNGFILES RELATIVE "${CMAKE_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}" "*.png")
    list(APPEND DOCFILES ${PNGFILES})

    message(STATUS "PNGFILES: " ${PNGFILES})

    file(GLOB_RECURSE BIBFILES RELATIVE "${CMAKE_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}" "*.bib")
    list(APPEND DOCFILES ${BIBFILES})

    message(STATUS "BIBFILES: " ${BIBFILES})

    file(GLOB_RECURSE PYFILES RELATIVE "${CMAKE_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}" "doc/paper/*.py")
    message(STATUS "PYFILES: " ${PYFILES})

    list(APPEND DOCFILES ${PYFILES})

    set(RSTFILES ${DOCFILES} CACHE INTERNAL "RST source file that we need to copy")

    set(RSTFILES ${DOCFILES} CACHE INTERNAL "RST source file that we need to copy")
    execute_process(
        COMMAND ${SPHINX_EXECUTABLE} --version
        OUTPUT_VARIABLE SPHINX_VERSION
    )
endif()

message(STATUS "**********************************************************************")
message(STATUS "************************* Dependencies found *************************")
message(STATUS "**********************************************************************")
if (CONDA_BUILD)
    message(STATUS "CONDA_BUILD: ${CONDA_BUILD}")
    if (APPLE)
	    set(Python_LIBRARIES "-undefined dynamic_lookup")
    else()
        if(EXISTS ${Python_LIBRARIES})
            message(${Python_LIBRARIES} " exits")
        else()
            message(${Python_LIBRARIES} "does not exits")
        endif()
    endif()
endif()

message(STATUS "THREADS             : ${THREADS} ${CMAKE_THREAD_LIBS_INIT} ${Boost_THREAD_LIBRARIES}")
message(STATUS "USE_OPENBLAS        : ${OPENBLAS_FOUND} OPENBLAS_CBLAS_FOUND :${OPENBLAS_CBLAS_FOUND} OpenBLAS_INCLUDE_DIR :${OpenBLAS_INCLUDE_DIR}")
message(STATUS "CHOLMOD_LIBRARIES   : ${CHOLMOD_LIBRARIES}")
message(STATUS "UMFPACK_LIBRARIES   : ${UMFPACK_LIBRARIES}")
message(STATUS "TRIANGLE_FOUND      : ${TRIANGLE_FOUND} Triangle_LIBRARIES: ${Triangle_LIBRARIES}")
message(STATUS "Python_EXECUTABLE   : ${Python_EXECUTABLE}" )
message(STATUS "Python_LIBRARY_DIRS : ${Python_LIBRARY_DIRS}")
message(STATUS "Python_LIBRARIES    : ${Python_LIBRARIES}" )
message(STATUS "Python_INCLUDE_DIRS : ${Python_INCLUDE_DIRS}")
message(STATUS "Boost_INCLUDE_DIR   : ${Boost_INCLUDE_DIR}")
message(STATUS "Boost_PYTHON_LIBRARY: ${Boost_PYTHON_LIBRARY}" )
message(STATUS "Python_NumPy        : ${Python_NumPy_FOUND} ${Python_NumPy_INCLUDE_DIRS} ver: ${Python_NumPy_VERSION}")
message(STATUS "CASTER_FOUND        : ${CASTER_FOUND} Caster: ${CASTER_EXECUTABLE}")
message(STATUS "PYGCCXML_FOUND      : ${PYGCCXML_FOUND} PYGCCXML: ${PYGCCXML_PATH}" )
message(STATUS "PYPLUSPLUS_FOUND    : ${PYPLUSPLUS_FOUND} PYPLUSPLUS_PATH: ${PYPLUSPLUS_PATH}" )
message(STATUS "")
message(STATUS "**********************************************************************")
message(STATUS "*************** Optional Dependencies found **************************")
message(STATUS "**********************************************************************")
message(STATUS "CPPUNIT             : ${CPPUNIT_FOUND} ${CPPUNIT_LIBRARIES}")
message(STATUS "Doxygen             : ${DOXYGEN_FOUND}")
message(STATUS "Sphinx              : ${SPHINX_FOUND} ${SPHINX_EXECUTABLE} ${SPHINX_VERSION}")
message(STATUS "**********************************************************************")

if (Python_EXECUTABLE
    AND TRIANGLE_FOUND
    AND Python_LIBRARIES
    AND Boost_PYTHON_LIBRARY
    AND Python_NumPy_FOUND
    AND CASTER_FOUND
    AND PYGCCXML_FOUND
    AND PYPLUSPLUS_FOUND
)
    set (PYGIMLI 1)
    message(STATUS "")
    message(STATUS "pygimli can be build. run: make pygimli")
    message(STATUS "**********************************************************************")
    message(STATUS "")
else()
    message (FATAL_ERROR "pygimli cannot be build due to some missing packages.
Check above for python, castxml, pygccxml, pyplusplus, boost_python, numpy")
endif ()

configure_file("${PROJECT_SOURCE_DIR}/core/config.cmake.h.in"
               "${PROJECT_BINARY_DIR}/config.cmake.h" )
add_definitions(-DHAVE_CONFIG_CMAKE_H)


add_custom_target(removedoc
   COMMAND
       ${CMAKE_MAKE_PROGRAM} clean -f ${CMAKE_CURRENT_BINARY_DIR}/doc/
)
add_custom_target(checkall
    DEPENDS
      doc
      #linkcheck
)
#add_dependencies(doc removedoc)

# Add if pytest found?
add_dependencies(checkall pgtest)
if (CPPUNIT_FOUND)
    add_dependencies(checkall gtest)
endif (CPPUNIT_FOUND)


################################################################################
# Define the configurable options
################################################################################
include(CMakeDependentOption)
cmake_dependent_option( GIMLI_BUILD_SHARED_LIBS_WITH_STATIC_DEPENDENCIES
"For WIN32 and APPLE where the dependency pack is used, prefer the static dependency libs over the shared/dynamic ones.  NOTE: On Windows you also need to be mindful of which C/C++ runtime setting has been used to compile the various
components - they must all match or you will get crashes, heap corruption and/or
other issues." FALSE "WIN32 OR APPLE" FALSE)

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

set(GIMLI_VER_INSTALL_SUBDIR "/${CMAKE_PROJECT_NAME}-${GIMLI_VERSION_MAJOR}" )
set(GIMLI_MODULE_INSTALL_DIR "${GIMLI_LIB_INSTALL_DIR}/${CMAKE_PROJECT_NAME}-${GIMLI_VERSION_MAJOR}.${GIMLI_VERSION_MINOR}" )
set(GIMLI_INCLUDE_INSTALL_DIR "include${GIMLI_VER_INSTALL_SUBDIR}" )
#set( GIMLI_PYTHON_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/pygimli)

################################################################################
# Add main header locations (for everything we build)
################################################################################
include_directories(${libgimli_SOURCE_DIR}/core/src)
include_directories("${PROJECT_BINARY_DIR}")

################################################################################
# descend into subdirs
################################################################################

if (LIBGIMLI_BUILD_TESTS)
     enable_testing()
endif(LIBGIMLI_BUILD_TESTS)

add_subdirectory(core/src)
add_subdirectory(core/tests EXCLUDE_FROM_ALL)

if (PYGIMLI)
    add_subdirectory(core/python EXCLUDE_FROM_ALL)
endif (PYGIMLI)
