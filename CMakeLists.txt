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

# Check if conda package is created
if(DEFINED ENV{CONDA_BUILD})
    message(STATUS "Conda package is being created.")
    set(CONDA_BUILD TRUE)
    set(Boost_INCLUDE_DIR "${CONDA_PREFIX}/include")
else()
    set(CONDA_BUILD FALSE)
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -pipe -ansi -Wall -Wno-long-long -Wno-unused-result -Wno-unused-variable")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-unused-value -Wno-strict-aliasing -Wno-unused-local-typedefs")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-int-in-bool-context")

    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -p -O2 -ansi -pedantic -fno-omit-frame-pointer -Wall")

    if (WIN32 OR CONDA_BUILD)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-deprecated-declarations")
    endif()

    if (NOT WIN32 AND ASAN)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
    endif()

    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-Os")
#     if (NOT WIN32)
#         add_definitions(-std=c++0x)
    add_definitions(-std=c++11)
    
    # add_definitions(-std=gnu++0x)
    # -Wl,--no-undefined
    # set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-Os -static-libgcc -Wl,-O2 -Wl,--as-needed -Wl,--sort-common")
    # set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-Os -static-libgcc")
    # set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-Os -Wl,--no-undefined -static-libgcc -Wl,-O2 -Wl,--as-needed -Wl,--sort-common")
    # set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS} -lbfd")

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-overloaded-virtual")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Qunused-arguments")

        if (APPLE)
            # If it's libc++ and you're on <= 10.8, you need to compile with clang++ -stdlib=libc++. If it's libstdc++ and you're on 10.9 or later, you need to compile with clang++ -stdlib=libstdc++.
            # add_compile_options(-stdlib=libstdc++)
            # add_compile_options(-stdlib=libc++)
        endif()
    else() # if gcc

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
    set(LIBGIMLI_VERSION "untaged")
endif()

set( PACKAGE_NAME  \"${PROJECT_NAME}\" )
set( PACKAGE_VERSION  \"${LIBGIMLI_VERSION}\" )
set( PACKAGE_BUGREPORT  \"carsten@pygimli.org\")
set( PACKAGE_AUTHORS  \"carsten@pygimli.org thomas@pygimli.org florian@pygimli.org\")

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
find_package(Threads REQUIRED)

if (PYVERSION)
    find_package(PythonInterp ${PYVERSION} REQUIRED)
else()
    find_package(PythonInterp)
endif()

# if (PYVERSION)
#     set(Python_ADDITIONAL_VERSIONS ${PYVERSION})
# endif()
# find_package(PythonInterp) # find version belongs to the active python command
#
# if(WIN32)
#     find_package(PythonLibs)
# endif(WIN32)

function(find_boost_from_distname_file)
    set (BOOST_DIST_NAME_FILE "${EXTERNAL_DIR}/.boost-py${PYTHON_VERSION_MAJOR}.dist")

    if (EXISTS ${BOOST_DIST_NAME_FILE})
        message(STATUS "reading ${BOOST_DIST_NAME_FILE}")
        FILE(READ "${EXTERNAL_DIR}/.boost-py${PYTHON_VERSION_MAJOR}.dist" BOOST_DIST_NAME)
        STRING(REGEX REPLACE "\n" "" BOOST_DIST_NAME "${BOOST_DIST_NAME}")
        get_filename_component(BOOST_ROOT "${EXTERNAL_DIR}/${BOOST_DIST_NAME}" ABSOLUTE)
    else()
        message(STATUS "there is no ${BOOST_DIST_NAME_FILE}")
    endif()

    if (EXISTS ${BOOST_ROOT})
        message(STATUS "BOOST_ROOT from ${EXTERNAL_DIR}/.boost-py${PYTHON_VERSION_MAJOR}.dist: ${BOOST_ROOT}")

        set(BOOST_ROOT ${BOOST_ROOT}  CACHE FILEPATH "Boost root distribution")

        set(Boost_INCLUDE_DIR ${BOOST_ROOT}/include)
        get_filename_component(Boost_INCLUDE_DIR "${Boost_INCLUDE_DIR}" ABSOLUTE)
        message(STATUS "Boost_INCLUDE_DIR, ${Boost_INCLUDE_DIR}")
    else()
        message(STATUS "BOOST_ROOT does not exist: ${BOOST_ROOT}")
    endif()
endfunction(find_boost_from_distname_file)

if (NOT USE_BOOST_THREAD)
    if (NOT WIN32)
        set (USE_BOOST_THREAD FALSE)
    else()
        set (USE_BOOST_THREAD TRUE)
    endif()
endif()

if (USE_BOOST_THREAD)
        message(STATUS "Trying to guess boost installation:")
    if (NOT Boost_INCLUDE_DIR)
        message(STATUS "no BOOST_INCLUDE_DIR")
        if (BOOST_ROOT)
            message(STATUS "Found boost root: ${BOOST_ROOT}")
            get_filename_component(BOOST_ROOT "${BOOST_ROOT}" ABSOLUTE)

            set(Boost_INCLUDE_DIR ${BOOST_ROOT}/include)
            message(STATUS "set BOOST_INCLUDE_DIR: ${BOOST_INCLUDE_DIR}")
        else()
            if (NOT WIN32)
                find_package(Boost 1.46.0 COMPONENTS system thread)
            endif()
            # we should think about cmake superbuild with ExternalProject_Add as the replacement for find_or_build_package (next generation)
            if (NOT Boost_INCLUDE_DIR)
                find_boost_from_distname_file()
            endif()
            find_or_build_package_check(Boost boost Boost_INCLUDE_DIR False)
            find_boost_from_distname_file()

            #exit()
        endif(BOOST_ROOT)
        message(STATUS "BOOST_ROOT, ${BOOST_ROOT}")
    endif(NOT Boost_INCLUDE_DIR)

    message(STATUS "Boost_INCLUDE_DIR, ${Boost_INCLUDE_DIR}")

    if (USE_BOOST_THREAD)
        find_package(Boost 1.46.0 COMPONENTS system thread REQUIRED)
        link_directories(${BOOST_LIBRARYDIR})
    endif()

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
endif()

if (NOT USE_BOOST_THREAD)
    set(Boost_THREAD_FOUND OFF)
    set(BOOST_BIND_FOUND OFF)
    set(USE_IPC OFF)
else()
    set(Boost_THREAD_FOUND ON)
    set(BOOST_BIND_FOUND ON)
endif()

# Find packages that CHOLMOD depends on
set(CMAKE_LIBRARY_PATH ${EXTERNAL_DIR}/lib $ENV{EXTERNAL_DIR}/lib ${CMAKE_LIBRARY_PATH})

if (J) # dummy to avoid error msg
endif()

set (UMFPACK_FOUND FALSE)

find_or_build_package(Triangle triangle LOCAL)
find_or_build_package(BLAS lapack)
find_or_build_package(LAPACK lapack)
find_or_build_package(CHOLMOD suitesparse)
find_package(UMFPACK)

if (";${BLAS_LIBRARIES};" MATCHES "openblas")
    message(STATUS "openblas is used: ${BLAS_LIBRARIES}")
    message(STATUS "openblas is used: ${PC_BLAS_INCLUDE_DIRS}")
    set(OPENBLAS_FOUND TRUE)
else()
    set(OPENBLAS_FOUND FALSE)
endif()

#find_package(openblas)

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
# Check for python stuff
################################################################################
if (WIN32)
    #if (NOT PYTHON_LIBRARY)
    #   get_filename_component(PYTHONPATH "${PYTHON_EXECUTABLE}" PATH )
    #    if (${ADDRESSMODEL} EQUAL 32 OR ${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} GREATER 34) # works at least for WinPython3.5
    #        set (PYTHON_LIBRARY ${PYTHONPATH}/libs/libpython${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.a)
    #    else()
    #        set (PYTHON_LIBRARY ${PYTHONPATH}/libs/libpython${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.dll.a)
    #    endif()
    # endif()
else()
	set(Python_ADDITIONAL_VERSIONS ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
	set(PythonLibs_FIND_VERSION ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
endif()

find_package(PythonLibs)

find_python_module(numpy)

if (Boost_PYTHON_LIBRARY)
    #set (Boost_PYTHON_LIBRARY ${Boost_PYTHON_LIBRARY})
else()
    #find_package(Boost 1.46.0 COMPONENTS python)

    # the name for py3 boost-python library lacks probably py3 suffix, which is different for different OS

    #debian /usr/lib/x86_64-linux-gnu/libboost_python-py34.so
    set(bp_version "${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
    find_package(Boost 1.46 COMPONENTS "python-py${bp_version}")
    if (${Boost_PYTHON-PY${bp_version}_FOUND})
        set(Boost_PYTHON_FOUND ${Boost_PYTHON-PY${bp_version}_FOUND})
        set(Boost_PYTHON_LIBRARY ${Boost_PYTHON-PY${bp_version}_LIBRARY_RELEASE})
    endif()

    if (NOT Boost_PYTHON_FOUND)
        # arch style /usr/lib/libboost_python3.so.1.64.0
        set(bp_version "${PYTHON_VERSION_MAJOR}")
        find_package(Boost 1.46 COMPONENTS "python${bp_version}")
        if (${Boost_PYTHON${bp_version}_FOUND})
            set(Boost_PYTHON_FOUND ${Boost_PYTHON${bp_version}_FOUND})
            set(Boost_PYTHON_LIBRARY ${Boost_PYTHON${bp_version}_LIBRARY_RELEASE})
        endif()
    endif()

    if (NOT Boost_PYTHON_FOUND)
        # gentoo: /usr/lib/libboost_python-3.4.so
        # should find the right one
        find_package(Boost 1.46 COMPONENTS python)
    endif()

    if (NOT Boost_PYTHON_FOUND)
        # gentoo: /usr/lib/libboost_python-3.4.so
        set(bp_version "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")
        # should find the right one
        find_package(Boost 1.46 COMPONENTS "python-${bp_version}")
        if (${Boost_PYTHON-${bp_version}_FOUND})
            set(Boost_PYTHON_FOUND ${Boost_PYTHON-${bp_version}_FOUND})
            set(Boost_PYTHON_LIBRARY ${Boost_PYTHON-${bp_version}_LIBRARY_RELEASE})
        endif()
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

if (PYTHON_EXECUTABLE)
    set (PYTHON_FOUND TRUE)
endif()

if (Boost_PYTHON_LIBRARY)
    set (Boost_PYTHON_FOUND TRUE)
endif()

find_package(Doxygen)

find_package(Sphinx 1.5)

if (SPHINX_FOUND)
    add_subdirectory(doc EXCLUDE_FROM_ALL)
    file(GLOB_RECURSE DOCFILES RELATIVE "${CMAKE_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}" "*.rst")
    file(GLOB_RECURSE PNGFILES RELATIVE "${CMAKE_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}" "*.png")
    list(APPEND DOCFILES ${PNGFILES})
    file(GLOB_RECURSE BIBFILES RELATIVE "${CMAKE_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}" "*.bib")
    list(APPEND DOCFILES ${BIBFILES})
    file(GLOB_RECURSE PYFILES RELATIVE "${CMAKE_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}" "doc/paper/*.py")
    list(APPEND DOCFILES ${PYFILES})

    set(RSTFILES ${DOCFILES} CACHE INTERNAL "RST source file that we need to copy")

    #foreach(file ${DOCFILES})
        #message(${file})
    #endforeach()

endif()

message(STATUS "**********************************************************************")
message(STATUS "************************* Dependencies found *************************")
message(STATUS "**********************************************************************")
if (CONDA_BUILD)
    message(STATUS "CONDA_BUILD: ${CONDA_BUILD}")
endif()
message(STATUS "THREADS            :${Threads_FOUND}    ${CMAKE_THREAD_LIBS_INIT}")
message(STATUS "USE_BOOST_THREAD   :${USE_BOOST_THREAD} ${Boost_THREAD_LIBRARIES}")
message(STATUS "CHOLMOD_LIBRARIES  :${CHOLMOD_LIBRARIES}")
message(STATUS "UMFPACK_LIBRARIES  :${UMFPACK_LIBRARIES}")
message(STATUS "TRIANGLE_FOUND     :${TRIANGLE_FOUND}    Triangle_LIBRARIES: ${Triangle_LIBRARIES}")
message(STATUS "PYTHON_FOUND       :${PYTHON_FOUND}    PYTHON_EXECUTABLE: ${PYTHON_EXECUTABLE}" )
message(STATUS "PYTHONLIBS_FOUND   :${PYTHONLIBS_FOUND}    PYTHON_LIBRARY: ${PYTHON_LIBRARY}" )
if(WIN32)
message(STATUS "PYTHON_INCLUDE_DIR :${PYTHON_INCLUDE_DIR}" )
endif(WIN32)
message(STATUS "Boost_INCLUDE_DIR  :${Boost_INCLUDE_DIR}")
message(STATUS "Boost_PYTHON_FOUND :${Boost_PYTHON_FOUND}    Boost_PYTHON_LIBRARY: ${Boost_PYTHON_LIBRARY}" )
message(STATUS "numpy_FOUND        :${numpy_FOUND}    PY_NUMPY: ${PY_NUMPY}" )
message(STATUS "CASTER_FOUND       :${CASTER_FOUND}    Caster: ${CASTER_EXECUTABLE}")
message(STATUS "PYGCCXML_FOUND     :${PYGCCXML_FOUND}    PYGCCXML: ${PYGCCXML_PATH}" )
message(STATUS "PYPLUSPLUS_FOUND   :${PYPLUSPLUS_FOUND}    PYPLUSPLUS_PATH: ${PYPLUSPLUS_PATH}" )
message(STATUS "")
message(STATUS "**********************************************************************")
message(STATUS "*************** Optional Dependencies found **************************")
message(STATUS "**********************************************************************")
message(STATUS "CPPUNIT            :${CPPUNIT_FOUND}    ${CPPUNIT_LIBRARIES}")
message(STATUS "Doxygen            :${DOXYGEN_FOUND}    ")
message(STATUS "Sphinx             :${SPHINX_FOUND}    ")
message(STATUS "**********************************************************************")

if (PYTHON_EXECUTABLE
    AND TRIANGLE_FOUND
    AND PYTHON_FOUND
    AND PYTHONLIBS_FOUND
    AND Boost_PYTHON_FOUND
    AND numpy_FOUND
    AND CASTER_FOUND
    AND PYGCCXML_FOUND
    AND PYPLUSPLUS_FOUND
)
    set (PYGIMLI 1)
    message(STATUS "")
    message (STATUS "pygimli can be build. run: make pygimli")
    message(STATUS "**********************************************************************")
    message(STATUS "")
    message(STATUS "")
else()
    message (FATAL_ERROR "pygimli cannot be build due to some missing packages.
Check above for python, castxml, pygccxml, pyplusplus, boost_python, numpy")
endif ()

configure_file("${PROJECT_SOURCE_DIR}/config.cmake.h.in" "${PROJECT_BINARY_DIR}/config.cmake.h" )
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
add_subdirectory(apps EXCLUDE_FROM_ALL)

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
