cmake_minimum_required(VERSION 2.8.6)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Sets the configuration to build (Release, Debug, etc...)")
endif()

project(libgimli)

set(LIBGIMLI_VERSION_MAJOR 0)
set(LIBGIMLI_VERSION_MINOR 9)
set(LIBGIMLI_VERSION_PATCH 1)
set(LIBGIMLI_VERSION ${LIBGIMLI_VERSION_MAJOR}.${LIBGIMLI_VERSION_MINOR}.${LIBGIMLI_VERSION_PATCH}) 

################################################################################
# Check for libs and other packages we might use.
################################################################################

find_package(Boost 1.48.0 COMPONENTS python system thread)
find_package(Doxygen)

################################################################################
# Define the configurable options
################################################################################
include(CMakeDependentOption)
cmake_dependent_option( CEGUI_BUILD_SHARED_LIBS_WITH_STATIC_DEPENDENCIES 
"For WIN32 and APPLE where the dependency pack is used, prefer the static dependency libs over the shared/dynamic ones.  NOTE: On Windows you also need to be mindful of which C/C++ runtime setting has been used to compile the various
components - they must all match or you will get crashes, heap corruption and/or
other issues." FALSE "WIN32 OR APPLE" FALSE)

################################################################################
# Add main header locations (for everything we build)
################################################################################
include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}/src
	${CMAKE_CURRENT_BINARY_DIR}/apps
)

################################################################################
# descend into subdirs
################################################################################
add_subdirectory(src)

if (DOXYGEN_FOUND)
    add_custom_target(html
        "${DOXYGEN_EXECUTABLE}"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen"
        COMMENT "Generating documentation" VERBATIM)
endif()

if (LIBGIMLI_BUILD_TESTS)
    enable_testing()

    add_subdirectory(tests)
endif()


