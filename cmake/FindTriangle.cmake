#
# Find Jonathan Shewchuk's Triangle include directory and library, and set these variables:
#
# TRIANGLE_FOUND
# Triangle_INCLUDE_DIRS
# Triangle_LIBRARIES

# Optional user supplied search path.
set (Triangle_PREFIX_PATH "" CACHE PATH "Directory to search Triangle header and library files")
message(STATUS "Checking for package 'TRIANGLE'")

message(STATUS "${EXTERNAL_DIR}")

# Find include directory.
find_path(Triangle_INCLUDE_DIR 
            NAMES
                triangle.h 
            PATHS
                ${Triangle_PREFIX_PATH}
                ${EXTERNAL_DIR}
                ${PROJECT_SOURCE_DIR}/external
                ${PROJECT_BINARY_DIR}/external
                /usr/local
                /usr
            PATH_SUFFIXES
                include
            )

find_library(Triangle_LIBRARIES
            NAMES
                triangle
            PATHS
                ${Triangle_PREFIX_PATH} 
                ${EXTERNAL_DIR}
                ${PROJECT_SOURCE_DIR}/external
                ${PROJECT_BINARY_DIR}/external
            PATH_SUFFIXES
                lib
            )

# Standard package handling
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Triangle 
                                    "Triangle could not be found."
                                    Triangle_LIBRARIES Triangle_INCLUDE_DIR
                                 )

mark_as_advanced(
    Triangle_INCLUDE_DIR
    Triangle_LIBRARIES
    )
