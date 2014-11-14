#
# Find the gccxml executable
#
# This module defines
# GCCXML_EXECUTABLE, If false, do not try to use gccxml

if(NOT EXTERNAL_DIR)
	set(EXTERNAL_DIR ${PROJECT_SOURCE_DIR}/external/)
endif()
#message(${GCCXML_EXECUTABLE})
#message(${EXTERNAL_DIR}/bin)
find_program(GCCXML_EXECUTABLE
  NAMES gccxml
        ../GCC_XML/gccxml
  PATHS [HKEY_CURRENT_USER\\Software\\Kitware\\GCC_XML;loc]
  		${EXTERNAL_DIR}/bin
        "$ENV{ProgramFiles}/GCC_XML"
        "C:/Program Files/GCC_XML"
		${GCCXML_ROOT}/bin
        "${PROJECT_SOURCE_DIR}/../../gccxml-bin/bin"
)

if (GCCXML_EXECUTABLE)
    message( STATUS "Found gccxml executable: ${GCCXML_EXECUTABLE}")
else()
    message( STATUS "NOT Found gccxml executable: cannot build pygimli.")
endif(GCCXML_EXECUTABLE)

mark_as_advanced(GCCXML)