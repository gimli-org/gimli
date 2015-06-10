#
# Find the gccxml executable
#
# This module defines
# GCCXML_EXECUTABLE, If false, do not try to use gccxml

find_program(GCCXML_EXECUTABLE
	NAMES
		gccxml
		../GCC_XML/gccxml
	PATHS
        ${EXTERNAL_DIR}/bin
	NO_DEFAULT_PATH
)

if (GCCXML_EXECUTABLE)
    set(GCCXML_FOUND TRUE)
    message( STATUS "Found gccxml executable: ${GCCXML_EXECUTABLE}")
else()
    message( STATUS "NOT Found gccxml executable: cannot build pygimli.")
endif(GCCXML_EXECUTABLE)

mark_as_advanced(GCCXML)
