#
# Find the gccxml executable
#
# This module defines
# GCCXML_EXECUTABLE, If false, do not try to use gccxml

find_program(CASTXML_EXECUTABLE
	NAMES
		castxml
	PATHS
        ${EXTERNAL_DIR}/bin
	NO_DEFAULT_PATH
)

if (CASTXML_EXECUTABLE)
    set(CASTXML_FOUND TRUE)
    message( STATUS "Found castxml executable: ${GCCXML_EXECUTABLE}")
else()
    message( STATUS "NOT Found castxml executable: cannot build pygimli.")
endif(CASTXML_EXECUTABLE)

mark_as_advanced(CASTXML)
