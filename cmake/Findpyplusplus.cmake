#
# Find the pygccxml python installation
#
# This module defines
# PYGCCXML_FOUND, If false, do not try to use pygccxml
# PYPLUSPLUS_FOUND, If false, do not try to use pyplusplus

find_path(PYPLUSPLUS_PATH pyplusplus/__init__.py
    HINTS 
		${EXTERNAL_DIR}
        ${PROJECT_SOURCE_DIR}/../../pyplusplus/
    NO_DEFAULT_PATH
)

if (PYPLUSPLUS_PATH)
    message( STATUS "Found pyplusplus path: ${PYPLUSPLUS_PATH}")
    set(PYPLUSPLUS_FOUND 1)
else (PYPLUSPLUS_PATH)
    message( STATUS "NOT Found pyplusplus: cannot build pygimli.")
endif(PYPLUSPLUS_PATH)

mark_as_advanced(PYPLUSPLUS)