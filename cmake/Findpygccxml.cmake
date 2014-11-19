#
# Find the pygccxml python installation
#
# This module defines
# PYGCCXML_FOUND, If false, do not try to use pygccxml
# PYPLUSPLUS_FOUND, If false, do not try to use pyplusplus

find_path(PYGCCXML_PATH pygccxml/__init__.py
    HINTS 
		${EXTERNAL_DIR}
        ${PROJECT_SOURCE_DIR}/../../pygccxml/
        /usr/lib/python2.7/site-packages/
        /usr/lib/python/site-packages/
        /usr/lib64/python2.7/site-packages/
)

if (PYGCCXML_PATH)
    message( STATUS "Found pygccxml path: ${PYGCCXML_PATH}")
    set(PYGCCXML_FOUND 1)
else (PYGCCXML_PATH)
    message( STATUS "NOT Found pygccxml: cannot build pygimli.")
endif(PYGCCXML_PATH)

mark_as_advanced(PYGCCXML)
