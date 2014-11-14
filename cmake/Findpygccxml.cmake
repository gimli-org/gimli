#
# Find the pygccxml python installation
#
# This module defines
# PYGCCXML_FOUND, If false, do not try to use pygccxml
# PYPLUSPLUS_FOUND, If false, do not try to use pyplusplus
if(NOT EXTERNAL_DIR)
	set(EXTERNAL_DIR ${PROJECT_SOURCE_DIR}/external/)
endif()

find_path(PYGCCXML_PATH pygccxml/__init__.py
    HINTS 
		${EXTERNAL_DIR}
        ${PROJECT_SOURCE_DIR}/../../pygccxml/
        /usr/lib/python2.7/site-packages/
        /usr/lib/python/site-packages/
        /usr/lib64/python2.7/site-packages/
)

find_path(PYPLUSPLUS_PATH pyplusplus/__init__.py
    HINTS 
		${EXTERNAL_DIR}
        ${PROJECT_SOURCE_DIR}/../../pyplusplus/
        /usr/lib/python/site-packages/
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

if (PYPLUSPLUS_PATH)
    message( STATUS "Found pyplusplus path: ${PYPLUSPLUS_PATH}")
    set(PYPLUSPLUS_FOUND 1)
else (PYPLUSPLUS_PATH)
    message( STATUS "NOT Found pyplusplus: cannot build pygimli.")
endif(PYPLUSPLUS_PATH)


mark_as_advanced(PYGCCXML)
mark_as_advanced(PYPLUSPLUS)