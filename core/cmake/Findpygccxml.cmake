#
# Find the pygccxml python installation
#
# This module defines
# PYGCCXML_FOUND, If false, do not try to use pygccxml
# PYPLUSPLUS_FOUND, If false, do not try to use pyplusplus

find_python_module(pygccxml)

if (pygccxml_FOUND)
    set(PYGCCXML_FOUND TRUE)
    set(PYGCCXML_PATH ${pygccxml_LOC})
else()
    find_path(PYGCCXML_PATH pygccxml/__init__.py
        HINTS
            ${EXTERNAL_DIR}
        NO_DEFAULT_PATH
    )

    if (PYGCCXML_PATH)
        message( STATUS "Found pygccxml path: ${PYGCCXML_PATH}")
        set(PYGCCXML_FOUND TRUE)
    else (PYGCCXML_PATH)
        message( STATUS "NOT Found pygccxml: we try to get a copy or cannot build pygimli.")
    endif(PYGCCXML_PATH)
    
endif()
    
if (PYGCCXML_PATH)
    message(STATUS "Found pygccxml path: ${PYGCCXML_PATH}")
    STRING(REGEX REPLACE "\\\\" "/" PYGCCXML_PATH ${PYGCCXML_PATH})
    message(STATUS "Found pygccxml path: ${PYGCCXML_PATH}")
endif()
mark_as_advanced(PYGCCXML)
