#
# Find the pygccxml python installation
#
# This module defines
# PYGCCXML_FOUND, If false, do not try to use pygccxml
# PYPLUSPLUS_FOUND, If false, do not try to use pyplusplus


find_python_module(pyplusplus)

if (pyplusplus_FOUND)
    set(PYPLUSPLUS_FOUND TRUE)
    set(PYPLUSPLUS_PATH ${pyplusplus_LOC})
else()

    find_path(PYPLUSPLUS_PATH pyplusplus/__init__.py
        HINTS
            ${EXTERNAL_DIR}
        NO_DEFAULT_PATH
    )

    if (PYPLUSPLUS_PATH)
        message( STATUS "Found pyplusplus path: ${PYPLUSPLUS_PATH}")
        set(PYPLUSPLUS_FOUND TRUE)
    else (PYPLUSPLUS_PATH)
        message( STATUS "NOT Found pyplusplus: we try to get a copy or cannot build pygimli.")
        find_package(Hg REQUIRED)

    endif(PYPLUSPLUS_PATH)
endif()

mark_as_advanced(PYPLUSPLUS)