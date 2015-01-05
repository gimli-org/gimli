
if (DOXYGEN_FOUND)
    set(DOXYGEN_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/doxygen")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/doxygen/doxyfile.in"
                    ${DOXYGEN_BUILD_DIR}/doxyfile)
    add_custom_target(doxygen
        "${DOXYGEN_EXECUTABLE}"
        WORKING_DIRECTORY "${DOXYGEN_BUILD_DIR}"
        COMMENT "Generating API documentation with doxygen" VERBATIM)

    add_dependencies(doxygen gimli)
endif()


set(SPHINX_SRC_DIR "${CMAKE_SOURCE_DIR}")

# configured documentation tools and intermediate build results
set(SPHINX_BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")

# HERE we excpect our rst. either copied or the original
set(SPHINX_SRCBUILD_DIR "${SPHINX_BINARY_BUILD_DIR}")

# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

# HTML output directory
set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

configure_file(
     "${SPHINX_SRC_DIR}/doc/conf.py.in"
     "${SPHINX_BINARY_BUILD_DIR}/conf.py"
     @ONLY)

configure_file(
     "${SPHINX_SRC_DIR}/doc/sidebar_gallery.py.in"
     "${SPHINX_BINARY_BUILD_DIR}/sidebar_gallery.py"
     @ONLY)
# this will ONLY cover files in the CMAKE_CURRENT_SOURCE_DIR .. maybee a cmake bug check with cmake 3
# we collect the RSTFILES in main path
# file(GLOB_RECURSE RSTFILES RELATIVE ${CMAKE_SOURCE_DIR} "*.rst")

ADD_CUSTOM_TARGET(sphinx_check_sources ALL)

foreach(file ${RSTFILES})
    add_custom_command(
            COMMAND
                cmake -E copy_if_different ${CMAKE_SOURCE_DIR}/${file} ${SPHINX_SRCBUILD_DIR}/${file}
            TARGET
                sphinx_check_sources
            VERBATIM
            COMMENT 
                "Updating documentation sources files into build path: ${file}"
    )
endforeach()

# add_custom_command(
#     COMMAND
#         python ${SPHINX_SRCBUILD_DIR}/doc/sidebar_gallery.py
#     TARGET
#         sphinx_check_sources
#     VERBATIM
#     COMMENT 
#         "calling: python sidebar_gallery.py"
# )

add_custom_target(sphinx
    ALL
        ${SPHINX_EXECUTABLE}
            -q
            -b html
            -c "${SPHINX_BINARY_BUILD_DIR}"
            -d "${SPHINX_CACHE_DIR}"
            "${SPHINX_SRCBUILD_DIR}"
            "${SPHINX_HTML_DIR}"
    COMMENT "Building HTML documentation with Sphinx")

#add_dependencies(sphinx pygimli)
#add_dependencies(sphinx doxygen)
add_dependencies(sphinx sphinx_check_sources)

add_custom_target(doc)
add_dependencies(doc doxygen)
add_dependencies(doc sphinx)
