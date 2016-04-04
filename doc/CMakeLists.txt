if (DOXYGEN_FOUND)

    set(DOXYGEN_BUILD_DIR    ${CMAKE_CURRENT_BINARY_DIR}/doxygen )
    set(DOXYFILE_IN          ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/doxyfile.in )
    set(DOXYFILE             ${DOXYGEN_BUILD_DIR}/doxyfile )
    set(DOXY_HTML_INDEX_FILE ${DOXYGEN_BUILD_DIR}/html/index.html )
    set(DOXY_TAG_FILE        ${DOXYGEN_BUILD_DIR}/gimli.tag )

    configure_file(${DOXYFILE_IN} ${DOXYFILE})

    # perform doxygen once or if DOXYFILE_IN or libgimli has been changed
    add_custom_command(
        OUTPUT
            ${DOXY_TAG_FILE}
        COMMAND
            ${DOXYGEN_EXECUTABLE} ${DOXYFILE}
        WORKING_DIRECTORY
            ${DOXYGEN_BUILD_DIR}
        COMMENT
            "Generating gimli html API documentation with doxygen"
        DEPENDS
            gimli ${DOXYFILE} ${DOXYFILE_IN}
    )

    add_custom_target(doxygen ALL
        DEPENDS ${DOXY_TAG_FILE}
    )
endif()

set(SPHINX_SRC_DIR         ${CMAKE_SOURCE_DIR}) # configured documentation tools and intermediate build results
set(SPHINX_BUILD_DIR       ${CMAKE_CURRENT_BINARY_DIR}/_build) # HERE we export our rst. either copied or the original
set(SPHINX_SRCBUILD_DIR    ${SPHINX_BUILD_DIR}) # Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR       ${CMAKE_CURRENT_BINARY_DIR}/_doctrees) # HTML output directory
set(SPHINX_HTML_DIR        ${CMAKE_CURRENT_BINARY_DIR}/html) # PDF output directory
set(SPHINX_PDF_DIR         ${CMAKE_CURRENT_BINARY_DIR}/pdf)

configure_file(
     "${SPHINX_SRC_DIR}/doc/conf_environment.py.in"
     "${SPHINX_BUILD_DIR}/doc/conf_environment.py"
     @ONLY)

# this will ONLY cover files in the CMAKE_CURRENT_SOURCE_DIR .. maybee a cmake bug check with cmake 3
# we collect the RSTFILES in main path
# file(GLOB_RECURSE RSTFILES RELATIVE ${CMAKE_SOURCE_DIR} "*.rst")

add_custom_target(sphinx_check_sources ALL)

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

add_custom_command(
    COMMAND
        cmake -E copy_if_different ${CMAKE_SOURCE_DIR}/doc/conf.py ${SPHINX_SRCBUILD_DIR}/doc
    TARGET
        sphinx_check_sources
    VERBATIM
        COMMENT
        "Updating documentation sources files into build path: conf.py"
)

# add_custom_target(doctests
#     ALL
#         ${SPHINX_EXECUTABLE} -b doctest  ${SPHINX_SRCBUILD_DIR}
#         ${SPHINX_CACHE_DIR}/doctests
#     COMMENT "Running doctests"
#     DEPENDS
#         sphinx_check_sources
#         sphinxapi
# )

add_custom_target(sphinxapi
    ALL
        # Copy build scripts for API doc generation
    COMMAND
        cmake -E copy_if_different ${SPHINX_SRC_DIR}/doc/apigen.py ${SPHINX_SRCBUILD_DIR}/doc
    COMMAND
        cmake -E copy_if_different ${SPHINX_SRC_DIR}/doc/build_api.py ${SPHINX_SRCBUILD_DIR}/doc
    COMMAND
        cmake -E copy_if_different ${SPHINX_SRC_DIR}/doc/sidebar_gallery.py ${SPHINX_SRCBUILD_DIR}/doc
    COMMAND
        python ${SPHINX_SRCBUILD_DIR}/doc/build_api.py
    COMMENT
        "Building pygimli api doc"
    DEPENDS
        doxygen
)

add_custom_target(sphinxhtml ALL
     WORKING_DIRECTORY
        ${SPHINX_BUILD_DIR}/doc
     COMMAND
        cmake -E copy_directory ${DOXYGEN_BUILD_DIR}/html ${SPHINX_HTML_DIR}/gimliapi
     COMMAND
        ${SPHINX_EXECUTABLE} -v
        -D abort_on_example_error=1
        -T
        -b html
        -c "${SPHINX_BUILD_DIR}/doc"
        -d "${SPHINX_CACHE_DIR}/html"
        "${SPHINX_BUILD_DIR}/doc"
        "${SPHINX_HTML_DIR}"
     COMMENT
         "Building HTML documentation with sphinx"
     DEPENDS
        sphinxapi
        sphinx_check_sources
)

add_custom_target(sphinxhtml-nogallery ALL
     WORKING_DIRECTORY
        ${SPHINX_BUILD_DIR}/doc
     COMMAND
        cmake -E copy_directory ${DOXYGEN_BUILD_DIR} ${SPHINX_HTML_DIR}/gimliapi
     COMMAND
        ${SPHINX_EXECUTABLE} -v
        -D plot_gallery=0
        -T
        -b html
        -c "${SPHINX_BUILD_DIR}/doc"
        -d "${SPHINX_CACHE_DIR}/html"
        "${SPHINX_BUILD_DIR}/doc"
        "${SPHINX_HTML_DIR}"
     COMMENT
         "Building HTML documentation with sphinx"
     DEPENDS
        sphinxapi
        sphinx_check_sources
)

add_custom_target(sphinxlatex ALL
    WORKING_DIRECTORY
        ${SPHINX_BUILD_DIR}/doc
    COMMAND
        ${SPHINX_EXECUTABLE} -v
        -D plot_gallery=1
        -T
        -q
        -b latex
        -c "${SPHINX_BUILD_DIR}/doc"
        -d "${SPHINX_CACHE_DIR}/pdf"
        "${SPHINX_SRCBUILD_DIR}/doc"
        "${SPHINX_PDF_DIR}"
    COMMENT
        "Building PDF documentation with sphinx"
    DEPENDS
        sphinx_check_sources
        sphinxapi
)

add_custom_target(sphinxpdf ALL
    WORKING_DIRECTORY
        ${SPHINX_PDF_DIR}
    COMMAND
        make
    COMMAND
        cmake -E copy_if_different gimli.pdf ${SPHINX_HTML_DIR}/GIMLi_Documentation.pdf
    DEPENDS
        sphinxlatex
)

add_custom_target(html
    DEPENDS
        sphinxhtml
)

add_custom_target(doc
    COMMENT
        "Building HTML, API and PDF Documentation"
    DEPENDS
        sphinxhtml
        sphinxpdf
)
