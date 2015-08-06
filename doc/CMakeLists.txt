if (DOXYGEN_FOUND)
    set( DOXYGEN_BUILD_DIR    ${CMAKE_CURRENT_BINARY_DIR}/doxygen )
    set( DOXYFILE_IN          ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/doxyfile.in )
    set( DOXYFILE             ${DOXYGEN_BUILD_DIR}/doxyfile )
    set( DOXY_HTML_INDEX_FILE ${DOXYGEN_BUILD_DIR}/html/index.html )
    set( DOXY_TAG_FILE        ${DOXYGEN_BUILD_DIR}/gimli.tag )

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

set( SPHINX_SRC_DIR         ${CMAKE_SOURCE_DIR}) # configured documentation tools and intermediate build results
set( SPHINX_BUILD_DIR       ${CMAKE_CURRENT_BINARY_DIR}/_build) # HERE we export our rst. either copied or the original

############################################################
#  The commands below build the docs in SPHINX_BUILD_DIR.i #
#       No copying of source files is necessary!           #
############################################################

add_custom_target(sphinxapi
    ALL
        cd ${SPHINX_SRC_DIR}/doc && make api
    COMMENT
        "Building pygimli api doc"
    DEPENDS
        doxygen
)

add_custom_target(sphinxpdf ALL
    COMMAND
    cd ${SPHINX_SRC_DIR}/doc && make pdf BUILDDIR=${SPHINX_BUILD_DIR} &&
    cd ${SPHINX_BUILD_DIR}/pdf && cp gimli.pdf ../html/GIMLi_Documentation.pdf
    COMMENT
        "Building PDF documentation with sphinx"
    DEPENDS
        sphinxapi
)

add_custom_target(sphinxhtml ALL
    COMMAND
        cd ${SPHINX_SRC_DIR}/doc && make html BUILDDIR=${SPHINX_BUILD_DIR}
    COMMENT
        "Building HTML documentation with sphinx"
    DEPENDS
        sphinxapi
        sphinxpdf
)

add_custom_target(doc
    COMMAND
        # redundant, but necessary for now. Builds html pages again with API and
        # sidebar gallery without running all examples again.
        cd ${SPHINX_SRC_DIR}/doc
        && make api
        && make html-nogallery BUILDDIR=${SPHINX_BUILD_DIR}
    COMMENT
        "Building HTML, API and PDF Documentation"
    DEPENDS
        sphinxhtml
)

add_custom_target(sphinxclean
    COMMAND
        cd ${SPHINX_SRC_DIR}/doc && make clean
    COMMENT
        "Cleaning PDF, HTML and API Documentation"
)

