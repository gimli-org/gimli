macro(gimli_add_python_module PYTHON_MODULE_NAME SOURCE_DIR EXTRA_LIBS )
    file( GLOB ${PYTHON_MODULE_NAME}_SOURCE_FILES ${SOURCE_DIR}/*.cpp )

    set_source_files_properties(${PYTHON_MODULE_NAME}_SOURCE_FILES
                                PROPERTIES GENERATED TRUE)

    include_directories(BEFORE ${SOURCE_DIR})
    include_directories(${PYTHON_INCLUDE_DIR})
    include_directories(${CMAKE_CURRENT_BINARY_DIR})
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/generated/)

    add_definitions(-DPYGIMLI)
    add_definitions(-DBOOST_PYTHON_NO_PY_SIGNATURES)

    add_library(${PYTHON_MODULE_NAME} MODULE ${${PYTHON_MODULE_NAME}_SOURCE_FILES})
    
    target_link_libraries(${PYTHON_MODULE_NAME} ${EXTRA_LIBS})
    target_link_libraries(${PYTHON_MODULE_NAME} ${PYTHON_LIBRARY})
    target_link_libraries(${PYTHON_MODULE_NAME} ${Boost_PYTHON_LIBRARY})
    target_link_libraries(${PYTHON_MODULE_NAME} gimli)

    set_target_properties(${PYTHON_MODULE_NAME} PROPERTIES PREFIX "")

    if (WIN32)
        set_target_properties(${PYTHON_MODULE_NAME} PROPERTIES SUFFIX ".pyd")
    endif()

    if (NOT APPLE AND GIMLI_INSTALL_WITH_RPATH)
        set_target_properties(${PYTHON_MODULE_NAME} PROPERTIES
            INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CEGUI_LIB_INSTALL_DIR}"
        )
    endif()

    if (CMAKE_COMPILER_IS_GNUCXX)
        set_target_properties(${PYTHON_MODULE_NAME} PROPERTIES COMPILE_FLAGS "-fvisibility=hidden")
    endif()

#     install(TARGETS ${PYTHON_MODULE_NAME} LIBRARY DESTINATION "${GIMLI_PYTHON_INSTALL_DIR}")
endmacro()


ADD_CUSTOM_TARGET(pygimlisrc ALL)
ADD_CUSTOM_COMMAND(
    COMMAND 
        ${PYTHON_EXECUTABLE}
    ARGS 
        ${CMAKE_CURRENT_SOURCE_DIR}/generate_pygimli_code.py
    TARGET
        pygimlisrc
    DEPENDS 
        ${CMAKE_CURRENT_SOURCE_DIR}/generate_pygimli_code.py ${CMAKE_CURRENT_SOURCE_DIR}/pygimli.h
   )

#HACK!! but we need this to force cmake looking for our generated binding files
ADD_CUSTOM_TARGET(fixpygimlisrc ALL)
ADD_CUSTOM_COMMAND(
    COMMAND 
        ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
    TARGET
        fixpygimlisrc
   )

ADD_CUSTOM_TARGET(buildpygimli ALL)
ADD_CUSTOM_COMMAND(
    COMMAND 
        make -j4 -C ${CMAKE_CURRENT_BINARY_DIR} pg
    TARGET
        buildpygimli
   )

file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/generated/dummy.cpp "")


gimli_add_python_module(_pygimli_ "${CMAKE_CURRENT_BINARY_DIR}/generated/" "" )

add_dependencies(fixpygimlisrc pygimlisrc)
add_dependencies(buildpygimli fixpygimlisrc)
add_dependencies(pygimli buildpygimli)

# for calling the build only
ADD_CUSTOM_TARGET(pg ALL)
add_dependencies(pg _pygimli_)


