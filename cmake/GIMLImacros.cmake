################################################################################
# Macro definitions used by GIMLI's cmake build
################################################################################
macro(add_python_module PYTHON_MODULE_NAME SOURCE_DIR EXTRA_LIBS OUTDIR)
    
    set(PYTHON_TARGET_NAME "_${PYTHON_MODULE_NAME}_")

    file( GLOB ${PYTHON_MODULE_NAME}_SOURCE_FILES ${SOURCE_DIR}/*.cpp )

    set_source_files_properties(${PYTHON_MODULE_NAME}_SOURCE_FILES
                                PROPERTIES GENERATED TRUE)

    include_directories(BEFORE ${SOURCE_DIR})
    include_directories(${PYTHON_INCLUDE_DIR})
    include_directories(${CMAKE_CURRENT_BINARY_DIR})
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/generated/)

    add_definitions(-DPYGIMLI)
    add_definitions(-DBOOST_PYTHON_NO_PY_SIGNATURES)

    add_library(${PYTHON_TARGET_NAME} MODULE ${${PYTHON_MODULE_NAME}_SOURCE_FILES})

    target_link_libraries(${PYTHON_TARGET_NAME} ${EXTRA_LIBS})
    target_link_libraries(${PYTHON_TARGET_NAME} ${PYTHON_LIBRARY})
    target_link_libraries(${PYTHON_TARGET_NAME} ${Boost_PYTHON_LIBRARY})

    set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES PREFIX "")

    if (WIN32)
        set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES SUFFIX ".pyd")
    endif()

    if (NOT APPLE AND BERT_INSTALL_WITH_RPATH)
        set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES
            INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${BERT_LIB_INSTALL_DIR}"
        )
    endif()

#     if (OUTDIR)
        set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES 
                            LIBRARY_OUTPUT_DIRECTORY_DEBUG ${OUTDIR})
        set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES 
                            LIBRARY_OUTPUT_DIRECTORY_RELEASE ${OUTDIR})
        set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES 
                            LIBRARY_OUTPUT_DIRECTORY_MINSIZEREL ${OUTDIR})
        set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES 
                            LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO ${OUTDIR})
#     endif(OUTDIR)
    
    if (CMAKE_COMPILER_IS_GNUCXX)
        set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES COMPILE_FLAGS "-fvisibility=hidden")
    endif()

    #--copy pattern files to build folder--
    set(PYTHON_IN_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
    set(PYTHON_OUT_PATH "${CMAKE_CURRENT_BINARY_DIR}")

    file(GLOB_RECURSE PYTHON_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" 
                    "${PYTHON_MODULE_NAME}/*.py" 
                    "${PYTHON_MODULE_NAME}/*.png" 
                    "${PYTHON_MODULE_NAME}/*.xrc"
                    "${PYTHON_MODULE_NAME}/*.fbp")

    
    foreach(file ${PYTHON_FILES})
        
        message ("${PYTHON_IN_PATH}/${file} ${PYTHON_OUT_PATH}/${file}")
        add_custom_command(
            OUTPUT "${PYTHON_OUT_PATH}/${file}"
            COMMAND cmake -E copy
                "${PYTHON_IN_PATH}/${file}"
                "${PYTHON_OUT_PATH}/${file}"
            DEPENDS "${PYTHON_IN_PATH}/${file}"
        )
        list( APPEND python_files_dest "${PYTHON_OUT_PATH}/${file}" )
    endforeach(file)

    add_custom_target(CopyPython ALL DEPENDS ${python_files_dest})
#----install-----------------------

    foreach(file ${PYTHON_FILES})
        get_filename_component( path_name "${file}" PATH )
        #file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pygimli/${file} DESTINATION ${path_name})
        install( FILES "${file}" DESTINATION ${path_name} )
    endforeach(file)

    install(TARGETS ${PYTHON_TARGET_NAME} LIBRARY DESTINATION "${PYTHON_MODULE_NAME}/")
endmacro()

