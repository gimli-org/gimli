################################################################################
# Macro definitions used by GIMLI's cmake build
################################################################################
macro(add_python_module PYTHON_MODULE_NAME SOURCE_DIR EXTRA_LIBS OUTDIR)
    
    set(PYTHON_TARGET_NAME "_${PYTHON_MODULE_NAME}_")

    file(GLOB ${PYTHON_MODULE_NAME}_SOURCE_FILES ${SOURCE_DIR}/*.cpp)
    list(SORT ${PYTHON_MODULE_NAME}_SOURCE_FILES)

    set_source_files_properties(${PYTHON_MODULE_NAME}_SOURCE_FILES
                                 PROPERTIES GENERATED TRUE)

    include_directories(BEFORE ${SOURCE_DIR})
    include_directories(${PYTHON_INCLUDE_DIR})
    include_directories(${CMAKE_CURRENT_BINARY_DIR})
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/generated/)

    add_definitions(-DPYGIMLI)
    add_definitions(-DBOOST_PYTHON_NO_PY_SIGNATURES)
	add_definitions(-DBOOST_PYTHON_USE_GCC_SYMBOL_VISIBILITY)

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
        
        #message ("${PYTHON_IN_PATH}/${file} ${PYTHON_OUT_PATH}/${file}")
        add_custom_command(
            OUTPUT "${PYTHON_OUT_PATH}/${file}"
            COMMAND cmake -E copy_if_different
                "${PYTHON_IN_PATH}/${file}"
                "${PYTHON_OUT_PATH}/${file}"
            DEPENDS "${PYTHON_IN_PATH}/${file}"
        )

#         add_custom_command(
#             OUTPUT "${PYTHON_OUT_PATH}/${file}"
#             COMMAND cmake -E copy
#                 "${PYTHON_IN_PATH}/${file}"
#                 "${PYTHON_OUT_PATH}/${file}"
#             DEPENDS "${PYTHON_IN_PATH}/${file}"
#         )
        list(APPEND python_files_dest "${PYTHON_OUT_PATH}/${file}")
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

function(find_python_module module)
    string(TOUPPER ${module} module_upper)
    if(NOT PY_${module_upper})
        if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
            set(${module}_FIND_REQUIRED TRUE)
        endif()
        # A module's location is usually a directory, but for binary modules
        # it's a .so file.
        execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c" 
            "import re, ${module}; print(re.compile('\\__init__.py.*').sub('',${module}.__file__))"
            RESULT_VARIABLE _${module}_status 
            OUTPUT_VARIABLE _${module}_location
            ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

        if(NOT _${module}_status)
            set(PY_${module_upper} ${_${module}_location} CACHE STRING 
                "Location of Python module ${module}")
            message(STATUS "Find python module ${module}")
            message(STATUS "Result: ${_${module}_status}")
            message(STATUS "Output: ${_${module}_location}")
        else(NOT _${module}_status)
            message(STATUS "Find python module ${module} fails:")
            message(STATUS "Result: ${_${module}_status}")
            message(STATUS "Output: ${_${module}_location}")
            message(STATUS "Maybe you can provide the location by stetting PY_${module_upper}")
        endif(NOT _${module}_status)

        find_package_handle_standard_args(${module} 
                                      FOUND_VAR ${module}_FOUND
                                      REQUIRED_VARS PY_${module_upper}
                                      )
    else(NOT PY_${module_upper})
        set(${module}_FOUND TRUE)
    endif(NOT PY_${module_upper})

    if (${module}_FOUND)
        set( ${module}_FOUND ${${module}_FOUND} CACHE INTERNAL ${module}_FOUND)
    endif()
   
endfunction(find_python_module)

macro(find_or_build_package package get_package)

    string(TOUPPER ${package} upper_package)
    find_or_build_package_check(${package} ${get_package} ${upper_package}_FOUND)
endmacro()

macro(find_or_build_package_check package get_package checkVar)
    find_package(${package})
    
    string(TOUPPER ${package} upper_package)
    string(TOLOWER ${package} lower_package)
    
    if (NOT ${checkVar})
        message(STATUS "${package} NOT found .. building version from foreign sources into ${THIRDPARTY_DIR}" )

        file(MAKE_DIRECTORY ${THIRDPARTY_DIR})
        
        if (J)
            set(ENV{PARALLEL_BUILD} ${J})
        endif()

        if (NOT get_package)
            set(get_package ${lower_package})
        endif()

        execute_process(
            COMMAND 
				bash ${PROJECT_SOURCE_DIR}/thirdParty/buildThirdParty.sh ${get_package}
            WORKING_DIRECTORY 
				${THIRDPARTY_DIR}
        )
		
		find_package(${package})
    else()
        message(STATUS "${package} found" )
    endif()
    
endmacro()


