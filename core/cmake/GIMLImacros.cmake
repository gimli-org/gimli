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
    include_directories(${Python_INCLUDE_DIRS})
    include_directories(${CMAKE_CURRENT_BINARY_DIR})
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/generated/)

    add_definitions(-DPYGIMLI)
    add_definitions(-DBOOST_PYTHON_NO_PY_SIGNATURES)
	add_definitions(-DBOOST_PYTHON_USE_GCC_SYMBOL_VISIBILITY)

    add_library(${PYTHON_TARGET_NAME} MODULE ${${PYTHON_MODULE_NAME}_SOURCE_FILES})

    target_link_libraries(${PYTHON_TARGET_NAME} ${EXTRA_LIBS})
    target_link_libraries(${PYTHON_TARGET_NAME} ${Python_LIBRARIES})
    target_link_libraries(${PYTHON_TARGET_NAME} ${Boost_PYTHON_LIBRARY})

    set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES PREFIX "")

    if (WIN32)
        set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES SUFFIX ".pyd")
    endif()

    #if (NOT APPLE AND BERT_INSTALL_WITH_RPATH)
    #    set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES
    #        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${BERT_LIB_INSTALL_DIR}"
    #    )
    #endif()

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
        set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES COMPILE_FLAGS "-fvisibility=hidden -Wno-unused-value")
        if (WIN32 AND ADDRESSMODEL EQUAL "64")
            set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES DEFINE_SYMBOL "MS_WIN64")
        endif()
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        #     using regular Clang or AppleClang
        set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES COMPILE_FLAGS "-fvisibility=hidden -Wno-unused-value")
    endif()

    #--copy pattern files to build folder--
    set(PYTHON_IN_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
    set(PYTHON_OUT_PATH "${CMAKE_BINARY_DIR}/package")

    file(GLOB_RECURSE PYTHON_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
                    "${PYTHON_MODULE_NAME}/*.py"
                    "${PYTHON_MODULE_NAME}/*.png"
                    "${PYTHON_MODULE_NAME}/*.xrc"
                    "${PYTHON_MODULE_NAME}/*.fbp")

    add_custom_target(copy_python ALL)

    foreach(file ${PYTHON_FILES})

        #message ("${PYTHON_IN_PATH}/${file} ${PYTHON_OUT_PATH}/${file}")
        add_custom_command(
            COMMAND
                cmake -E copy_if_different
                ${PYTHON_IN_PATH}/${file}
                ${PYTHON_OUT_PATH}/${file}
            DEPENDS "${PYTHON_IN_PATH}/${file}"
            TARGET
                copy_python
            VERBATIM
            COMMENT
                "Updating python file: ${file}"
        )
    endforeach(file)


#----install-----------------------
#     foreach(file ${PYTHON_FILES})
#         get_filename_component( path_name "${file}" PATH )
#         #file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pygimli/${file} DESTINATION ${path_name})
#         install( FILES "${file}" DESTINATION ${path_name} )
#     endforeach(file)
#
#     install(TARGETS ${PYTHON_TARGET_NAME} LIBRARY DESTINATION "${PYTHON_MODULE_NAME}/")
endmacro()

function(find_python_module module)
    string(TOUPPER ${module} module_upper)
    if(NOT PY_${module_upper})
        if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
            set(${module}_FIND_REQUIRED TRUE)
        endif()
        # A module's location is usually a directory, but for binary modules
        # it's a .so file.
        execute_process(COMMAND "${Python_EXECUTABLE}" "-c"
            "import re, ${module}; print(re.compile('\\__init__.py.*').sub('',${module}.__file__))"
            RESULT_VARIABLE _${module}_status
            OUTPUT_VARIABLE _${module}_location
            OUTPUT_STRIP_TRAILING_WHITESPACE)

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
        set( ${module}_LOC ${_${module}_location} CACHE INTERNAL ${module}_LOC)
    endif()

endfunction(find_python_module)

macro(findBuildTools)
    message(STATUS "checking for some build tools ...")
    #unzip try cmake -E tar
    #find_package(Tar REQUIRED)  ${CMAKE_COMMAND} -E tar "cfvz"
    find_program(PATCH_TOOL NAMES patch  REQUIRED)
    find_program(SED_TOOL NAMES sed REQUIRED)
    find_package(Wget REQUIRED)
    find_package(Git REQUIRED)
endmacro(findBuildTools)

macro(find_or_build_package package get_package)

    set (extraMacroArgs ${ARGN})

    # Did we get any optional args?
    list(LENGTH extraMacroArgs numExtraArgs)
    if (${numExtraArgs} GREATER 0)
        list(GET extraMacroArgs 0 optionalArg)
        set(foceLocal True)
    else()
        set(foceLocal False)
    endif()

    string(TOUPPER ${package} upper_package)
    find_or_build_package_check(${package} ${get_package} ${upper_package}_FOUND ${foceLocal})
endmacro()

macro(find_or_build_package_check package get_package checkVar forceLocal)

    message(STATUS "** Find or build ${package} at: ${checkVar}")
    find_package(${package})
    message(STATUS "Found: ${${package}_FOUND}")

    string(TOUPPER ${package} upper_package)
    string(TOLOWER ${package} lower_package)

    set (FORCE_LOCAL_REBUILD 0)

    message(STATUS "Local build ${package} forced: ${forceLocal}")

    if ($ENV{CLEAN})
        if(${forceLocal} OR ${package}_LOCAL)
            set(FORCE_LOCAL_REBUILD 1)
            set(ENV{CLEAN} 1)
            message(STATUS "Rebuild forced for: ${package}")
        endif()
    endif()

    if (NOT ${checkVar} OR ${FORCE_LOCAL_REBUILD})

        findBuildTools()

        message(STATUS "building ${package} from foreign sources into ${THIRDPARTY_DIR}" )

        file(MAKE_DIRECTORY ${THIRDPARTY_DIR})

        if (J)
            set(ENV{PARALLEL_BUILD} ${J})
        endif()

        if (NOT get_package)
            set(get_package ${lower_package})
        endif()

        execute_process(
            COMMAND
				bash ${PROJECT_SOURCE_DIR}/core/scripts/buildThirdParty.sh ${get_package}
            WORKING_DIRECTORY
				${THIRDPARTY_DIR}
        )

        message(STATUS "checking again for ${package} ...")
		find_package(${package})
        message(STATUS "Found: ${${package}_FOUND}")

        set(${package}_LOCAL 1 CACHE INTERNAL "this package was build local")
    else()
        message(STATUS "** Find or build ${package} done.")
    endif()

endmacro()


