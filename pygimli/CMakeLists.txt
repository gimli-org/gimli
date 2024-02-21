
cprint(green "** START prepare whl build")

set(TARGET_NAME pygimli)
file(GLOB_RECURSE PYGIMLI_FILES 
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
    ${CMAKE_CURRENT_SOURCE_DIR}/*.py)

message(STATUS "Copying python files to ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}/${TARGET_NAME}")
foreach (filename ${PYGIMLI_FILES})
    #message(STATUS "Copying ${filename} to ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}/${TARGET_NAME}/${target}")
    configure_file("${filename}" "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}/${filename}" COPYONLY)
endforeach (filename)

CONFIGURE_FILE(setup.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/setup.py @ONLY
               )
CONFIGURE_FILE(getwheelname.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/getwheelname.py @ONLY
               )
execute_process(COMMAND 
                    "${Python_EXECUTABLE}" getwheelname.py 
                RESULT_VARIABLE STATUS
                OUTPUT_VARIABLE PY_WHEELFILE
                WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
                )
                
if (STATUS EQUAL 1)
    message( FATAL_ERROR "Can't determine wheel filenname. Is wheel installed? python -m pip install wheel")
endif()

message(STATUS "Detected Python wheel filename: ${PY_WHEELFILE}")

add_custom_target(PyPackageCopyLibs)

add_dependencies(PyPackageCopyLibs copy_libs)

add_custom_target(copy_libs)
    foreach(file 
            $<TARGET_FILE:${Boost_PYTHON_LIBRARY}>
            $<TARGET_FILE:gimli>
            $<TARGET_FILE:_pygimli_>
            )
    add_custom_command(
            COMMAND
                cmake -E copy_if_different
                ${file}
                ${CMAKE_BINARY_DIR}/${PYTHON_MODULE_NAME}/${TARGET_NAME}/${TARGET_NAME}/core/libs
            DEPENDS ${file}
            TARGET copy_libs
            VERBATIM
            COMMENT
                "Copy libraries: ${file} ${CMAKE_BINARY_DIR}/${TARGET_NAME}/${TARGET_NAME}/core/libs"
    )
    endforeach(file)

add_custom_target(PyPackageBuild
                COMMAND 
                    "${Python_EXECUTABLE}" -m pip wheel .
                WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
                COMMENT 
                    "Building python wheel package"
)

add_custom_target(whl)
add_dependencies(PyPackageBuild PyPackageCopyLibs)
add_dependencies(whl PyPackageBuild)

add_custom_target(PyPackageWhlUpload
                COMMAND 
                    "${Python_EXECUTABLE}" -m twine upload --repository testpypi ${PY_WHEELFILE}
                WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/${TARGET_NAME}"
                COMMENT 
                    "Upload python wheel package to 'testpypi' repository"
)
add_custom_target(whldeploy)
add_dependencies(PyPackageWhlUpload PyPackageBuild)
add_dependencies(whldeploy PyPackageWhlUpload)


# add package test
set(PYTHON_VIRTUAL_ENV "${CMAKE_CURRENT_BINARY_DIR}/pytest_venv")
if (NOT EXISTS "${PYTHON_VIRTUAL_ENV}")
    message(STATUS "Creating virtual python environment for testing [${PYTHON_VIRTUAL_ENV}]")
    execute_process(COMMAND "${Python_EXECUTABLE}" -m venv "${PYTHON_VIRTUAL_ENV}")
else ()
    message(STATUS "Using virtual python environment for testing [${PYTHON_VIRTUAL_ENV}]")
endif ()

# update the environment with VIRTUAL_ENV variable (mimic the activate script) and search again
set(ENV{VIRTUAL_ENV} "${PYTHON_VIRTUAL_ENV}")
set(Python3_FIND_VIRTUALENV ONLY)
unset(Python3_EXECUTABLE)
find_package(Python3 COMPONENTS Interpreter)

message(STATUS "Found python executable for building: ${Python_EXECUTABLE}")
message(STATUS "Found python executable for testing (venv): ${Python3_EXECUTABLE}")

add_custom_target(PyPackageTestInstall
        #DEPENDS PyPackageBuild
        COMMAND "${Python3_EXECUTABLE}" -m pip install "${CMAKE_BINARY_DIR}/${TARGET_NAME}/${PYPACKAGE_NAME}"
        COMMENT "Installing python package in virtual test environment"
        )
add_custom_target(PyPackageTestRun
        #DEPENDS PyPackageBuild
        COMMAND "${Python3_EXECUTABLE}" -c "import pygimli as pg; pg.version()"
        VERBATIM
        COMMENT "Running pygimli installation in virtual test environment."
        )

add_custom_target(whltest)
add_dependencies(PyPackageTestInstall whl)
add_dependencies(PyPackageTestRun PyPackageTestInstall)
add_dependencies(whltest PyPackageTestRun)


#get_target_property(PYPACKAGE_NAME PyPackageBuild Package_Name)

# add_custom_target(PyTests ALL
#         DEPENDS PyPackageBuild
#         COMMAND "${Python3_EXECUTABLE}" -m pip install "${CMAKE_BINARY_DIR}/${TARGET_NAME}/${PYPACKAGE_NAME}"
#         COMMENT "Installing python package in virtual test environment"
#         )

# add_test(NAME package:wheel
#         COMMAND "${Python3_EXECUTABLE}" ${CMAKE_CURRENT_SOURCE_DIR}/test_sample.py
#         )
