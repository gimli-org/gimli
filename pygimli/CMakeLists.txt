## PRE Build
cprint(green "********** Prepare wheel build ***********")

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
execute_process(COMMAND "${Python_EXECUTABLE}" getwheelname.py 
                RESULT_VARIABLE STATUS
                OUTPUT_VARIABLE PY_WHEELFILE
                WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
                )
                
if (STATUS EQUAL 1)
    message( FATAL_ERROR "Can't determine wheel filenname. Is wheel installed? python -m pip install wheel")
endif()

message(STATUS "Detected Python wheel filename: ${PY_WHEELFILE}")
string(REPLACE "linux" "manylinux_2_28" PY_WHEELFILE_OUTPLATFORM ${PY_WHEELFILE})
message(STATUS "Desired Python wheel filename after audtitwheel repair: ${PY_WHEELFILE_OUTPLATFORM}")

# add virtual env for whl package test
set(PYTHON_VIRTUAL_ENV "${CMAKE_CURRENT_BINARY_DIR}/wheeltest_venv")
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


## POST Build

## create the wheel: make whl
add_custom_target(whlCopyLibs)
set(WHEEL_LIBRARIES_DIR ${CMAKE_BINARY_DIR}/${PYTHON_MODULE_NAME}/${TARGET_NAME}/${TARGET_NAME}/core/libs/ )

add_custom_command(TARGET whlCopyLibs
    COMMAND ${CMAKE_COMMAND} -E make_directory ${WHEEL_LIBRARIES_DIR}
    COMMENT "Create target dir to copy libraries: ${WHEEL_LIBRARIES_DIR}"
)
foreach(file 
        #$<TARGET_FILE:${Boost_PYTHON_LIBRARY}>
        #$<TARGET_FILE:gimli>
        $<TARGET_FILE:_pygimli_>
    )
    add_custom_command(TARGET whlCopyLibs
            COMMAND ${CMAKE_COMMAND} -E copy_if_different ${file} ${WHEEL_LIBRARIES_DIR}
            DEPENDS ${file}
            COMMENT "Copy libraries: ${file} ${CMAKE_BINARY_DIR}/${TARGET_NAME}/${TARGET_NAME}/core/libs"
    )
endforeach(file)
add_custom_target(whlPackageBuild DEPENDS whlCopyLibs
    COMMAND "${Python_EXECUTABLE}" -m pip wheel --wheel-dir=${CMAKE_CURRENT_BINARY_DIR}/wheelhouse .
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Building python wheel package"
)
add_custom_target(whlPackageRepair DEPENDS whlPackageBuild
    COMMAND "auditwheel" repair ./wheelhouse/${PY_WHEELFILE}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Building python wheel package at: ${CMAKE_CURRENT_BINARY_DIR}"
)
add_custom_target(whl DEPENDS whlPackageRepair)


## test the wheel: make whlTest
add_custom_target(whlTestInstall DEPENDS whl
    COMMAND "${Python3_EXECUTABLE}" -m pip uninstall -y ${TARGET_NAME}
    COMMAND "${Python3_EXECUTABLE}" -m pip install ${TARGET_NAME} --find-links=${CMAKE_CURRENT_BINARY_DIR}/wheelhouse
    COMMENT "Installing ${TARGET_NAME} in virtual test environment"
)
add_custom_target(whlTestRun DEPENDS whlTestInstall
    COMMAND ${Python3_EXECUTABLE} -c "import pygimli as pg; pg.version()"
    VERBATIM
    COMMENT "Running pygimli installation in virtual test environment. ${Python3_EXECUTABLE} -c \'import pygimli as pg; pg.version()\'"
)
add_custom_target(whlTest DEPENDS whlTestRun)


## deploy the wheel: make whlDeploy
add_custom_target(whlDeploy DEPENDS whlTest
    COMMAND "${Python_EXECUTABLE}" -m twine upload --repository testpypi ${CMAKE_CURRENT_BINARY_DIR}/wheelhouse/${PY_WHEELFILE_OUTPLATFORM}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Upload python wheel package ${PY_WHEELFILE_OUTPLATFORM} to 'testpypi' repository"
)
