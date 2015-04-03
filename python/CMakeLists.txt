ADD_CUSTOM_TARGET(pggen ALL)
ADD_CUSTOM_COMMAND(
    COMMAND
        #${CMAKE_COMMAND} -E env "PYTHONPATH=$ENV{PYTHONPATH}${PYGCCXML_PATH}" ${PYTHON_EXECUTABLE} will work on cmake3.1
        ${PYTHON_EXECUTABLE}
    ARGS
        ${CMAKE_CURRENT_SOURCE_DIR}/generate_pygimli_code.py --extra-path=${PYGCCXML_PATH} --extra-includes=${Boost_INCLUDE_DIR} --gccxml=${GCCXML_EXECUTABLE}
    TARGET
        pggen
    DEPENDS
        #set(ENV{PYTHONPATH} $ENV{PYTHONPATH}${PYGCCXML_PATH})
        #${CMAKE_CURRENT_SOURCE_DIR}/generate_pygimli_code.py
        #${CMAKE_CURRENT_SOURCE_DIR}/pygimli.h
        #ENV PYTHONPATH $ENV{PYTHONPATH}${PYGCCXML_PATH} APPEND
   )

#HACK!! but we need this to force cmake looking for our generated binding files
ADD_CUSTOM_TARGET(fixpygimlisrc ALL)
ADD_CUSTOM_COMMAND(
    COMMAND
        ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
    TARGET
        fixpygimlisrc
   )

ADD_CUSTOM_TARGET(pgbuild ALL)
ADD_CUSTOM_COMMAND(
    COMMAND
        test -z "$J" && ${CMAKE_MAKE_PROGRAM} -j 1 -C ${CMAKE_BINARY_DIR} pg || ${CMAKE_MAKE_PROGRAM} -j $J -C ${CMAKE_BINARY_DIR} pg
    TARGET
        pgbuild
   )

find_file (HAVE_DUMMY dummy.cpp PATHS ${CMAKE_CURRENT_BINARY_DIR}/generated/)
if (NOT HAVE_DUMMY)
    file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/generated/dummy.cpp "")
endif()


# if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANGXX)
#     set (EXTRA_OPTIONS "-Wno-uninitialized -Wno-reorder
#                         -Wno-unused-but-set-variable")
#
#     set_source_files_properties(${SOURCE_FILES}
#         PROPERTIES
#         COMPILE_FLAGS  "${EXTRA_OPTIONS}")
# endif()

add_python_module(pygimli "${CMAKE_CURRENT_BINARY_DIR}/generated/" "gimli" "${CMAKE_CURRENT_SOURCE_DIR}/pygimli/")

add_dependencies(pggen gimli)
add_dependencies(fixpygimlisrc pggen)
add_dependencies(pgbuild fixpygimlisrc)
add_dependencies(pygimli pgbuild)

#ADD_CUSTOM_TARGET(pgrebuild ALL)
#add_dependencies(pgrebuild pggen _pygimli_)

# for calling the build only
ADD_CUSTOM_TARGET(pg ALL)
add_dependencies(pg _pygimli_)
include_directories(${PY_NUMPY}/core/include/)

set(SCRIPTS pytripatch.py; meshmerge.py; meshconvert.py; pygi; datamerge.py)

foreach(program ${SCRIPTS})
    install(PROGRAMS apps/${program} DESTINATION bin)
endforeach(program)

if (WIN32)
	file(COPY ${Boost_PYTHON_LIBRARIES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
	install(PROGRAMS ${Boost_PYTHON_LIBRARIES} DESTINATION bin)
endif()
