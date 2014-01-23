
ADD_CUSTOM_TARGET(pggen ALL)
ADD_CUSTOM_COMMAND(
    COMMAND 
        ${PYTHON_EXECUTABLE}
    ARGS 
        ${CMAKE_CURRENT_SOURCE_DIR}/generate_pygimli_code.py --extra-includes=${Boost_INCLUDE_DIR}
    TARGET
        pggen
    DEPENDS 
        #${CMAKE_CURRENT_SOURCE_DIR}/generate_pygimli_code.py 
        #${CMAKE_CURRENT_SOURCE_DIR}/pygimli.h
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
        make -j4 -C ${CMAKE_CURRENT_BINARY_DIR} pg
    TARGET
        pgbuild
   )

file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/generated/dummy.cpp "")


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

# for calling the build only
ADD_CUSTOM_TARGET(pg ALL)
#add_dependencies(pg pgbuild)
add_dependencies(pg _pygimli_)

set(SCRIPTS pytripatch.py; meshmerge.py; meshconvert.py; pygi; datamerge.py)

foreach(program ${SCRIPTS})
    install(PROGRAMS apps/${program} DESTINATION bin)
endforeach(program)

if (WIN32)
	file(COPY ${Boost_PYTHON_LIBRARIES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
	install(PROGRAMS ${Boost_PYTHON_LIBRARIES} DESTINATION bin)
endif()
