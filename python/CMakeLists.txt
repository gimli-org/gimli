
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

add_python_module(pygimli "${CMAKE_CURRENT_BINARY_DIR}/generated/" "gimli" "${CMAKE_CURRENT_SOURCE_DIR}/pygimli/")

add_dependencies(fixpygimlisrc pygimlisrc)
add_dependencies(buildpygimli fixpygimlisrc)
add_dependencies(pygimli buildpygimli)

# for calling the build only
ADD_CUSTOM_TARGET(pg ALL)
add_dependencies(pg _pygimli_)

set(SCRIPTS pytripatch.py; meshmerge.py; meshconvert.py; pygi; datamerge.py)

foreach(program ${SCRIPTS})
    install(PROGRAMS apps/${program} DESTINATION bin)
endforeach(program)

if (WIN32)
	file(COPY ${Boost_PYTHON_LIBRARIES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
	install(PROGRAMS ${Boost_PYTHON_LIBRARIES} DESTINATION bin)
endif()
