#link_directories (${libgimli_BINARY_DIR}/src) 

set( PROGRAMS polyfit harmfit )

foreach(program ${PROGRAMS})
    add_executable(${program} ${program}.cpp)
    target_link_libraries(${program} gimli)
    install(TARGETS ${program}
            DESTINATION bin)
endforeach(program)