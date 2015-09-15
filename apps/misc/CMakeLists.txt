#link_directories (${libgimli_BINARY_DIR}/src) 

set( PROGRAMS bms2vtk invlinearmat )


foreach(program ${PROGRAMS})
    add_dependencies(apps ${program})
    add_executable(${program} ${program}.cpp)
    target_link_libraries(${program} gimli)
    install(TARGETS ${program}
            DESTINATION bin)
endforeach(program)