#link_directories (${libgimli_BINARY_DIR}/src) 

set( PROGRAMS ttinv
              ttmod
#ttresolution
#ttinvoffset 
)

foreach(program ${PROGRAMS})
    add_dependencies(apps ${program})
    add_executable(${program} ${program}.cpp)
    target_link_libraries(${program} gimli)
    #target_link_libraries(polyfit ${libgimli_TARGET_NAME} gimli)
endforeach(program)