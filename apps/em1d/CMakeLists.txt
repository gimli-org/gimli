#link_directories (${libgimli_BINARY_DIR}/src) 

set( PROGRAMS dc1dinv dc1dsmooth fdem1dinv fdem1dsmooth mrs1dblock mrs1dsmooth
    mt1dinv mt1dsmooth )

foreach(program ${PROGRAMS})
    add_executable(${program} ${program}.cpp)
    target_link_libraries(${program} gimli)
    #target_link_libraries(polyfit ${libgimli_TARGET_NAME} gimli)
endforeach(program)