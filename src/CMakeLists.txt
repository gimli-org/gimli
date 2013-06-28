set( libgimli_TARGET_NAME gimli )

include_directories(${Boost_INCLUDE_DIR})

link_directories(${libgimli_BINARY_DIR})

file (GLOB SOURCE_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.cpp)
file (GLOB HEADER_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.h)

###########################################################################
#                       SHARED LIBRARY SET UP
###########################################################################
add_library(${libgimli_TARGET_NAME} SHARED ${SOURCE_FILES} ${HEADER_FILES_VAR})

include (GenerateExportHeader)

generate_export_header( ${libgimli_TARGET_NAME}
             BASE_NAME ${libgimli_TARGET_NAME}
             EXPORT_MACRO_NAME DLLEXPORT
             EXPORT_FILE_NAME ${libgimli_TARGET_NAME}_export.h
             STATIC_DEFINE LIBGIMLI_BUILT_AS_STATIC
)

if (WIN32)
    target_link_libraries(${libgimli_TARGET_NAME} psapi)
elseif (UNIX AND NOT APPLE)
    target_link_libraries(${libgimli_TARGET_NAME} ${CMAKE_DL_LIBS})
endif(WIN32)

if (TRIANGLE_FOUND)
    include_directories(${Triangle_INCLUDE_DIR})
    target_link_libraries(${libgimli_TARGET_NAME} ${Triangle_LIBRARIES})
endif (TRIANGLE_FOUND)

if (CHOLMOD_FOUND)
    include_directories(${CHOLMOD_INCLUDE_DIRS})
    target_link_libraries(${libgimli_TARGET_NAME} ${CHOLMOD_LIBRARIES})
endif (CHOLMOD_FOUND)


#install(TARGETS ${libgimli_TARGET_NAME} 
#        LIBRARY 
#        DESTINATION ${LIBRARY_INSTALL_DIR})

install(FILES gimli.h
              ${PROJECT_BINARY_DIR}/${libgimli_TARGET_NAME}_export.h 
        DESTINATION ${INCLUDE_INSTALL_DIR})


