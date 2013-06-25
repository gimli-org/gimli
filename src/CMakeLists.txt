set( GIMLI_TARGET_NAME gimli )

include_directories(
	${LIBGIMLI_SOURCE_DIR} ../src
	${Boost_INCLUDE_DIR}
)

link_directories(${LIBGIMLI_BINARY_DIR})

file (GLOB SOURCE_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.cpp)
file (GLOB HEADER_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.h)

###########################################################################
#                       SHARED LIBRARY SET UP
###########################################################################
add_library(${GIMLI_TARGET_NAME} SHARED ${SOURCE_FILES} ${HEADER_FILES_VAR})

include (GenerateExportHeader)
generate_export_header( ${GIMLI_TARGET_NAME}
             BASE_NAME ${GIMLI_TARGET_NAME}
             EXPORT_MACRO_NAME DLLEXPORT
             EXPORT_FILE_NAME ${GIMLI_TARGET_NAME}_export.h
             STATIC_DEFINE LIBGIMLI_BUILT_AS_STATIC
)

if (WIN32)
    target_link_libraries(${GIMLI_TARGET_NAME} psapi)
elseif (UNIX AND NOT APPLE)
    target_link_libraries(${GIMLI_TARGET_NAME} ${CMAKE_DL_LIBS})
endif()

#install(TARGETS ${GIMLI_TARGET_NAME} DESTINATION ${LIBRARY_INSTALL_DIR})

#install(FILES
#    gimli.h
#    ${PROJECT_BINARY_DIR}/${GIMLI_TARGET_NAME}_export.h DESTINATION ${INCLUDE_INSTALL_DIR}
#)


