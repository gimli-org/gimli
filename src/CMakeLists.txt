set( libgimli_TARGET_NAME gimli )

link_directories(${libgimli_BINARY_DIR})

file (GLOB SOURCE_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.cpp)
file (GLOB HEADER_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.h)

if (NOT MSVC)
	list(REMOVE_ITEM SOURCE_FILES wingetopt.cpp)
endif(NOT MSVC)

###########################################################################
#                       SHARED LIBRARY SET UP
###########################################################################
add_library(${libgimli_TARGET_NAME} SHARED ${SOURCE_FILES} ${HEADER_FILES})

include (GenerateExportHeader)

generate_export_header( ${libgimli_TARGET_NAME}
             BASE_NAME ${libgimli_TARGET_NAME}
             EXPORT_MACRO_NAME DLLEXPORT
             EXPORT_FILE_NAME ${libgimli_TARGET_NAME}_export.h
             STATIC_DEFINE LIBGIMLI_BUILT_AS_STATIC
)

if (TRIANGLE_FOUND)
    include_directories(${Triangle_INCLUDE_DIR})
    target_link_libraries(${libgimli_TARGET_NAME} ${Triangle_LIBRARIES})
endif (TRIANGLE_FOUND)

if (CHOLMOD_FOUND)
    include_directories(${CHOLMOD_INCLUDE_DIRS})
    target_link_libraries(${libgimli_TARGET_NAME} ${CHOLMOD_LIBRARIES})
endif (CHOLMOD_FOUND)

if (UMFPACK_FOUND)
    include_directories(${UMFPACK_INCLUDE_DIRS})
    target_link_libraries(${libgimli_TARGET_NAME} ${UMFPACK_LIBRARIES})
endif (UMFPACK_FOUND)

if (Boost_THREAD_FOUND)
    target_link_libraries(${libgimli_TARGET_NAME} ${Boost_THREAD_LIBRARIES})
endif (Boost_THREAD_FOUND)

if (Boost_SYSTEM_FOUND)
    target_link_libraries(${libgimli_TARGET_NAME} ${Boost_SYSTEM_LIBRARIES})
endif (Boost_SYSTEM_FOUND)

if (WIN32)
    target_link_libraries(${libgimli_TARGET_NAME} psapi)

	file(COPY ${Boost_THREAD_LIBRARIES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
	file(COPY ${Boost_SYSTEM_LIBRARIES} DESTINATION ${CMAKE_BINARY_DIR}/bin)

	install(PROGRAMS ${Boost_THREAD_LIBRARIES} DESTINATION bin)
	install(PROGRAMS ${Boost_SYSTEM_LIBRARIES} DESTINATION bin)

	if(${BLAS_LIBRARIES})
		file(COPY ${BLAS_LIBRARIES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
		if (LAPACK_LIBRARIES)
			file(COPY ${LAPACK_LIBRARIES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
		endif()
		install(PROGRAMS ${BLAS_LIBRARIES} DESTINATION bin)
		install(PROGRAMS ${LAPACK_LIBRARIES} DESTINATION bin)
	endif(${BLAS_LIBRARIES})

elseif (UNIX AND NOT APPLE)
	if (READPROC_FOUND)
        target_link_libraries(${libgimli_TARGET_NAME} ${CMAKE_DL_LIBS})
		include_directories(${READPROC_INCLUDE_DIR})
		target_link_libraries(${libgimli_TARGET_NAME} ${READPROC_LIBRARIES})
	endif (READPROC_FOUND)
	install(TARGETS ${libgimli_TARGET_NAME} LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})	
endif(WIN32)

#install(FILES gimli.h ${PROJECT_BINARY_DIR}/${libgimli_TARGET_NAME}_export.h 
#        DESTINATION ${INCLUDE_INSTALL_DIR})


