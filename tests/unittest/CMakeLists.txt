
file (GLOB unitTest_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.h)

#PURIFYCOMMAND:FILEPATH=c:/Progra~1/Rational/common/purify.exe
#MEMORYCHECK_COMMAND:FILEPATH=/home/kitware/local/bin/valgrind
#set(WITH_MEMCHECK TRUE)
#set(WITH_COVERAGE TRUE)

add_executable(gimliUnitTest unitTest.cpp ${unitTest_HEADERS})

if (CPPUNIT_FOUND)
    add_definitions( -DHAVE_UNITTEST )
	include_directories(${CPPUNIT_INCLUDE_DIR})
    target_link_libraries(gimliUnitTest ${CPPUNIT_LIBRARIES})
endif (CPPUNIT_FOUND)

if (CHOLMOD_FOUND)
    target_link_libraries(gimliUnitTest ${CPPUNIT_LIBRARIES})
endif (CHOLMOD_FOUND)

target_link_libraries(gimliUnitTest ${Boost_THREAD_LIBRARY})
target_link_libraries(gimliUnitTest gimli)
target_link_libraries(gimliUnitTest ${CHOLMOD_LIBRARIES})

add_test(unitTest unittest.sh)

add_dependencies(check gimliUnitTest)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/unittest.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
#file(COPY ${CMAKE_CURRENT_BINARY_DIR}/../../bin/gimliUnitTest DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

INSTALL_FILES(${CMAKE_CURRENT_SOURCE_DIR} unittest.sh)
