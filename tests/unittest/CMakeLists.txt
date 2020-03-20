file (GLOB unitTest_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.h)

#PURIFYCOMMAND:FILEPATH=c:/Progra~1/Rational/common/purify.exe
#MEMORYCHECK_COMMAND:FILEPATH=/home/kitware/local/bin/valgrind
#set(WITH_MEMCHECK TRUE)
#set(WITH_COVERAGE TRUE)

add_executable(gimliUnitTest unitTest.cpp ${unitTest_HEADERS})

if (CPPUNIT_FOUND)
    add_definitions(-DHAVE_UNITTEST)
    include_directories(${CPPUNIT_INCLUDE_DIR})
    target_link_libraries(gimliUnitTest ${CPPUNIT_LIBRARIES})


    #if (CHOLMOD_FOUND)
    #    target_link_libraries(gimliUnitTest ${CPPUNIT_LIBRARIES})
    #endif (CHOLMOD_FOUND)

    # target_link_libraries(gimliUnitTest ${Boost_THREAD_LIBRARY}) # maybe obsolete 200103
    target_link_libraries(gimliUnitTest gimli)
    target_link_libraries(gimliUnitTest pthread)
    #target_link_libraries(gimliUnitTest ${CHOLMOD_LIBRARIES})

    add_dependencies(check gimliUnitTest)

    ADD_CUSTOM_TARGET(gtest)
    ADD_CUSTOM_COMMAND(
        COMMAND
            ${CMAKE_BINARY_DIR}/bin/gimliUnitTest
        TARGET
            gtest
    )

    add_dependencies(gimliUnitTest gimli)
    add_dependencies(gtest gimliUnitTest)


    add_test(
        NAME
            gtest:gimli_unit_tests
        COMMAND
            #cmake "--build . --target gtest"
            make gtest
        )


    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/unittest.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    #file(COPY ${CMAKE_CURRENT_BINARY_DIR}/../../bin/gimliUnitTest DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    INSTALL_FILES(${CMAKE_CURRENT_SOURCE_DIR} unittest.sh)

endif()
