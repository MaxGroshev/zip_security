cmake_minimum_required (VERSION 3.21)

#-----------------------------------------------------------------------------------------

project (encryptor_tests)

#-----------------------------------------------------------------------------------------

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 20)

#-----------------------------------------------------------------------------------------

#Unit_tests
set (UNIT_TESTS
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

#-----------------------------------------------------------------------------------------


add_executable            (encryptor_unit_tests ${UNIT_TESTS} ${DEBUG_UTILS})
target_include_directories(encryptor_unit_tests PRIVATE${GTEST_INCLUDE_DIRS}
                                      ${CMAKE_CURRENT_SOURCE_DIR}/../)
target_link_libraries     (encryptor_unit_tests ${GTEST_BOTH_LIBRARIES} debug_utils)
