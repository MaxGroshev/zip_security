cmake_minimum_required (VERSION 3.21)

#-----------------------------------------------------------------------------------------

project (tests)

#-----------------------------------------------------------------------------------------

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 20)

#-----------------------------------------------------------------------------------------

#Unit_tests
set (UNIT_TESTS
    ./unit_tests/main.cpp)

# End_to_end_tests
set (E2E_TESTS
    ./end_to_end_tests/e2e_test.py)

#-----------------------------------------------------------------------------------------

add_executable            (unit_tests ${UNIT_TESTS} ${DEBUG_UTILS})
target_include_directories(unit_tests PRIVATE${GTEST_INCLUDE_DIRS}
                                      ${CMAKE_LZW_SRC_ROOT_DIRECTORY}/compressor/include/
                                      ${CMAKE_LZW_SRC_ROOT_DIRECTORY}/utils)
target_link_libraries     (unit_tests ${GTEST_BOTH_LIBRARIES} debug_utils)
