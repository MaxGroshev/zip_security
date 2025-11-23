cmake_minimum_required(VERSION 3.21)

#-----------------------------------------------------------------------------------------

project(my_lzw)

#-----------------------------------------------------------------------------------------

set(TESTING_ENABLED ON CACHE BOOL [FORCE])

#-----------------------------------------------------------------------------------------

set(CMAKE_BUILD_TYPE "Release")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS  -fsanitize=address,leak)
    message("Debug_version")
endif ()

#-----------------------------------------------------------------------------------------


set(CMAKE_LZW_SRC_ROOT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(compressor/)
add_subdirectory(encryptor/)
add_subdirectory(debug_utils/)

if (TESTING_ENABLED)
    add_subdirectory(tests/)
endif()

# unset(TESTING_ENABLED CACHE)


