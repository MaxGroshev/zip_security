
cmake_minimum_required(VERSION 3.22.1)

project("zipper")

add_library(${CMAKE_PROJECT_NAME} SHARED
        native-lib.cpp)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lzw/debug_utils)
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/lzw
        ${CMAKE_CURRENT_SOURCE_DIR}/lzw/utils)

target_link_libraries(${CMAKE_PROJECT_NAME}
        debug_utils
        android
        log)