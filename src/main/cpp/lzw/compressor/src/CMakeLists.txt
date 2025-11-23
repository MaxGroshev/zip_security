set(SRCS
    lzw.cpp
    utils.cpp
)


add_library(src ${SRCS})
target_include_directories(src PRIVATE ${CMAKE_BINARY_DIR})

