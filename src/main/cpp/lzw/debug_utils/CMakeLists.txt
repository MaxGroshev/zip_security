cmake_minimum_required(VERSION 3.21)

#-----------------------------------------------------------------------------------------

project(debug_utils)

#-----------------------------------------------------------------------------------------

set(SOURCE
    ./error_control/error_control.cpp
    ./debug_utils.hpp)

#-----------------------------------------------------------------------------------------

add_library(${PROJECT_NAME} ${SOURCE})

target_include_directories(${PROJECT_NAME} PUBLIC
                            ../debug_utils/
                            ./logs/
                            ./time_control/
                            ./error_control/
                            )




