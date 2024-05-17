cmake_minimum_required(VERSION 3.24)

project(lab_lib)

file(GLOB lib_files "*.cpp" )

add_library(lab_lib SHARED ${lib_files})

set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../dll)