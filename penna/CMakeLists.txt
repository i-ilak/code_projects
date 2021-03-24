cmake_minimum_required(VERSION 2.8)

project(penna)

set(CMAKE_CXX_FLAGS "-march=native -O3")

add_subdirectory(src)
add_subdirectory(test)

enable_testing()

add_executable(main main.cpp)

target_link_libraries(main penna)
