# This needs to be in both the top-level and the test CMake file.
enable_testing()

include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/src)

add_executable(genome_default_ctor genome_default_ctor.cpp)
add_executable(copy_ctor copy_ctor.cpp)

target_link_libraries(genome_default_ctor penna)
target_link_libraries(copy_ctor penna)

add_test(NAME default_ctor_test COMMAND genome_default_ctor)
add_test(NAME copy_test COMMAND copy_ctor)