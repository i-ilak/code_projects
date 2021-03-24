add_library(penna STATIC genome.cpp animal.cpp population.cpp fishing_pop.cpp)

target_include_directories(penna PUBLIC ./)