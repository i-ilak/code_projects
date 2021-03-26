ADD_LIBRARY(naive_sol STATIC naive_sol.cpp n_body_solver.cpp)

TARGET_INCLUDE_DIRECTORIES(naive_sol PUBLIC ./)