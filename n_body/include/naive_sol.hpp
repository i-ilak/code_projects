#ifndef NAIVE_SOL_HPP
#define NAIVE_SOL_HPP

#include <chrono>
#include <vector>

#include "stelar_object.hpp"
#include "integrate.hpp"
#include "n_body_solver.hpp"

using planet_container_type = std::vector<StellarObject>;

void perform_simulation(planet_container_type const & ,
                        double const & , int const & , double const & ,
                        int const &);

#endif // NAIVE_SOL_HPP