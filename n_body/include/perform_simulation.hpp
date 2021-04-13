#ifndef NAIVE_SOL_HPP
#define NAIVE_SOL_HPP

#include <chrono>
#include <vector>

#include "stelar_object.hpp"
#include "integrate.hpp"
#include "n_body_solver.hpp"

#include <H5Cpp.h>
#include <omp.h>

#include <iostream>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

#include "point_type.hpp"
#include "stelar_object.hpp"

// we simulate 5 planets and the sun
const size_t n = 6;
const double gravitational_constant = 2.95912208286e-4;

using planet_container_type = std::vector<StellarObject>;

void perform_simulation(planet_container_type const & ,
                        double const & , int const & , double const & ,
                        int const &);

void perform_ode_int_simulation(planet_container_type const & objects,
                        double const & T, int const & N);

#endif // NAIVE_SOL_HPP