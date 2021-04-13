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

/* Wrapper function that thakes an initial distribution of planets in objects and performs a 
 * simulation for some time T in N steps with gravitational constant G. 
 * PRE:
 *  - objects:  contains a vector of StellarObject-type
 *  - T:        for which time-span should the simulation run
 *  - N:        how many steps are made
 *  - G:        values of the gravitational constants (specifies unit system)
 *  - k:        integer that specifies integration method
 *                  1.  Explicit Euler
 *                  2.  Explicit Midpoint
 *                  3.  Velocity Verlet
 * 
 * POST:
 *  Creates a .hfd5 file with the data of the simulation
 */
void perform_simulation(planet_container_type const & objects,
                        double const & T, int const & N, 
                        double const & G, int const & k);

// Same as above, we just don't need to specify the method here. 
void perform_ode_int_simulation(planet_container_type const & objects,
                        double const & T, int const & N);

#endif // NAIVE_SOL_HPP