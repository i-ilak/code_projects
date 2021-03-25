// We will adapt here a streamed down version of the following code:
// https://github.com/headmyshoulder/odeint-v2/blob/master/examples/solar_system.cpp

#ifndef ODE_INT_SOL_HPP
#define ODE_INT_SOL_HPP

#include "stelar_object.hpp"

#include <boost/numeric/odeint.hpp>
#include <boost/array.hpp>
#include "point_type.hpp"               // a sophisticated way of definig a point in space

using point_type = point<double, 3>;    
const std::size_t number_of_planets = 6;
using container_type = boost::array<point_type, number_of_planets>;

void perform_ode_int_simulation(std::vector<StellarObject> const & objects,
                                double const & T, int const & N, double const & G){


using stepper_type = boost::numeric::odeint::symplectic_rkn_sb3a_mclachlan<container_type>;

}

#endif // ODE_INT_SOL_HPP