// The following code is pretty much verbatim copied from 
// https://github.com/headmyshoulder/odeint-v2/blob/master/examples/solar_system.cpp
// We use it to check our generated solutions. Some small modifications are done.

/* Boost libs/numeric/odeint/examples/solar_system.cpp
 Copyright 2010-2012 Karsten Ahnert
 Copyright 2011 Mario Mulansky
 Solar system example for Hamiltonian stepper
 Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

// We will wrap the code to make it callable in a similar way as the other ode-solvers.
// We also modify the streaming_observer slightly to be able to write the data to a 
// HDF5 document in the end.

#include <H5Cpp.h>
#include <omp.h>

#include <iostream>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

#include "point_type.hpp"
#include "stelar_object.hpp"

//[ container_type_definition
// we simulate 5 planets and the sun
const size_t n = 6;

typedef point< double , 3 > point_type;
typedef boost::array< point_type , n > container_type;
typedef boost::array< double , n > mass_type;
//]

// added container_type for planets
using planet_container_type = std::vector<StellarObject>;


//[ coordinate_function
const double gravitational_constant = 2.95912208286e-4;

struct solar_system_coor;
//]


//[ momentum_function
struct solar_system_momentum;
//]


//[ some_helpers
point_type center_of_mass( const container_type &x , const mass_type &m );

//[ streaming_observer
struct streaming_observer;


void perform_ode_int_simulation(planet_container_type const & objects,
                        double const & T, int const & N);