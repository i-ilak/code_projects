#ifndef INTEGRATE_HPP
#define INTEGRATE_HPP

#include <string>
#include <eigen3/Eigen/Dense>
#include <utility>
#include <iostream>
#include <cmath>
#include <fstream>

/*
 Define some useful types. All the names are related to physicsl quantities and
 should be interpreted in this context. 
 */
using mass_t = double const;
using name_t = std::string const;
using vector_t = Eigen::Vector3d;
using phase_t = Eigen::VectorXd;    // type used for elements of "phase space"

// The following is in relation to how ODE's are defined. If we are given the
// most general form of an ODE, z' = f(t,z), where t is a constant (usually
// called time) and z an element of the phase space, then rhs_t is the type of
// f. 
using rhs_t = phase_t (*)(double const &, phase_t const &);
using step_t = phase_t (*)(rhs_t , phase_t const &, double, double); // redundant an can probably be eliminated 

/*
 The following are two, well-known solution schemes for ODE's. For both the
 following is true.
 PRE:
    rhs_t f:            characterizes the ODE. 
    phase_t const z0:   initial conditions of the problem in form of phase-space
                        coordinates.
    double const T:     time interval you want to simulate, i.e. has to be a
                        positive real number.
    unsigned const N:   number of individual steps you want to take in simulation.
                        The bigger N is the more precise the result will be and  
                        the longer the computation will take. 
 */
Eigen::MatrixXd explicit_euler(rhs_t, phase_t const &, double const &, unsigned const &);
Eigen::MatrixXd explicit_midpoint(rhs_t, phase_t const &, double const &, unsigned const &);

#endif // INTEGRATE_HPP