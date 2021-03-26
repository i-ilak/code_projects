#ifndef INTEGRATE_HPP
#define INTEGRATE_HPP

#include <string>
#include <eigen3/Eigen/Dense>
#include <utility>
#include <iostream>
#include <cmath>
#include <fstream>

/*
 Define some useful types. All the names are related to physical quantities and
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
template <class S>
using step_t = phase_t (*)(S, phase_t const &, double, double);
using integration_t = Eigen::MatrixXd (*)(rhs_t , phase_t const &, double const &, unsigned const &);

/*
 * The following are two, well-known solution schemes for ODE's. For both the
 * following is true.
 * PRE:
 *  - f:    characterizes the ODE. 
 *  - z0:   initial conditions of the problem in form of phase-space
 *          coordinates.
 *  - T:    time interval you want to simulate, i.e. has to be a
 *          positive real number.
 *  - N:    number of individual steps you want to take in simulation.
 *          The bigger N is the more precise the result will be and  
 *          the longer the computation will take. 
 */

template <typename RHS>
phase_t explicit_euler_step(RHS rhs, phase_t const & z0, double t0, double dt){
    return z0 + dt * rhs(t0,z0);
}
template <typename RHS>
phase_t explicit_midpoint_step(RHS rhs, phase_t const & z0, double t0, double dt){
    phase_t tmp = z0;
    tmp += (.5 * dt) * rhs(t0, z0);
    return z0 + dt * rhs(t0 + .5*dt, tmp);
}
template <typename RHS>
phase_t velocity_verlet_step(RHS rhs, phase_t const & z0, double t0, double dt) {
    int const N = z0.size()/2/3;
    phase_t x0 = z0.head(3*N);
    phase_t v0 = z0.tail(3*N);
    
    phase_t x1 = x0 + dt*v0 + dt*dt/2. * rhs(t0,z0).tail(3*N);
    phase_t tmp(6*N);
    tmp << x1, v0;
    phase_t v1 = v0 + dt/2. * (rhs(t0,z0)+rhs(t0+dt,tmp)).tail(3*N);
    phase_t result(6*N);
    result << x1,v1;
    return result;
}


template<typename RHS,
         typename STEP>
Eigen::MatrixXd integrate(STEP step, RHS rhs, phase_t const & z0,
                          double const & time_interval, unsigned const & bins){
    Eigen::MatrixXd z(bins+1, z0.size());
    z.fill(0);
    double const t0 = 0;
    double const dt = time_interval/bins;
    z.row(0) = z0;
    for(std::size_t i = 0; i<bins; ++i) {
        z.row(i + 1) = step(rhs, z.row(i), t0 + i * dt, dt);
    }
    phase_t t(bins+1);
    for(std::size_t k = 0; k<t.size(); ++k){
        t(k) = (k+1) * dt;
    }
    z.conservativeResize(Eigen::NoChange, z.cols()+1);
    z.col(z.cols()-1) = t;
    return z;
}

template <typename RHS>
Eigen::MatrixXd explicit_euler(RHS rhs, phase_t const &z0, double const & T, unsigned const & N){
    return integrate(explicit_euler_step<RHS>,rhs,z0,T,N);
}

template <typename RHS>
Eigen::MatrixXd explicit_midpoint(RHS rhs, phase_t const & z0, double const & T, unsigned const & N){
    return integrate(explicit_midpoint_step<RHS>,rhs,z0,T,N);
}

template <typename RHS>
Eigen::MatrixXd velocity_verlet(RHS rhs, phase_t const &z0, double const & T, unsigned const & N){
    return integrate(velocity_verlet_step<RHS>,rhs,z0,T,N);
}
#endif // INTEGRATE_HPP