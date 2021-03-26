#include "integrate.hpp"

/*
 * Function that perform the integration forward in time. 
 * PRE:
 *  - step:     specifies the method used to perform the integration
 *  - rhs:      specifies the differential equation, i.e. the ODE has 
 *              the form dz = rhs(t,z).
 *  - z0:       initial conditions
 *  - time_interval:
 *              gives the timer-interval over which we integrate the ODE
 *  - bins:     number of steps we take in the simulation. Higher number implies
 *              higher precision, but also higher computation time. 
 * POST:
 *  - Each row of the resulting matrix contains the phase space coordinates at one
 *    fixed time.  
 */

/*

*/

/*
 Explicit Euler method.     For more details see: Wikipedia https://en.wikipedia.org/wiki/Euler_method
 Explicit Midpoint method.  For more detials see: Wikipedia https://en.wikipedia.org/wiki/Midpoint_method
 PRE:
    - f:    characterizes the ODE. 
    - z0:   initial conditions of the problem in form of phase-space
            coordinates.
    - T:    time interval you want to simulate, i.e. has to be a
            positive real number.
    - N:    number of individual steps you want to take in simulation.
            The bigger N is the more precise the result will be and  
            the longer the computation will take. 
 POST:
    - row k contains the phase space coordinates at time k*T/N.
 */
/*
Eigen::MatrixXd explicit_euler(rhs_t rhs, phase_t const &z0, double const & T, unsigned const & N){
    return integrate(explicit_euler_step,rhs,z0,T,N);
}
*/
/*
Eigen::MatrixXd explicit_midpoint(rhs_t rhs, phase_t const &z0, double const & T, unsigned const & N){
    return integrate(explicit_midpoint_step,rhs,z0,T,N);
}

Eigen::MatrixXd velocity_verlet(rhs_t rhs, phase_t const &z0, double const & T, unsigned const & N){
    return integrate(velocity_verlet_step,rhs,z0,T,N);
}
*/