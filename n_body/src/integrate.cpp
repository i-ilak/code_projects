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
Eigen::MatrixXd integrate(step_t step, rhs_t rhs, phase_t const & z0,
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


/*
 For all methods we divide the computation into steps. Each step gets perform N times. 
 The following function implement this step. Later, these will be called and summed up.
 */

/*
  PRE:
    - rhs:      f(t,z) for a given ODE.
    - z0:       initial condition in phase space coordinates. 
    - t0:       time at which we start the step
    - dt:       time step, usually defined as time interval of the simulation
                divided by number of total steps taken in the simulation.
  POST:
    - z+dz:     phase space type coordinates
 */
phase_t explicit_euler_step(rhs_t rhs, phase_t const & z0, double t0, double dt){
    return z0 + dt * rhs(t0,z0);
}

phase_t explicit_midpoint_step(rhs_t rhs, phase_t const & z0, double t0, double dt){
    phase_t tmp = z0;
    tmp += (.5 * dt) * rhs(t0, z0);
    return z0 + dt * rhs(t0 + .5*dt, tmp);
}

phase_t velocity_verlet_step(rhs_t rhs, phase_t const & z0, double t0, double dt) {
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
Eigen::MatrixXd explicit_euler(rhs_t rhs, phase_t const &z0, double const & T, unsigned const & N){
    return integrate(explicit_euler_step,rhs,z0,T,N);
}

Eigen::MatrixXd explicit_midpoint(rhs_t rhs, phase_t const &z0, double const & T, unsigned const & N){
    return integrate(explicit_midpoint_step,rhs,z0,T,N);
}

Eigen::MatrixXd velocity_verlet(rhs_t rhs, phase_t const &z0, double const & T, unsigned const & N){
    return integrate(velocity_verlet_step,rhs,z0,T,N);
}